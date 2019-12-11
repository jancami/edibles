import numpy as np
from scipy.interpolate import CubicSpline
import astropy.constants as cst
from sherpa.models.model import ArithmeticModel
from sherpa.models.parameter import Parameter

from edibles.edibles.functions.voigt import voigtAbsorptionLine

__all__ = (
    "Cont1D",
    "VoigtAbsorptionLine",
    "LinkedWavelengthVoigtAbsorptionLine",
    "Sightline",
    "KnownVelocityLine",
)


class Cont1D(ArithmeticModel):
    """
    Create a spline through a set number of anchor points.

    Attributes
    ----------
    y1, y2, ... y8 : float64
        y value of anchor points 
    n_points : int
        number of anchor points to put spline.

    """

    def __init__(self, name="Cont_flux"):
        """
        Parameters
        ----------
        y1, y2, ... y7, y8 : float64
            initial y_points to fit
        n_points : int
            number of anchor points to put spline.

        """

        self.y1 = Parameter(name, "y1", 1.0, frozen=True)
        self.y2 = Parameter(name, "y2", 1.0, frozen=True)
        self.y3 = Parameter(name, "y3", 1.0, frozen=True)
        self.y4 = Parameter(name, "y4", None, frozen=True)
        self.y5 = Parameter(name, "y5", None, frozen=True)
        self.y6 = Parameter(name, "y6", None, frozen=True)
        self.y7 = Parameter(name, "y7", None, frozen=True)
        self.y8 = Parameter(name, "y8", None, frozen=True)

        ArithmeticModel.__init__(
            self,
            name,
            (self.y1, self.y2, self.y3, self.y4, self.y5, self.y6, self.y7, self.y8),
        )

    def calc(self, pars, x, *args, **kwargs):

        """
        This function fits a continuum to data separated into n sections
        where the x and y-values are the median of each section using a cubic spline

        INPUT
        -----
        x:
            [ndarray]               wavelength grid (angstroms)
        pars:
            y1 - y8:
                [floats]            input y_points to fit spline

        OUTPUT:
        -------
        y_spline:
            [ndarray]               continuum flux value array
        """

        n_points = 0
        for i in range(len(pars)):

            if np.isnan(pars[i]):
                pars[i] = None

            if pars[i] is not None:
                n_points += 1

        n_piece = n_points - 1

        # split x & y arrays into n_piece*2 sections
        x_sections = np.array_split(x, n_piece * 2)

        # initialize list of points to spline fit
        x_points = [np.min(x)]

        # loop through every other section (1, 3, 5...)
        # make n_piece+1 points to fit a spline through
        # create a spline point on edge of each piece
        for i in range(1, len(x_sections), 2):
            # set x_point
            x_point = np.max(x_sections[i])

            if i == range(len(x_sections))[-1]:
                x_point = np.max(x)

            x_points.append(x_point)

        y_points = []
        for i in range(len(pars)):
            if pars[i] is not None:
                y_points.append(pars[i])
            else:
                break

        spline = CubicSpline(x_points, y_points)
        y_spline = spline(x)

        return y_spline


class VoigtAbsorptionLine(ArithmeticModel):
    """
    Class for modelling Voigt absorption with astronomical parameters.



    """

    def __init__(self, name="voigtabsorptionline"):

        self.lam_0 = Parameter(name, "lam_0", 5000.0, frozen=False, min=0.0)
        self.b = Parameter(name, "b", 3.5, frozen=False, min=1e-12)
        self.d = Parameter(name, "d", 0.0005, frozen=False, min=0)
        self.N = Parameter(name, "N", 999, frozen=True, hidden=True, min=0.0)
        self.f = Parameter(name, "f", 999, frozen=True, hidden=True, min=0.0)
        self.tau_0 = Parameter(name, "tau_0", 0.1, frozen=False, min=0.0)

        ArithmeticModel.__init__(
            self, name, (self.lam_0, self.b, self.d, self.N, self.f, self.tau_0)
        )

    def calc(self, pars, x, *args, **kwargs):
        """

        Parameters
        ----------
        lam : float64
            Wavelength grid
        lam_0 : float64
            Central wavelength
        b : float64
            Gaussian standard deviation
        d : float64
            Damping parameter
        N : float64
            Column density
        f : float64
            Oscillator strength
        tau_0 : float64
            Scaling parameter

        Returns
        -------
        ndarray
            Voigt absorption line transmission
        """

        lam_0, b, d, N, f, tau_0 = pars

        if N != 999:
            transmission = voigtAbsorptionLine(lam=x, lam_0=lam_0, b=b, d=d, N=N, f=f)
        else:
            transmission = voigtAbsorptionLine(
                lam=x, lam_0=lam_0, b=b, d=d, tau_0=tau_0
            )

        return transmission


class LinkedWavelengthVoigtAbsorptionLine(ArithmeticModel):
    """
    Identical to VoigtAbsorptionLine but splits lam_0 into two variables, where lam_0 -> lam_0 * k
    Used for linking lam_0 to another model and changing k
    """

    def __init__(self, name="voigtabsorptionline"):
        # lambda' = lambda * (c + v')/(c + v) = lambda * k
        self.k = Parameter(name, "k", 1.00005, frozen=False, min=1e-12)
        self.lam_0 = Parameter(name, "lam_0", 5000.0, frozen=False, min=0.0)
        self.b = Parameter(name, "b", 3.5, frozen=False, min=1e-12)
        self.d = Parameter(name, "d", 0.0005, frozen=False, min=0)
        self.N = Parameter(name, "N", 999, frozen=True, hidden=True, min=0.0)
        self.f = Parameter(name, "f", 999, frozen=True, hidden=True, min=0.0)
        self.tau_0 = Parameter(name, "tau_0", 0.1, frozen=False, min=0.0)

        ArithmeticModel.__init__(
            self, name, (self.k, self.lam_0, self.b, self.d, self.N, self.f, self.tau_0)
        )

    def calc(self, pars, x, *args, **kwargs):
        """

        Parameters
        ----------
        lam : float64
            Wavelength grid
        lam_0 : float64
            Central wavelength
        b : float64
            Gaussian standard deviation
        d : float64
            Damping parameter
        N : float64
            Column density
        f : float64
            Oscillator strength
        tau_0 : float64
            Scaling parameter

        Returns
        -------
        ndarray
            Voigt profile line absorption

        """

        k, lam_0, b, d, N, f, tau_0 = pars

        if N != 999:
            transmission = voigtAbsorptionLine(
                lam=x, lam_0=k * lam_0, b=b, d=d, N=N, f=f
            )
        else:
            transmission = voigtAbsorptionLine(
                lam=x, lam_0=k * lam_0, b=b, d=d, tau_0=tau_0
            )

        return transmission


class Sightline:
    """
    An object containing continuum and all Voigt models.


    A sightline with multiple sources of lines that share the same b 
    and d parameters. Groups can be Stellar, Interstellar, Telluric, 
    or subsets of each.

    Attributes
    ----------
    star_name : str
        name of the target star
    model : Object
        Combination of continuum and 0 or more Voigt lines
    sources : dict
        Dict of all the sources in the sightline
    lines : dict
        Dict of all the lines in the sightline
    current_source : str
        The current source of lines
    resolution : float64
        The resolution of the data


    """

    def __init__(self, star_name, cont):
        """
        Parameters
        ----------
        star_name : str
            Name of the target star
        cont : Object
            Initial Continuum model of the spectrum

        """

        self.star_name = star_name
        self.model = cont

        self.sources = {}
        self.lines = {}
        self.lines["all"] = []
        self.current_source = None
        self.resolution = 80000

    def dupSource(self, old_source_name, new_source_name, k):
        """
        Used for easily linking the wavelengths of one source to another in a new
        source, creating each of the lines automatically.

        k is related to the velocity relative to the two sources. k=0.99995 will
        decrease the wavelength. All other parameters are intialized to be the 
        same as the other source.


        Parameters
        ----------
        old_source_name : str
            Name of source to be duplicated
        new_source_name : str
            Name of new source
        k : float64
            relative wavelength shift?


        """
        b = self.sources[old_source_name].b.val
        d = self.sources[old_source_name].d.val
        self.addSource(new_source_name, b, d, k)

        for i, line in enumerate(self.lines[old_source_name]):
            new_line = VoigtAbsorptionLine(name=new_source_name + "-line " + str(i))
            new_line.lam_0 = self.init_line.k * line.lam_0
            new_line.b = self.init_line.b
            new_line.d = self.init_line.d
            new_line.tau_0 = line.tau_0.val

            self.model *= new_line
            self.lines[new_source_name].append(new_line)
            self.lines["all"].append(new_line)

    def addSource(self, source_name, b, d, k=None):
        """
        A method to add a source of absorption lines to the sightline object.

        All lines from this source share the same b and d parameters. 

        ### ???
        Do not use k unless lam_0 is linked to another model. It is redundant



        Parameters
        ----------
        source_name : str
            Name of source 
        b : float
            Gaussian standard deviation
        d : float
            Damping parameter (lifetime broadening?)
        k : float64
            Relative wavelength shift

        """
        if k is None:
            init_line = VoigtAbsorptionLine(name=source_name)
            # This is an "initializer line" for the purpose of linking
            # the b and d parameters of subsequent lines - VALUES WILL BE HIDDEN

        else:
            init_line = LinkedWavelengthVoigtAbsorptionLine(name=source_name)
            init_line.k = k

        init_line.lam_0 = 5000  # arbitrary
        init_line.lam_0.frozen = True  # frozen to decrease fit params - arbitrary
        # init_line.lam_0.hidden=True
        init_line.b = b
        # init_line.b.hidden=True
        init_line.d = d
        # init_line.d.hidden=True
        init_line.tau_0 = 0.0  # MUST BE ZERO
        init_line.tau_0.frozen = True  # MUST BE FROZEN
        init_line.tau_0.hidden = True

        self.init_line = init_line
        self.model *= init_line

        self.sources[source_name] = init_line
        self.lines[source_name] = []  # REMOVES WHAT IS ALREADY THERE
        self.current_source = source_name

    def setSource(self, source_name):
        """
        Changes which source a line is being added to when addLine is called

        Parameters
        ----------
        source_name : str
            Name of the source to switch to.
            Source must already be initialized!

        """
        if source_name in self.sources.keys():
            self.current_source = source_name
            self.init_line = self.sources[self.current_source]
        else:
            raise NameError("Source {} was not initialized".format(source_name))

    def addLine(self, name, lam_0, tau_0):
        """
        Creates an instance of a VoigtAbsorptionLine object

        Parameters
        ----------
        name : str
            Name of line
        lam_0 : float64
            Central wavelength of line
        tau_0 : float64
            Oprical depth of center of line

        """
        line = VoigtAbsorptionLine(name=name)
        line.lam_0 = lam_0

        # Attempting to stop line crossover
        if len(self.lines["all"]) != 0:

            val_list = []
            for old_line in self.lines["all"]:
                val_list.append(old_line.lam_0.val)
            sorted_val_list = np.sort(val_list)

            try:  # MAX
                max_we_want = min(sorted_val_list[sorted_val_list > line.lam_0.val])
                for old_line in self.lines["all"]:
                    if max_we_want == old_line.lam_0.val:
                        max_line = old_line
                line.lam_0.max = (
                    max_line.lam_0.val - 0.25 * max_line.lam_0.val / self.resolution
                )
            except ValueError:
                pass

            try:  # MIN
                min_we_want = max(sorted_val_list[sorted_val_list < line.lam_0.val])
                for old_line in self.lines["all"]:
                    if min_we_want == old_line.lam_0.val:
                        min_line = old_line
                line.lam_0.min = (
                    min_line.lam_0.val + 0.25 * min_line.lam_0.val / self.resolution
                )
            except ValueError:
                pass

        line.b = self.init_line.b
        line.d = self.init_line.d
        line.tau_0 = tau_0

        self.model *= line
        self.lines[self.current_source].append(line)
        self.lines["all"].append(line)


class KnownVelocityLine(ArithmeticModel):
    """Voigt function for modelling absorption, with KNOWN astronomical parameters.

    Attributes
    ----------
    pars:
        v_cloud:
            [float64]  (km/s)       Velocity of cloud
        b:
            [float64]  (km/s)       Gaussian standard deviation
        d:
            [float64]  (units)      Damping parameter
        N:  
            [float64]  (1/cm^2)     Column density
        f:
            [float64]  (unitless)  Oscillator strength
        lab_lam_0:
            [float64]  (Angstroms)  Lab rest wavelength (AIR)
    """

    def __init__(self, name="Known Velocity Line"):
        self.v_cloud = Parameter(name, "v_cloud", 0.0, frozen=False, min=-200, max=200)
        self.b = Parameter(name, "b", 3.5, frozen=False, min=1e-12)
        self.d = Parameter(name, "d", 0.0005, frozen=False, min=1e-12)
        self.N = Parameter(name, "N", 999, frozen=False, hidden=False, min=0)
        self.f = Parameter(name, "f", 999, frozen=True, hidden=False, min=0)
        self.lab_lam_0 = Parameter(name, "lab_lam_0", 5000, frozen=True)

        ArithmeticModel.__init__(
            self, name, (self.v_cloud, self.b, self.d, self.N, self.f, self.lab_lam_0)
        )

    def calc(self, pars, x, *args, **kwargs):
        """
        INPUT:

        x:
            [float64]  (Angstroms)  Wavelength grid
        pars:
            v_cloud:
                [float64]  (km/s)       Velocity of cloud
            b:
                [float64]  (km/s)       Gaussian standard deviation
            d:
                [float64]  (units)      Damping parameter
            N:  
                [float64]  (1/cm^2)     Column density
            f:
                [float64]  (unitless)  Oscillator strength
            lab_lam_0:
                [float64]  (Angstroms)  Lab rest wavelength (AIR)

        OUTPUT:

        line:
            [ndarray]    Voigt profile
        """

        v_cloud, b, d, N, f, lab_lam_0 = pars

        lam_0 = lab_lam_0 * (1.0 + v_cloud / cst.c.to("km/s").value)
        # print(v_cloud)

        tau_0 = None

        transmission = voigtAbsorptionLine(
            x, lam_0=lam_0, b=b, d=d, tau_0=tau_0, N=N, f=f
        )

        return transmission
