import numpy as np
from scipy.special import wofz
from scipy.interpolate import CubicSpline
import astropy.constants as cst

from sherpa.models.model import ArithmeticModel
from sherpa.models.parameter import Parameter

from edibles.functions.parameter_converter import param_convert
from edibles.functions.voigtMathematical import voigt_math

from edibles.functions.voigt import voigtAbsorptionLine

__all__ = ('Cont1D', 'Voigt1D', 'AstroVoigt1D', 'VoigtAbsorptionLine')


class VoigtAbsorptionLine(ArithmeticModel):
    """Class for modelling Voigt absorption with astronomical parameters.



    Attributes
    ----------


    """

    def __init__(self, name='voigtabsorptionline'):

        self.lam_0 = Parameter(name, 'lam_0', 5000., frozen=False, min=0.0)
        self.b = Parameter(name, 'b', 3.5, frozen=False, min=1e-12)
        self.d = Parameter(name, 'd', 0.0005, frozen=False, min=0)
        self.N = Parameter(name, 'N', 999, frozen=True, hidden=True, min=0.0)
        self.f = Parameter(name, 'f', 999, frozen=True, hidden=True, min=0.0)
        self.tau_0 = Parameter(name, 'tau_0', 0.1, frozen=False, min=0.0)

        ArithmeticModel.__init__(self, name,
                                 (self.lam_0, self.b, self.d, self.N, self.f, self.tau_0))

    def calc(self, pars, x, *args, **kwargs):
        '''

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
        '''

        lam_0, b, d, N, f, tau_0 = pars

        if N != 999:
            transmission = voigtAbsorptionLine(lam=x, lam_0=lam_0, b=b, d=d, N=N, f=f)
        else:
            transmission = voigtAbsorptionLine(lam=x, lam_0=lam_0, b=b, d=d, tau_0=tau_0)

        return transmission


class LinkedWavelengthVoigtAbsorptionLine(ArithmeticModel):
    """
    Identical to VoigtAbsorptionLine but splits lam_0 into two variables, where lam_0 -> lam_0 * k
    Used for linking lam_0 to another model and changing k
    """
    def __init__(self, name='voigtabsorptionline'):
        # lambda' = lambda * (c + v')/(c + v) = lambda * k
        self.k = Parameter(name, 'k', 1.00005, frozen=False, min=1e-12)
        self.lam_0 = Parameter(name, 'lam_0', 5000., frozen=False, min=0.0)
        self.b = Parameter(name, 'b', 3.5, frozen=False, min=1e-12)
        self.d = Parameter(name, 'd', 0.0005, frozen=False, min=0)
        self.N = Parameter(name, 'N', 999, frozen=True, hidden=True, min=0.0)
        self.f = Parameter(name, 'f', 999, frozen=True, hidden=True, min=0.0)
        self.tau_0 = Parameter(name, 'tau_0', 0.1, frozen=False, min=0.0)

        ArithmeticModel.__init__(self, name,
                                 (self.k, self.lam_0, self.b, self.d, self.N, self.f, self.tau_0))

    def calc(self, pars, x, *args, **kwargs):
        '''

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

        '''

        k, lam_0, b, d, N, f, tau_0 = pars

        if N != 999:
            transmission = voigtAbsorptionLine(lam=x, lam_0=k * lam_0, b=b, d=d, N=N, f=f)
        else:
            transmission = voigtAbsorptionLine(lam=x, lam_0=k * lam_0, b=b, d=d, tau_0=tau_0)

        return transmission


class Sightline:
    '''
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
    clouds : dict
        Dict of all the sources in the sightline
    lines : dict
        Dict of all the lines in the sightline
    current_cloud : str
        The current source of lines
    resolution : float64
        The resolution of the data


    '''
    def __init__(self, star_name, cont):

        self.star_name = star_name
        self.model = cont

        self.clouds = {}
        self.lines = {}
        self.lines['all'] = []
        self.current_cloud = None
        self.resolution = 80000

    def dupCloud(self, old_cloud_name, new_cloud_name, k):
        """
        Used for easily linking the wavelengths of one cloud to another in a new
        cloud, creating each of the lines automatically.

            k is related to the velocity relative to the two clouds. k=0.99995 will decrease the wavelength.
            All other parameters are intialized to be the same as the other cloud.

        INPUT:
            old_cloud_name:     [string]
            new_cloud_name:     [string]
            k:                  [float]
        """
        b = self.clouds[old_cloud_name].b.val
        d = self.clouds[old_cloud_name].d.val
        self.addCloud(new_cloud_name, b, d, k)

        for i, line in enumerate(self.lines[old_cloud_name]):
            new_line = VoigtAbsorptionLine(name=new_cloud_name + '-line ' + str(i))
            new_line.lam_0 = self.init_line.k * line.lam_0
            new_line.b = self.init_line.b
            new_line.d = self.init_line.d
            new_line.tau_0 = line.tau_0.val

            self.model *= new_line
            self.lines[new_cloud_name].append(new_line)
            self.lines['all'].append(new_line)

    def addCloud(self, cloud_name, b, d, k=None):
        '''This is an "initializer line" for the purpose of linking 
            the b and d parameters of subsequent lines - VALUES WILL BE HIDDEN

            DONT use k unless lam_0 is linked to another model. It is redundant
        '''
        if k is None:
            init_line = VoigtAbsorptionLine(name=cloud_name)
        else:
            init_line = LinkedWavelengthVoigtAbsorptionLine(name=cloud_name)
            init_line.k = k

        init_line.lam_0 = 5000  # arbitrary
        init_line.lam_0.frozen = True  # frozen to decrease fit params - arbitrary
        # init_line.lam_0.hidden=True
        init_line.b = b
        # init_line.b.hidden=True
        init_line.d = d
        # init_line.d.hidden=True
        init_line.tau_0 = 0.0           # MUST BE ZERO
        init_line.tau_0.frozen = True   # MUST BE FROZEN
        init_line.tau_0.hidden=True

        self.init_line = init_line
        self.model *= init_line

        self.clouds[cloud_name] = init_line
        self.lines[cloud_name] = []  # REMOVES WHAT IS ALREADY THERE
        self.current_cloud = cloud_name

    def setCloud(self, cloud_name):
        """
        Changes which cloud a line is being added to when addLine is called

        """
        if cloud_name in self.clouds.keys():
            self.current_cloud = cloud_name
            self.init_line = self.clouds[self.current_cloud]
        else:
            raise NameError("Cloud {} was not initialized".format(cloud_name))

    def addLine(self, name, lam_0, tau_0):
        '''Creates an instance of a VoigtAbsorptionLine object

        INPUT:
            name:       [string]
            lam_0       [float]
            tau_0       [float]
        OUTPUT:
            line model  [object instance]
        '''
        line = VoigtAbsorptionLine(name=name)
        line.lam_0 = lam_0

        # Attempting to stop line crossover
        if len(self.lines['all']) != 0:

            val_list=[]
            for old_line in self.lines['all']:
                val_list.append(old_line.lam_0.val)
            sorted_val_list = np.sort(val_list)

            try:   # MAX
                max_we_want = min(sorted_val_list[sorted_val_list > line.lam_0.val])
                for old_line in self.lines['all']:
                    if max_we_want == old_line.lam_0.val:
                        max_line = old_line 
                line.lam_0.max = max_line.lam_0.val - 0.25*max_line.lam_0.val/self.resolution
            except ValueError:
                pass

            try:   # MIN
                min_we_want = max(sorted_val_list[sorted_val_list < line.lam_0.val])
                for old_line in self.lines['all']:
                    if min_we_want == old_line.lam_0.val:
                        min_line = old_line
                line.lam_0.min = min_line.lam_0.val + 0.25*min_line.lam_0.val/self.resolution
            except ValueError:
                pass

        line.b = self.init_line.b
        line.d = self.init_line.d
        line.tau_0 = tau_0

        self.model *= line
        self.lines[self.current_cloud].append(line)
        self.lines['all'].append(line)


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

    def __init__(self, name='Known Velocity Line'):
        self.v_cloud = Parameter(name, 'v_cloud', 0.0, frozen=False, min=-200, max=200)
        self.b = Parameter(name, 'b', 3.5, frozen=False, min=1e-12)
        self.d = Parameter(name, 'd', 0.0005, frozen=False, min=1e-12)
        self.N = Parameter(name, 'N', 999, frozen=False, hidden=False, min=0)
        self.f = Parameter(name, 'f', 999, frozen=True, hidden=False, min=0)
        self.lab_lam_0 = Parameter(name, 'lab_lam_0', 5000, frozen=True)

        ArithmeticModel.__init__(self, name,
                                 (self.v_cloud, self.b, self.d, self.N, self.f, self.lab_lam_0))

    def calc(self, pars, x, *args, **kwargs):
        '''
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
        '''

        v_cloud, b, d, N, f, lab_lam_0 = pars

        lam_0 = lab_lam_0 * (1. + v_cloud / cst.c.to('km/s').value)
        # print(v_cloud)

        tau_0 = None

        transmission = voigtAbsorptionLine(x, lam_0=lam_0, b=b, d=d, tau_0=tau_0, N=N, f=f)

        return transmission


class Cont1D(ArithmeticModel):
    """

    A spline continuum.

    Attributes
    ----------

    y1, y2, ... y7, y8:
        initial y_points to fit
    n_points
        number of segments to break spectrum into.

    """

    def __init__(self, name='Cont_flux'):

        self.y1 = Parameter(name, 'y1', 1.0, frozen=True)
        self.y2 = Parameter(name, 'y2', 1.0, frozen=True)
        self.y3 = Parameter(name, 'y3', 1.0, frozen=True)
        self.y4 = Parameter(name, 'y4', None, frozen=True)
        self.y5 = Parameter(name, 'y5', None, frozen=True)
        self.y6 = Parameter(name, 'y6', None, frozen=True)
        self.y7 = Parameter(name, 'y7', None, frozen=True)
        self.y8 = Parameter(name, 'y8', None, frozen=True)

        ArithmeticModel.__init__(self, name,
                                 (self.y1, self.y2, self.y3, self.y4, self.y5, self.y6,
                                  self.y7, self.y8))

    def calc(self, pars, x, *args, **kwargs):

        '''
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
        '''

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


class Voigt1D(ArithmeticModel):
    """
    Voigt function for modeling absorption.

    Attributes
    ----------

    x:
        [ndarray]  Data grid
    cent:
        [float]    Peak of the Voigt profile
    alpha:
        [float]    Gaussian HWHM component
    gamma:
        [float]    Lorentzian HWHM component
    scaling:
        [float]    Height of the peak

    """

    def __init__(self, name='voigt1d'):
        self.cent = Parameter(name, 'cent', 5000., frozen=True)
        self.alpha = Parameter(name, 'alpha', 0.05, frozen=False, min=0)
        self.gamma = Parameter(name, 'gamma', 0.0005, frozen=False, min=0)
        self.scaling = Parameter(name, 'scaling', 1.0, frozen=True, min=0)

        ArithmeticModel.__init__(self, name, (self.cent, self.alpha, self.gamma, self.scaling))

    def calc(self, pars, x, *args, **kwargs):
        """
        Return the Voigt line shape centered at cent with Lorentzian component HWHM gamma
        and Gaussian component HWHM alpha.

        INPUT
        -----

        x:
            [ndarray]  Data grid
        cent:
            [float]    Peak of the Voigt profile
        alpha:
            [float]    Gaussian HWHM component
        gamma:
            [float]    Lorentzian HWHM component
        scaling:
            [float]    Height of the peak

        OUTPUT
        ------

        y:
            [ndarray]  Flux data for given inputs

        """

        (cent, alpha, gamma, scaling) = pars
        sigma = alpha / np.sqrt(2. * np.log(2.))

        z = (x - cent + 1j * gamma) / (sigma * np.sqrt(2.))
        y = scaling * wofz(z).real / (sigma * np.sqrt(2. * np.pi))

        return -y


class AstroVoigt1D(ArithmeticModel):
    """Voigt function for modelling absorption, with astronomical parameters.

    Attributes
    ----------

    x:
        [ndarray]  Data grid
    cent:
        [float]    Peak of the Voigt profile
    alpha:
        [float]    Gaussian HWHM component
    gamma:
        [float]    Lorentzian HWHM component
    scaling:
        [float]    Height of the peak
    

    """

    def __init__(self, name='astrovoigt1d'):
        self.cent = Parameter(name, 'cent', 5000., frozen=True)
        self.alpha = Parameter(name, 'b_eff', 3.5, frozen=False, min=0)
        self.gamma = Parameter(name, 'Gamma', 6.0e07, frozen=False, min=0)
        self.scaling = Parameter(name, 'scaling', 1.0, frozen=True, min=0)

        ArithmeticModel.__init__(self, name,
                                 (self.cent, self.b_eff, self.Gamma, self.scaling))

    def calc(self, pars, x, *args, **kwargs):
        '''
        INPUT:

        x:
            [ndarray]    Wavelength grid

        pars:
            cent:
                [float]      Central wavelength
            b_eff:
                [float]      Velocity width [km/s]
            Gamma:
                [float]      Lorentzian HWHM component * 4pi
            scaling:
                [float]      Height of the line

        OUTPUT:
        y:
            [ndarray]    Voigt profile
        '''

        # cent, b_eff, Gamma, scaling = pars
        new_params = param_convert(pars)

        cent, alpha, gam, scaling = new_params

        y = voigt_math(x, cent, alpha, gam, scaling)

        return -y
