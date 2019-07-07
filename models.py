import numpy as np
from scipy.special import wofz
from scipy.interpolate import CubicSpline


from sherpa.models.model import ArithmeticModel
from sherpa.models.parameter import Parameter

from edibles.functions.parameter_converter import param_convert
from edibles.functions.voigtMathematical import voigt_math

from edibles.voigt import voigtAbsorptionLine

__all__ = ('Cont1D', 'Voigt1D', 'AstroVoigt1D', 'VoigtAbsorptionLine')



class VoigtAbsorptionLine(ArithmeticModel):
    """Voigt function for modelling absorption, with astronomical parameters.

    Attributes
    ----------

    lam:
        [float64]  (Angstroms)  Wavelength grid

    pars:
        lam_0:
            [float64]  (Angstroms)  Central wavelength
        b:
            [float64]  (km/s)       Gaussian standard deviation
        d:
            [float64]  (units)      Damping parameter

        Choose:
            N:        [float64]  (units)      Column density
            f:        [float64]  (units)      Oscillator strength
            ========================  OR  ========================
            tau_0:    [float64]  (units)      Scaling parameter, default = 0.1


    """

    def __init__(self, name='voigtabsorptionline'):

        self.lam_0 = Parameter(name, 'lam_0', 5000., frozen=False)
        self.b = Parameter(name, 'b', 3.5, frozen=False, min=1e-12)
        self.d = Parameter(name, 'd', 0.0005, frozen=False, min=1e-12)
        self.N = Parameter(name, 'N', 999, frozen=True, hidden=True, min=0)
        self.f = Parameter(name, 'f', 999, frozen=True, hidden=True, min=0)
        self.tau_0 = Parameter(name, 'tau_0', 0.1, frozen=False, min=1e-12)

        ArithmeticModel.__init__(self, name, 
                                (self.lam_0, self.b, self.d, self.N, self.f, self.tau_0))

    def calc(self, pars, x, *args, **kwargs):
        '''
        INPUT:

        lam:
            [float64]  (Angstroms)  Wavelength grid

        pars:
            lam_0:
                [float64]  (Angstroms)  Central wavelength
            b:
                [float64]  (km/s)       Gaussian standard deviation
            d:
                [float64]  (units)      Damping parameter

            Choose:
                N:        [float64]  (units)      Column density
                f:        [float64]  (units)      Oscillator strength
                ========================  OR  ========================
                tau_0:    [float64]  (units)      Optical Depth at peak, default = 0.1

        OUTPUT:
        line:
            [ndarray]    Voigt profile
        '''

        lam_0, b, d, N, f, tau_0 = pars


        if N != 999:
            transmission = voigtAbsorptionLine(x, lam_0, b, d, N, f)
        else:
            transmission = voigtAbsorptionLine(x, lam_0, b, d, tau_0)


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

    def __init__(self, name='continuum_flux_value'):


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
        x_sections = np.array_split(x, n_piece*2)

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



        z = (x - cent + 1j*gamma)/ (sigma*np.sqrt(2.))
        y = scaling * wofz(z).real / (sigma*np.sqrt(2.*np.pi))


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


