from sherpa.models.model import ArithmeticModel
from sherpa.models.parameter import Parameter

from edibles.functions.parameter_converter import param_convert
from edibles.functions.voigtMathematical import voigt_math


# __all__ = ('Cont1D', )




def _astrovoigt(pars, x):
    '''
    INPUT:
    x:         [ndarray]    Wavelength grid

    pars:
        cent:      [float]      Central wavelength
        b_eff:     [float]      Velocity width [km/s]
        Gamma:     [float]      Lorentzian HWHM component * 4pi
        scaling: [float]      Height of the line

    OUTPUT:
    y:     [ndarray]    Voigt profile
    '''


    # cent, b_eff, Gamma, scaling = pars
    new_params = param_convert(pars)

    cent, alpha, gam, scaling = new_params

    y = voigt_math(x, cent, alpha, gam, scaling)

    return -y




class AstroVoigt1D(ArithmeticModel):
    """A one-dimensional continuum spline.

    The model parameters are:

    x:         [ndarray]  Data grid
    cent:      [float]    Peak of the Voigt profile
    alpha:     [float]    Gaussian HWHM component
    gamma:     [float]    Lorentzian HWHM component
    scaling:   [float]    Height of the peak
    

    """

    def __init__(self, name='astrovoigt1d'):

        # self.flux    = Parameter(name, 'flux', None, alwaysfrozen=True, hidden=True)
        # self.y_pts   = Parameter(name, 'y_pts', None)

        self.cent = Parameter(name, 'cent', 5000., frozen=True)
        self.alpha = Parameter(name, 'b_eff', 3.5, frozen=False, min=0)
        self.gamma = Parameter(name, 'Gamma', 6.0e07, frozen=False, min=0)
        self.scaling = Parameter(name, 'scaling', 1.0, frozen=True, min=0)


        ArithmeticModel.__init__(self, name,
            (self.cent, self.b_eff, self.Gamma, self.scaling))

    def calc(self, pars, x, *args, **kwargs):
        """Evaluate the model"""

        # If given an integrated data set, use the center of the bin
        

        return _astrovoigt(pars, x)