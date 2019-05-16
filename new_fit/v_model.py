import numpy as np
from scipy.special import wofz

from sherpa.models.model import ArithmeticModel
from sherpa.models.parameter import Parameter

__all__ = ('Cont1D', )




def _voigt(pars, x):
    """
    Return the Voigt line shape centered at cent with Lorentzian component HWHM gamma
    and Gaussian component HWHM alpha.

    Input:  [Angstroms]

    x:         [ndarray]  Data grid
    cent:      [float]    Peak of the Voigt profile
    alpha:     [float]    Gaussian HWHM component
    gamma:     [float]    Lorentzian HWHM component
    scaling:   [float]    Height of the peak

    OUTPUT:

    y:         [ndarray]  Flux data for given inputs

    """

    (cent, alpha, gamma, scaling) = pars
    sigma = alpha / np.sqrt(2. * np.log(2.))



    z = (x - cent + 1j*gamma)/ (sigma*np.sqrt(2.))
    y = scaling * wofz(z).real / (sigma*np.sqrt(2.*np.pi))


    return -y




class Voigt1D(ArithmeticModel):
    """A one-dimensional continuum spline.

    The model parameters are:

    x:         [ndarray]  Data grid
    cent:      [float]    Peak of the Voigt profile
    alpha:     [float]    Gaussian HWHM component
    gamma:     [float]    Lorentzian HWHM component
    scaling:   [float]    Height of the peak
    

    """

    def __init__(self, name='voigt1d'):

        # self.flux    = Parameter(name, 'flux', None, alwaysfrozen=True, hidden=True)
        # self.y_pts   = Parameter(name, 'y_pts', None)

        self.cent = Parameter(name, 'cent', 5000., frozen=True)
        self.alpha = Parameter(name, 'alpha', 0.05, frozen=False, min=0)
        self.gamma = Parameter(name, 'gamma', 0.0005, frozen=False, min=0)
        self.scaling = Parameter(name, 'scaling', 1.0, frozen=True, min=0)

        # self.y1 = Parameter(name, 'y1', 1.0, frozen=True)
        # self.y2 = Parameter(name, 'y2', 1.0, frozen=True)
        # self.y3 = Parameter(name, 'y3', 1.0, frozen=True)
        # self.y4 = Parameter(name, 'y4', None, frozen=True)
        # self.y5 = Parameter(name, 'y5', None, frozen=True)
        # self.y6 = Parameter(name, 'y6', None, frozen=True)
        # self.y7 = Parameter(name, 'y7', None, frozen=True)
        # self.y8 = Parameter(name, 'y8', None, frozen=True)

        # self.n_piece = Parameter(name, 'n_piece', 2, min=0, hard_min=0, frozen=True)

        ArithmeticModel.__init__(self, name,
            (self.cent, self.alpha, self.gamma, self.scaling))

    def calc(self, pars, x, *args, **kwargs):
        """Evaluate the model"""

        # If given an integrated data set, use the center of the bin
        

        return _voigt(pars, x)