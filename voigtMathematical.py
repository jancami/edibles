import numpy as np
from scipy.special import wofz


def voigt_math(nu, nu0, alpha, gamma):
    """
    Return the Voigt line shape centered at nu0 with Lorentzian component HWHM gamma
    and Gaussian component HWHM alpha.

    Input:
    nu:     [ndarray]  Data grid
    nu0:    [float]    Peak of the Voigt profile
    alpha:  [float]    Gaussian HWHM component
    gamma:  [float]    Lorentzian HWHM component
    """

    sigma = alpha / np.sqrt(2 * np.log(2))

    return np.real(wofz(((nu - nu0) + 1j*gamma)/sigma/np.sqrt(2))) / sigma/np.sqrt(2*np.pi)
