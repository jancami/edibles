import numpy as np
from scipy.special import wofz


def voigt_math(x, cent, alpha, gamma):
    """
    Return the Voigt line shape centered at cent with Lorentzian component HWHM gamma
    and Gaussian component HWHM alpha.

    Input:
    x:      [ndarray]  Data grid
    cent:   [float]    Peak of the Voigt profile
    alpha:  [float]    Gaussian HWHM component
    gamma:  [float]    Lorentzian HWHM component
    """

    sigma = alpha / np.sqrt(2 * np.log(2))

    return np.real(wofz(((x - cent) + 1j*gamma)/sigma/np.sqrt(2))) / sigma/np.sqrt(2*np.pi)
