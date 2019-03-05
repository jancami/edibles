import numpy as np
from scipy.special import wofz


def voigt_math(x, cent, a1, g1):
    """
    Return the Voigt line shape centered at cent with Lorentzian component HWHM g1
    and Gaussian component HWHM a1.

    Input:  [Angstroms]

    x:      [ndarray]  Data grid
    cent:   [float]    Peak of the Voigt profile
    a1:     [float]    Gaussian HWHM component
    g1:     [float]    Lorentzian HWHM component
    """

    sigma = a1 / np.sqrt(2 * np.log(2))

    return np.real(wofz(((x - cent) + 1j*g1)/sigma/np.sqrt(2))) / sigma/np.sqrt(2*np.pi)
