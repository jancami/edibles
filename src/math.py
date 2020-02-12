import numpy as np
from scipy.special import wofz
import astropy.constants as cst

##### Voigt Calculation #####
def voigtMath(x, alpha, gamma):
    """
    Returns a NORMALIZED Voigt profile whose peak area is unity.
    The Voigt profile peaks at zero with Lorentzian component HWHM gamma
    and Gaussian component HWHM alpha.
    Voigt profile created using scipy.special.wofz, which returns the value
    of a Faddeeva function.
    WARNING
    scipy.special.wofz is not compaible with np.float128 type parameters.
    Units for x, alpha, and gamma not specified but should be the same

    :param x: wavelength grid centered at zero
    :type x: ndarray of float64
    :param alpha: Gaussian HWHM component, in AA
    :type alpha: float
    :param gamma: Lorentzian HWHM component, in AA
    :type gamma: float

    :return: standard_voigt, y points of a Voigt profile as defined by inputs
    :rtype: ndarray of float
    """
    sigma = alpha / np.sqrt(2 * np.log(2))
    standard_voigt = np.real(wofz((x + 1j*gamma)/sigma/np.sqrt(2))) / sigma/np.sqrt(2*np.pi)

    return standard_voigt


def voigtAbsorption(lam, lam_0, b, d, area = 1.0):
    """
    Returns a Voigt absorption defined by lam_0, b, d, and peak area of "area"
    Can be in flux (e.g. for DIBs) or optical depth unit (e.g. for known ISM line).
    Need to be "exponentialed" if in optical depth unit.
    :param lam: Wavelength grid, in AA
    :type lam: ndarray of float
    :param lam_0: central wavelength, in AA
    :type lam_0: float
    :param b: Gaussian standard deviation, in km/s
    :type b: float
    :param d: Lorentzian damping parameter, in AA
    :type d: float
    :param area: peak area of the profile, in AA, default to 1.0
    :type area: float

    :return: voigt_absorption, y points of voigt absorption as defined by lam_0, b, d, and area
    :rtype: ndarray
    """

    # convert lam & lam_0 to x
    x = lam - lam_0

    # convert b to sigma, then alpha
    sigma = b * lam_0 / cst.c.to('km/s').value
    alpha = sigma * np.sqrt(2. * np.log(2.))

    # convert d to gamma -- [ depends on what units we want to use ]
    # Currently, we are using the Lorentzian HWHM. This can easily be changed...
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gamma = d
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # create y data from voigtMath
    y = voigtMath(x, alpha, gamma)

    voigt_absorption = - y * area
    return voigt_absorption


    # Calculate tau_0
    #tau_0 = np.pi * (cst.e.esu.value) ** 2 * Nf * (1e-8 * lam_0) ** 2 / (
     #           cst.m_e.to('g').value * (cst.c.to('cm/s').value) ** 2)  # cm
    # Convert cm to angstroms
    #tau_0 *= 1e8


def voigtOpticalDepthAbsorption(lam, lam_0, b, d, tau_0=0.1, N=None, f=None):
    """
    Returens an absorption whose OPTICAL DEPTH has Voigt profile, as defined by lam_0, b, d.
    Total Optical depth defined by tau_0 or N the column density and f the oscillator strength.
    N and f, if set, will override tau_0.
    The calculation is tau-based, which is more straight-forward

    :param lam: wavelength grid, in AA
    :type lam: ndarray of float
    :param lam_0: central wavelength, in AA
    :type lam_0: float
    :param b: Gaussian standard deviation, in km/s
    :type b: float
    :param d: Lorentzian damping parameter, in AA
    :type d: float
    :param tau_0: TOTAL optical depth of the absorption, overrode by N and f, default to 1.0
    :type tau_0: float
    :param N: column density, in cm^-2
    :type N: float
    :param f: oscillator strength
    :type f: float

    :return: transmission, ypoints of absorption whose optical depth is Vogit profile defined by inputs
    :rtype: ndarray
    """

    # create tau_0
    if (N is not None) and (f is not None):
        Nf = N * f
        tau_0 = np.pi * (cst.e.esu.value) ** 2 * Nf * (1e-8 * lam_0) ** 2 / (
                    cst.m_e.to('g').value * (cst.c.to('cm/s').value) ** 2)  # cm
        # Convert cm to angstroms
        tau_0 *= 1e8

    voigt_opticaldepth = voigtAbsorption(lam, lam_0, b, d, tau_0)

    transmission = np.exp(voigt_opticaldepth)

    return transmission


if __name__ == '__main__':
    print("Math modular for EDIBLES spectrum")
