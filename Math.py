import numpy as np
from scipy.special import wofz
import astropy.constants as cst

##### Voigt Calculation #####
def voigtMath(x, alpha, gamma):
    """
    Returns a NORMALIZED Voigt profile whose peak area is unity.
    The Voigt profile peaks at zero with Lorentzian component HWHM gamma
    plus Gaussian component HWHM alpha.
    Voigt profile created using scipy.special.wofz, which returns the value
    of a Faddeeva function.
    WARNING
    scipy.special.wofz is not compaible with np.float128 type parameters.

    Args:
        x (float64): Dimensionless point/array, centered at zero
        alpha (float64): Gaussian HWHM component, in AA
        gamma (float64): Lorentzian HWHM component, in AA

    OUTPUT:
    standard_voigt:  [ndarray]       a Voigt profile centered at zero with peak area of unity
    """

    sigma = alpha / np.sqrt(2 * np.log(2))

    standard_voigt = np.real(wofz((x + 1j*gamma)/sigma/np.sqrt(2))) / sigma/np.sqrt(2*np.pi)

    return standard_voigt

def voigtAbsorption(lam, lam_0, b, d, area = 1.0):
    """
    returns a Voigt absorption with "EW" of area, and centered at lam_0;
    b, and d stands for the Gaussian and Lorentzian components in Voigt profile;
    The output can be directly used in a DIB component,
    or treated as optical depth profile and "exponentialed" before output

    Args:
        lam (float64): input wavelength grid, in AA
        lam_0 (float64): central wavelength, in AA
        b (float64): Gaussian standard deviation, in km/s
        d (float64): Lorentzian damping parameter, in AA
        area (float 64): peak area of the profile, in AA

    OUTPUT:
    voigt_absorption:  [ndarray]       voigt profile at lam_0, b, d, with peak area of area.
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
    returns an absorption line whose optical depth has Voigt profile.
    the profile centers at lam_0, with b and d for Gaussian and Lorentzian components.
    the total optical can be given by tau_0, or by N and f for know line

    Args:
        lam (float64): input wavelength grid, in AA
        lam_0 (float64): central wavelength, in AA
        b (float64): Gaussian standard deviation, in km/s
        d (float64): Lorentzian damping parameter, in AA
        tau_0 (float64): TOTAL optical depth of the absorption, overrode if N and f are set
        N (float64): Column density, in cm^-2 (??)
        f (float64): Oscillator strength

    OUTPUTS:
    transmission:   [ndarray]       flux array of light transmission
    """
    # create tau_0
    if (N is not None) and (f is not None):
        NF = N * f
        tau_0 = np.pi * (cst.e.esu.value) ** 2 * Nf * (1e-8 * lam_0) ** 2 / (
                    cst.m_e.to('g').value * (cst.c.to('cm/s').value) ** 2)  # cm
        # Convert cm to angstroms
        tau_0 *= 1e8

    voigt_opticaldepth = voigtAbsorption(lam, lam_0, b, d, tau_0)

    transmission = np.exp(voigt_opticaldepth)

    return transmission

