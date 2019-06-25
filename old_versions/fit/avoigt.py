# +
# NAME:
#     avoigt
#
# PURPOSE:
#     defining a voigt profile for using in fitting and modeling
#     astronomical features.
#
#
# INPUTS:
#
#   x         : Wavelength grid in Angstroms
#   v_cloud   : The cloud component velocity in km/s
#   f         : Oscillator strength
#   logN      : log of Column density in units of cm^-2
#   b_eff     : Velocity width of the Voigt profile in km/s
#               the output b_eff should be and is in km/s but since we
#               convert all scales to cm, we should convert b to cm/s
#   gam       : Radiation damping constant, or Einstein constant (A_ul)
#   z         : The redshift of the observed wavelength grid
#
#
#
# OUTPUT:
#     voigt_model : return a voigt model based on the input parameters
#
# AUTHOR:
#     Amin Farhang
#     University of Western Ontario
#
# DATE:
#     V1.1
#     Aug 15 2018
#     V1.2
#     Oct 1 2018
#
# Copyright:
#     The Copyright for these codes belongs to EDIBLES survey.
#     For use it please refer to survey P.I. Nick Cox, Jan Cami or
#     send an email to Amin Farhang.
#     farhang.amin@gmail.com
#
# +
from __future__ import print_function
import sys
import numpy as np
from scipy.special import wofz
from scipy.signal import fftconvolve, gaussian
from astropy import constants as cst


# =========================================================================
#   Voigt model -- the astronomical version of Voigt function
# =========================================================================
def voigt(x, lambda_peak=None, b_eff=None, log_N=None, gamma=None, osc_freq=None, resolving_power=None):

    # check negative wavelength
    if any(nn < 0 for nn in x):
        print('The wavelength contain negative values !!!')
        sys.exit()


    # check existense of lambda_peak
    if lambda_peak is None:
        print('The lambda_peak is not defined !!!')
        sys.exit()
    else: lambda_peak = np.float(lambda_peak)


    # check existense of b_eff
    if b_eff is None:
        print('The b_eff is not defined !!!')
        sys.exit()

    # check the oscillator frequencty
    if osc_freq is None:
        osc_freq = 1.0

    # check the gamma
    if gamma is None:
        gamma = 1.0e06

    # check the resolving power (R)
    if resolving_power is None:
        resolving_power = 80000


    # --------------------
    # define the constants
    # --------------------
    # instrumental resolution defined with c/R and determine the velocity resolution in (km/s)
    resolution = cst.c.to('km/s').value / np.float(resolving_power)

    sigma0 = 0.02654
    # resolving power of 0.1 km/s
    R = cst.c.to('cm/s').value * 1.0e-5 / 0.1
    delta_lambda = lambda_peak / R
    xmin = min(x)
    xmax = max(x)
    x_nonbroad = np.arange(xmin, xmax, delta_lambda)
    # using this method for making grid increase fitting time so much
    # x_nonbroad = mg.make_grid(xmin, xmax, resolution=R)
    x_profile = np.array(x_nonbroad)

    # convert vel_cloud to central wavelength
    central_wave = lambda_peak


    # Calculate intrinsic profile
    nu = cst.c.to('cm/s').value / (x_profile * 1.e-8)
    nu0 = cst.c.to('cm/s').value / (central_wave * 1.e-8)
    delta_nu = nu - nu0
    delta_nu_D = (b_eff*1.e5) * nu / cst.c.to('cm/s').value
    prf = 1.0 / ((np.pi**0.5) * delta_nu_D)
    Z_xval = delta_nu / delta_nu_D
    Z_gval = gamma / (4 * np.pi * delta_nu_D)
    vgt = prf * wofz(Z_xval + 1j*Z_gval).real
    tau = (10**log_N) * sigma0 * osc_freq * vgt
    voigt_model = np.exp(-tau) - 1

    # broad the intrinsic line profile by the instrumental resolution
    pxs = np.diff(x_profile)[0] / x_profile[0] * cst.c.to('km/s').value
    fwhm_instrumental = resolution
    sigma_instrumental = fwhm_instrumental / 2.35482 / pxs
    LSF = gaussian( int(len(x_profile)/2), sigma_instrumental)
    LSF = LSF / LSF.sum()
    y_profile = fftconvolve(voigt_model, LSF, 'same')


    # interpolate into the observed wavelength grid
    obs_profile = np.interp(x, x_profile, y_profile)

    return obs_profile
