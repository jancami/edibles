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
#   lambda0   : The central wavelength for known atoms/molecules
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
#
# Copyright:
#     The Copyright for these codes belongs to EDIBLES survey.
#     For use it please refer to survey P.I. Nick Cox, Jan Cami or
#     send an email to Amin Farhang.
#     farhang.amin@gmail.com
#
# +

import numpy as np

# =========================================================================
#   Voigt Profile Approximation from T. Tepper-Garcia 2006, 2007
# =========================================================================
def H(a, x):
    P = x**2
    H0 = np.exp(-x**2)
    Q = 1.5/x**2
    H_a_x = H0 - (a/np.sqrt(np.pi)/P) * (H0*H0 * (4.*P*P + 7.*P + 4. + Q) - Q - 1)
    return H_a_x


# =========================================================================
#   Voigt model -- the astronomical version of Voigt function
# =========================================================================
def voigt(x, v_cloud, lambda0, z, b_eff, log_N, f_osc, gam):

    c = 2.99792e10         # cm/s
    m_e = 9.1094e-28       # g
    e = 4.8032e-10         # cgs units

    # convert v_cloud to central wavelength
    cw = (v_cloud / 299792.458) * lambda0 + lambda0

    # Calculate Profile
    C_a = (np.sqrt(np.pi) * e**2 * f_osc * (cw*1.e-8)) / (m_e * c * (b_eff* 1.e5))
    aa = ((cw*1.e-8) * gam) / (4. * np.pi * (b_eff* 1.e5))

    dl_D = (b_eff * 1.e5)/c * cw
    x = x / (z + 1.0)
    xx = (x - cw)/dl_D + 0.00001

    tau = np.float64(C_a) * 10**log_N * H(aa, xx)
    voigt_model = np.exp(-tau) - 1.0
    return voigt_model
