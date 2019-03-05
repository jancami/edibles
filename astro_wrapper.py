import numpy as np
from astropy import constants as cst
from scipy.special import wofz


def voigtMathToAstro(nu, nu0, alpha, gamma):
	'''
	INPUT:
	nu:    [ndarray]	Frequency grid
	nu0:   [float] 		Central frequency
	alpha: [float]		Gaussian HWHM component
	gamma: [float]		Lorentzian HWHM component


	OUTPUT:
	tau: [ndarray] Optical depth array


	'''

	sigma = alpha / np.sqrt(2 * np.log(2))
	delta_nu = nu - nu0
	delta_nu_D = sigma * np.sqrt(2)
	Gamma = 4 * np.pi * gamma
	b = cst.c.value / nu * delta_nu_D

	vgt = 1. / (delta_nu_D * np.sqrt(np.pi)) \
		* wofz((delta_nu / delta_nu_D) + 1j * (Gamma / (4 * np.pi * delta_nu_D)))

	tau = (10**log_N) * sigma0 * osc_freq * vgt
	voigt_model = np.exp(-tau) - 1

	# # Amin's
	# nu = cst.c.to('cm/s').value / (x * 1.e-8)
	# nu0 = cst.c.to('cm/s').value / (cent * 1.e-8)
	# delta_nu = nu - nu0
	# delta_nu_D = (b_eff*1.e5) * nu / cst.c.to('cm/s').value
	# prf = 1.0 / ((np.pi**0.5) * delta_nu_D)
	# Z_xval = delta_nu / delta_nu_D
	# Z_gval = gamma / (4 * np.pi * delta_nu_D)
	# vgt = prf * wofz(Z_xval + 1j*Z_gval).real
	# tau = (10**log_N) * sigma0 * osc_freq * vgt
	# voigt_model = np.exp(-tau) - 1

	# # From Jan's book
	# N = '?'	  	# column density
	# g1 = '?'  	# statistical weight of level 1
	# B12 = '?'  	# Einstein absorption coefficient
	# E1 = '?'  	# Energy of level 1
	# T = '?'  	# Temperature

	# def P(T):
	# 	'''
	# 	Partition function at temp T
	# 	'''
	# 	return '?'

	# tau = N * (cst.h.value * nu0 / (4*np.pi)) 							\
	# 	* g1 * B12 * (np.exp(-E1 / (cst.k.value * T)) / P(T)) 			\
	# 	* (1 - np.exp(-cst.h.value * nu0 / (cst.k.value * T))) * nu


	return voigt_model
