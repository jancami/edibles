import numpy as np
from astropy import constants as cst

import sys
import os
path = os.getcwd()
os.chdir('..')
sys.path.append(os.getcwd())
os.chdir(path)

from edibles.frequency_convert import converter
from edibles.voigtMathematical import voigt_math




def voigt_astro(x, cent, b_eff, Gamma):
	'''
	INPUT:
	x:     [ndarray]	Wavelength grid
	cent:  [float] 		Central wavelength
	alpha: [float]		Gaussian HWHM component
	gamma: [float]		Lorentzian HWHM component


	OUTPUT:
	tau: [ndarray] Optical depth array


	'''

	nu, nu0 = converter(x, cent)

	# b_eff -> delta_nu_D -> sigma -> alpha

	delta_nu_D = nu0 * b_eff / cst.c.value
	sigma = delta_nu_D / np.sqrt(2)
	alpha = sigma * np.sqrt(2 * np.log(2))


	gamma = Gamma / (4 * np.pi)


	print(alpha)
	print(gamma)

	y = voigt_math(x, cent, alpha, gamma)

	# delta_nu = nu - nu0

	return y
