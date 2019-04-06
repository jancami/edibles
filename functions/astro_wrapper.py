import numpy as np
from astropy import constants as cst
from edibles.functions.frequency_convert import converter
from edibles.functions.voigtMathematical import voigt_math
from edibles.functions.freq_convert_2 import unit_converter


def voigt_astro(x, cent, b_eff, Gamma):
	'''
	INPUT:
	x:     [ndarray]	Wavelength grid
	cent:  [float] 		Central wavelength
	b_eff: [float]		Velocity width [km/s]
	Gamma: [float]		Lorentzian HWHM component * 4pi

	OUTPUT:
	x:	   [ndarray]	Wavelength grid in km/s
	y: 	   [ndarray] 	Voigt profile
	'''

	nu, nu0 = converter(x, cent)

	# nu = unit_converter(x, units='Angstroms')
	# nu0 = unit_converter(cent, units='Angstroms')

	gamma = Gamma / (4 * np.pi) * 1e-10


	# b_eff -> delta_nu_D -> sigma -> alpha

	delta_nu_D = nu0 * b_eff / cst.c.to('km/s').value  # freq * km/s / km/s
	sigma = delta_nu_D / np.sqrt(2)
	alpha_Hz = sigma * np.sqrt(2 * np.log(2))


	alpha = alpha_Hz*cst.c.value/(nu0**2*1.e-10)
	# x, cent, alpha, gamma2 = converter(nu, nu0, alpha_Hz, gamma, units='Hz')


	print('b_eff: ', b_eff)
	print('Gamma: ', Gamma)
	print('alpha: ', alpha)
	print('gamma: ', gamma)

	y = voigt_math(x, cent, alpha, gamma)

	return x, y
