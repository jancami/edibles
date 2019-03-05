from astropy import constants as cst

def convert_to_nu(x, cent):
	'''
	Short function that converts wavelength to frequency

	INPUT:

	x:      [ndarray]  Wavelength data grid
    cent:   [float]    Wavelength peak of the Voigt profile
	
	OUTPUT:

	nu:		[ndarray]  Frequency data grid
	nu0:	[float]    Frequency peak of the Voigt profile
	'''

	nu = cst.c.value / (x * 1.e-10)
	nu0 = cst.c.value / (cent * 1.e-10)

	return nu, nu0
