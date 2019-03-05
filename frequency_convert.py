from astropy import constants as cst

def convert_to_nu(x, cent, alpha, gamma):
	'''
	Short function that converts wavelength to frequency

	INPUT:  [Angstroms]

	x:      [ndarray]  Wavelength data grid
    cent:   [float]    Wavelength peak of the Voigt profile
	

	OUTPUT: [Hz]

	nu:		[ndarray]  Frequency data grid
	nu0:	[float]    Frequency peak of the Voigt profile
	'''

	nu = cst.c.value / (x * 1.e-10)
	nu0 = cst.c.value / (cent * 1.e-10)
	a1 = alpha*nu0**2/cst.c.value*1.e-10
	g1 = gamma*nu0**2/cst.c.value*1.e-10

	return nu, nu0, a1, g1
