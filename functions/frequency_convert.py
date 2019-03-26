from astropy import constants as cst

def converter(x, cent, alpha=None, gamma=None, units=None):
	'''
	Short function that converts wavelength to frequency

	INPUT:  [Angstroms]

	x:      [ndarray]  Wavelength data grid
	cent:   [float]    Wavelength peak of the Voigt profile
	type:   [str] 	   Angstroms or Hz
	

	OUTPUT: [Hz]

	nu:		[ndarray]  Frequency data grid
	nu0:	[float]    Frequency peak of the Voigt profile
	'''

	a1 = None
	g1 = None


	nu = cst.c.value / (x * 1.e-10)
	nu0 = cst.c.value / (cent * 1.e-10)

	if alpha is not None:
		if units == 'Angstroms':
			a1 = alpha*nu0**2/cst.c.value*1.e-10
		if units == 'Hz':
			a1 = alpha*cst.c.value/(nu0**2*1.e-10)


	if gamma is not None:

		if units == 'Angstroms':
			g1 = gamma*nu0**2/cst.c.value*1.e-10

		if units == 'Hz':
			g1 = gamma*cst.c.value/(nu0**2*1.e-10)

	
	if a1 is None:
		if g1 is None:
			return nu, nu0

		# a1 is None, g1 is not None
		return nu, nu0, g1

	# beyond here, a1 is not None
	if g1 is None:
		return nu, nu0, a1


	return nu, nu0, a1, g1





	# try:
	# 	return nu, nu0, a1, g1
	# except UnboundLocalError:
	# 	return nu, nu0, a1
	# except UnboundLocalError:
	# 	return nu, nu0, g1
	# except UnboundLocalError:
	# 	return nu, nu0
