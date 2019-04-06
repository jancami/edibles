import astropy.units as u

def unit_converter(array, units):
	'''
	INPUT
	array: [ndarray OR float OR integer] value or listr of values to convert
	units: [str] must be 'Angstroms' or 'Hz'

	OUTPUT
	converted_array: [ndarray OR float OR int] converted value or list of values

	'''

	if units == 'Hz':
		converted_array = (array * u.Hz).to(u.AA, equivalencies=u.spectral()).value


	if units == 'Angstroms':
		converted_array = (array * u.AA).to(u.Hz, equivalencies=u.spectral()).value

	return converted_array