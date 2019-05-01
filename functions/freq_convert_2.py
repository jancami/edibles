import astropy.units as u

def unit_converter(param, units):
    '''
    INPUT
    param: [ndarray OR float OR integer] value or list of values to convert
    units: [str] must be 'Angstroms' or 'Hz'

    OUTPUT
    converted_param: [ndarray OR float OR int] converted value or list of values
    '''

    if units == 'Hz':
        converted_param = (param * u.Hz).to(u.AA, equivalencies=u.spectral()).value

    if units == 'Angstroms':
        converted_param = (param * u.AA).to(u.Hz, equivalencies=u.spectral()).value

    return converted_param