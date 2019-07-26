
def vac2air_morton(vacw):
    """ Convert vacuum wavelengths in Angstroms to air wavelengths.
    
    This uses the relation from Morton 1991, ApJS, 77, 119. Only valid
    for wavelengths > 2000 Ang.  Use this for compatibility with older
    spectra that may have been corrected using the (older) Morton
    relation.  The Ciddor relation used in vac2air_ciddor() is claimed
    to be more accurate at IR wavelengths.
    """
    temp = (1e4 / vacw) ** 2
    airw = 1. / (1. + 6.4328e-5 + 2.94981e-2/(146 - temp) +
                 2.5540e-4/(41 - temp)) * vacw
    return airw


""" example:

wtc = vac2air_morton(wt*10.0) #make sure wave is in \AA.

"""
