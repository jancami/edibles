def barycorrect_spectrum(wave_array,fits_header):

    #import astropy.constants
    #c = constants.c
    c = 299792.458 # speed of light (km/s)

    vbary = fits_header['ESO QC VRAD BARYCOR']
    wave_array = wave_array + (vbary / c)*wave_array

    return wave_array
