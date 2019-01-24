from astropy.io import fits
import numpy as np

def fits2ascii(input_fits):

    """ usage:
    fits2ascii("an_edibles_spectrum.fits")
    """

    """ note: probably needs some error catching checks to see if FITS file is valid """

    hdu = fits.open(input_fits)
    spec_flux = hdu[0].data
    crval1 = hdu[0].header["CRVAL1"]
    cdelt1 = hdu[0].header["CDELT1"]
    nwave = len(spec_flux)
    wave = np.arange(0, nwave, 1)
    spec_wave = (wave) * cdelt1 + crval1
    #d=(np.round(spec_wave,4),np.round(spec_flux,6))
    d=(spec_wave,spec_flux)
    np.savetxt(input_fits[:-5]+".ascii",np.array(d).T, fmt=['%12.4f', '%14.6f'])
    #return np.array(d).T
    print("created ascii spectrum: "+input_fits[:-5]+".ascii")
