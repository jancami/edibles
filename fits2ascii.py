from astropy.io import fits
import glob
import os

def fits2ascii(input_fits):

    """ usage:
    spec_wave, spec_flux = fits2ascii("a_edibles_spectrum.fits")
    """

    """ note: probably needs some checks to see if FITS file is valid """

    hdu=fits.open(input_fits)
    spec_flux=hdu[0].data
    crval1=hdu[0].header["CRVAL1"]
    cdelt1=hdu[0].header["CDELT1"]
    nwave     = len(spec_flux)
    wave      = np.arange(0, nwave, 1)
    spec_wave = (wave) * cdelt1 + crval1
    #d=(np.round(spec_wave,4),np.round(spec_flux,6))
    d=(spec_wave,spec_flux)
    #np.savetxt(output_dir+spectrum[:-5]+".ascii",np.array(d).T, fmt=['%12.4f', '%14.6f'])
    return np.array(d).T
