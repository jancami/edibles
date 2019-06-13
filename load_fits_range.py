import numpy as np
from astropy.io import fits


def load_fits_range(path, xmin, xmax):

    assert xmin < xmax, 'xmin must be less than xmax'


    hdu = fits.open(path)
    spec_flux = hdu[0].data
    crval1 = hdu[0].header["CRVAL1"]
    cdelt1 = hdu[0].header["CDELT1"]
    nwave = len(spec_flux)
    grid = np.arange(0, nwave, 1)
    spec_wave = (grid) * cdelt1 + crval1

    # create data subset
    min_idx = (np.abs(spec_wave - xmin)).argmin()
    max_idx = (np.abs(spec_wave - xmax)).argmin()
    wave = spec_wave[min_idx:max_idx]
    flux = spec_flux[min_idx:max_idx]

    data = (wave, flux)

    return data



