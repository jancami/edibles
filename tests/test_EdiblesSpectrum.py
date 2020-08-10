import astropy
import datetime
import numpy as np
from edibles.utils.edibles_spectrum import EdiblesSpectrum


def testEdiblesSpectrum(filename="tests/HD170740_w860_redl_20140915_O12.fits"):

    # Spectrum information
    sp = EdiblesSpectrum(filename=filename, noDATADIR=True)
    assert isinstance(sp.header, astropy.io.fits.header.Header)
    assert isinstance(sp.target, str)
    assert isinstance(sp.date, str)
    assert isinstance(sp.datetime, datetime.datetime)
    assert isinstance(sp.v_bary, float)
    assert isinstance(sp.wave_units, str)
    assert isinstance(sp.flux_units, str)

    # Raw
    assert isinstance(sp.raw_wave, np.ndarray)
    assert isinstance(sp.raw_bary_wave, np.ndarray)
    assert isinstance(sp.raw_flux, np.ndarray)
    assert len(sp.raw_wave) == len(sp.raw_bary_wave)
    assert len(sp.raw_wave) == len(sp.raw_flux)

    assert isinstance(sp.raw_grid, np.ndarray)
    assert len(sp.raw_grid) == 200443  # print(len(sp.raw_grid))

    assert isinstance(sp.raw_sky_wave, np.ndarray)
    assert isinstance(sp.raw_sky_flux, np.ndarray)
    assert len(sp.raw_sky_wave) == len(sp.raw_sky_flux)

    assert isinstance(sp.wave, np.ndarray)
    assert isinstance(sp.bary_wave, np.ndarray)
    assert isinstance(sp.flux, np.ndarray)


    # getSpectrum
    xmin = 7660
    xmax = 7680
    sp.getSpectrum(xmin=xmin, xmax=xmax)
    assert xmin == sp.xmin
    assert xmax == sp.xmax
    assert isinstance(sp.wave, np.ndarray)
    assert isinstance(sp.flux, np.ndarray)
    assert len(sp.wave) == len(sp.flux)
    assert np.min(sp.wave) > sp.xmin
    assert np.max(sp.wave) < sp.xmax

    assert isinstance(sp.bary_wave, np.ndarray)
    assert isinstance(sp.bary_flux, np.ndarray)
    assert len(sp.bary_wave) == len(sp.bary_flux)
    assert np.min(sp.bary_wave) > sp.xmin
    assert np.max(sp.bary_wave) < sp.xmax

    assert isinstance(sp.grid, np.ndarray)
    assert isinstance(sp.interp_flux, np.ndarray)
    assert isinstance(sp.interp_bary_flux, np.ndarray)
    assert len(sp.grid) == len(sp.interp_flux)
    assert len(sp.grid) == len(sp.interp_bary_flux)
    assert np.min(sp.grid) > sp.xmin
    assert np.max(sp.grid) < sp.xmax

    assert isinstance(sp.sky_wave, np.ndarray)
    assert isinstance(sp.sky_flux, np.ndarray)
    assert len(sp.sky_wave) == len(sp.sky_flux)
    assert np.min(sp.sky_wave) > sp.xmin
    assert np.max(sp.sky_wave) < sp.xmax

    # shift
    zoom_xmin = 7661
    zoom_xmax = 7679
    shift = 0.05
    sp.shift(shift=shift, zoom_xmin=zoom_xmin, zoom_xmax=zoom_xmax)

    assert isinstance(sp.wave, np.ndarray)
    assert isinstance(sp.flux, np.ndarray)
    assert len(sp.wave) == len(sp.flux)
    assert np.min(sp.wave) > sp.xmin
    assert np.max(sp.wave) < sp.xmax

    assert isinstance(sp.bary_wave, np.ndarray)
    assert isinstance(sp.bary_flux, np.ndarray)
    assert len(sp.bary_wave) == len(sp.bary_flux)
    assert np.min(sp.bary_wave) > sp.xmin
    assert np.max(sp.bary_wave) < sp.xmax

    assert isinstance(sp.grid, np.ndarray)
    assert isinstance(sp.interp_flux, np.ndarray)
    assert isinstance(sp.interp_bary_flux, np.ndarray)
    assert len(sp.grid) == len(sp.interp_flux)
    assert len(sp.grid) == len(sp.interp_bary_flux)
    assert np.min(sp.grid) > sp.xmin
    assert np.max(sp.grid) < sp.xmax

    assert isinstance(sp.sky_wave, np.ndarray)
    assert isinstance(sp.sky_flux, np.ndarray)
    assert len(sp.sky_wave) == len(sp.sky_flux)
    assert np.min(sp.sky_wave) > sp.xmin
    assert np.max(sp.sky_wave) < sp.xmax







if __name__ == "__main__":

    filename = "HD170740_w860_redl_20140915_O12.fits"
    testEdiblesSpectrum(filename=filename)
