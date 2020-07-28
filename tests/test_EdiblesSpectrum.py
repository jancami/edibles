import astropy
import pytest
import numpy as np
from edibles.utils.edibles_spectrum import EdiblesSpectrum


def testEdiblesSpectrum(filename="tests/HD170740_w860_redl_20140915_O12.fits"):

    sp = EdiblesSpectrum(filename=filename, noDATADIR=True)
    assert isinstance(sp.header, astropy.io.fits.header.Header)
    assert isinstance(sp.target, str)
    assert isinstance(sp.date, str)
    assert isinstance(sp.v_bary, float)

    print(type(sp.wave))
    assert isinstance(sp.wave, np.ndarray)
    assert isinstance(sp.wave_units, str)
    assert isinstance(sp.bary_wave, np.ndarray)
    assert isinstance(sp.flux, np.ndarray)
    assert isinstance(sp.flux_units, str)


    with pytest.raises(AssertionError):
        assert sp.getSpectrum()

    xmin = 7660
    xmax = 7680
    sp.getSpectrum(xmin=xmin, xmax=xmax)
    assert isinstance(sp.wave, np.ndarray)
    assert isinstance(sp.flux, np.ndarray)

    assert np.min(sp.wave) > xmin
    assert np.max(sp.wave) < xmax



if __name__ == "__main__":

    filename = "HD170740_w860_redl_20140915_O12.fits"
    testEdiblesSpectrum(filename=filename)
