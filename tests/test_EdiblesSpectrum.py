import pandas
import astropy
import numpy as np
from edibles.edibles.utils.edibles_spectrum import EdiblesSpectrum


def testEdiblesSpectrum():

    sp = EdiblesSpectrum("tests/testdata/HD170740_w860_redl_20140915_O12.fits", noDATADIR=True)
    assert isinstance(sp.header, astropy.io.fits.header.Header)
    assert isinstance(sp.target, str)
    assert isinstance(sp.date, str)
    assert isinstance(sp.v_bary, float)

    assert isinstance(sp.df, pandas.core.frame.DataFrame)
    assert isinstance(sp.wave, pandas.core.frame.Series)
    assert isinstance(sp.wave_units, str)
    assert isinstance(sp.bary_wave, pandas.core.frame.Series)
    assert isinstance(sp.flux, pandas.core.frame.Series)
    assert isinstance(sp.flux_units, str)

    subset1 = sp.getSpectrum()
    assert isinstance(subset1, pandas.core.frame.DataFrame)
    assert isinstance(subset1.wave, pandas.core.frame.Series)
    assert isinstance(subset1.flux, pandas.core.frame.Series)
    assert sp.df.equals(subset1)

    xmin = 7660
    xmax = 7680
    subset2 = sp.getSpectrum(xmin=xmin, xmax=xmax)
    assert isinstance(subset2, pandas.core.frame.DataFrame)
    assert isinstance(subset2.wave, pandas.core.frame.Series)
    assert isinstance(subset2.flux, pandas.core.frame.Series)

    assert np.min(subset2.wave) > xmin
    assert np.max(subset2.wave) < xmax


if __name__ == "__main__":
    testEdiblesSpectrum()
