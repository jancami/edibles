import numpy as np
from astropy.io import fits
import astropy.constants as cst
import matplotlib.pyplot as plt
import pandas as pd

from edibles import DATADIR


class EdiblesSpectrum:
    """
    This class takes a spectrum file from EDIBLES,
    reads the header and data, and creates a DataFrame.

    The class will also contain a set of methods to operate on the data.

    :param filename: Name of the file, starting with the target
    :type filename: str
    :param noDATADIR: If true, DATADIR will not be added to the front of the filename
    :type noDATADIR: Bool
    :param header: The header of the FITS file from the target observation
    :type header: Object (astropy.io.fits.header.Header)
    :param target: The name of the target
    :type target: str
    :param date: Date of the target observation
    :type date: str
    :param v_bary: Barycentric velocity of the target star
    :type v_bary: float
    :param df: Pandas array containing geocentric and barycentric wavelength, and flux
    :type df: Pandas array (pandas.core.series.Series)
    :param wave: The wavelength grid for the spectrum, geocentric reference frame
    :type wave: Pandas array (pandas.core.series.Series)
    :param wave_units: The units of the wavelength array
    :type wave_units: str
    :param bary_wave: The wavelength grid for the spectrum, barycentric reference frame
    :type bary_wave: Pandas array (pandas.core.series.Series)
    :param flux: The flux data for the spectrum
    :type flux: Pandas array (pandas.core.series.Series)
    :param flux_units: The units of the flux data
    :type flux_units: str

    """

    def __init__(self, filename, noDATADIR=False):
        """
        Filename is relative to the DR3 directory
        """
        self.filename = DATADIR + filename

        if noDATADIR is True:
            self.filename = filename

        self.loadSpectrum()

    def loadSpectrum(self):
        # Assume the file is a DR3 product here.
        hdu = fits.open(self.filename)
        self.header = hdu[0].header
        self.target = self.header["OBJECT"]
        self.date = self.header["DATE-OBS"]

        self.flux = hdu[0].data
        crval1 = self.header["CRVAL1"]
        cdelt1 = self.header["CDELT1"]
        lenwave = len(self.flux)
        grid = np.arange(0, lenwave, 1)
        self.wave = (grid) * cdelt1 + crval1
        self.v_bary = self.header["HIERARCH ESO QC VRAD BARYCOR"]
        self.bary_wave = self.wave + (self.v_bary / cst.c.to("km/s").value) * self.wave

        d = {
            "wave": self.wave,
            "bary_wave": self.bary_wave,
            "flux": self.flux,
        }
        self.df = pd.DataFrame(data=d)

        self.wave_units = "AA"
        self.flux_units = "arbitrary"


    def getSpectrum(self, xmin=None, xmax=None):
        """Function to update the wavelength region held in an EdiblesSpectrum object.

        Args:
            xmin (float): minimum wavelength (Optional)
            xmax (float): Maximum wavelength (Optional)

        """

        assert xmin is not None, "xmin is not defined"
        assert xmax is not None, "xmax is not defined"
        assert xmin < xmax, "xmin must be less than xmax"

        idx = np.where(np.logical_and(self.wave > xmin, self.wave < xmax))
        idxmin = np.min(idx)
        idxmax = np.max(idx)

        self.wave = self.wave[idxmin:idxmax]
        self.bary_wave = self.bary_wave[idxmin:idxmax]
        self.flux = self.flux[idxmin:idxmax]
        self.df = self.df[idxmin:idxmax]


        return self.df


if __name__ == "__main__":
    filename = "/HD170740/RED_860/HD170740_w860_redl_20140915_O12.fits"
    sp = EdiblesSpectrum(filename)
    print(sp.target)
    print("Barycentric Velocity is", sp.v_bary)
    plt.plot(sp.wave, sp.flux, label="Geocentric")
    plt.show()

    sp.getSpectrum(xmin=7660, xmax=7680)
    plt.plot(sp.wave, sp.flux, label="Geocentric Subset")
    plt.plot(sp.bary_wave, sp.flux, label="Barycentric Subset")
    plt.legend()
    plt.show()
