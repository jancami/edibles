import numpy as np
from astropy.io import fits
import astropy.constants as cst
import matplotlib.pyplot as plt
import pandas as pd

from edibles.edibles import DATADIR


class EdiblesSpectrum:
    """
    This class takes a spectrum file from EDIBLES,
    reads the header and data, and creates a DataFrame.

    The class will also contain a set of methods to operate on the data.

    :param filename: Name of the file, starting with the target
    :type filename: str
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

    def __init__(self, filename):
        """
        Filename is relative to the DR3 directory
        """
        self.filename = DATADIR + filename
        self.loadSpectrum()

    def loadSpectrum(self):
        # Assume the file is a DR3 product here.
        hdu = fits.open(self.filename)
        self.header = hdu[0].header
        self.target = self.header["OBJECT"]
        self.date = self.header["DATE-OBS"]

        flux = hdu[0].data
        crval1 = self.header["CRVAL1"]
        cdelt1 = self.header["CDELT1"]
        lenwave = len(flux)
        grid = np.arange(0, lenwave, 1)
        wave = (grid) * cdelt1 + crval1
        self.v_bary = self.header["HIERARCH ESO QC VRAD BARYCOR"]
        bary_wave = wave + (self.v_bary / cst.c.to("km/s").value) * wave

        d = {
            "wave": wave.tolist(),
            "bary_wave": bary_wave.tolist(),
            "flux": flux.tolist(),
        }
        self.df = pd.DataFrame(data=d)

        self.wave = self.df["wave"]
        self.wave_units = "AA"

        self.bary_wave = self.df["bary_wave"]

        self.flux = self.df["flux"]
        self.flux_units = "arbitrary"

    def getSpectrum(self, xmin=None, xmax=None):
        """
        Function to get the wavelength and flux arrays of a particular target.
        If xmin/xmax are not called, the data for the entire spectrum will be returned.

        Args:
            xmin (float): minimum wavelength (Optional)
            xmax (float): Maximum wavelength (Optional)
            bary (bool): Barycentric rest frame, default=False

        Returns:
            ndarray: wavelength grid
            ndarray: flux grid

        """

        if (xmin is not None) and (xmax is not None):
            assert xmin < xmax, "xmin must be less than xmax"

            df_subset = self.df[self.df["wave"].between(xmin, xmax)]

            return df_subset

        return self.df


if __name__ == "__main__":
    filename = "/HD170740/RED_860/HD170740_w860_redl_20140915_O12.fits"
    sp = EdiblesSpectrum(filename)
    print(sp.target)
    print("Barycentric Velocity is", sp.v_bary)
    plt.figure()
    plt.plot(sp.wave, sp.flux, label="Geocentric")
    plt.figure()

    subset = sp.getSpectrum(xmin=7660, xmax=7680)
    plt.plot(subset["wave"], subset["flux"], label="Geocentric Subset")
    plt.plot(subset["bary_wave"], subset["flux"], label="Barycentric Subset")
    plt.legend()
    plt.show()
