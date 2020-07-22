import numpy as np
from astropy.io import fits
import astropy.constants as cst
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from datetime import datetime

from edibles import DATADIR


class EdiblesSpectrum:
    """
    This class takes a spectrum file from EDIBLES,
    reads the header and data, and creates a DataFrame.

    The class will also contain a set of methods to operate on the data.

    Args:
        filename (str): Name of the file, starting with the target
        noDATADIR (bool): If true, DATADIR will not be added to the front of the filename

    Attributes:
        header (astropy.io.fits.header.Header): The header of the FITS file from the observation
        target (str): The name of the target
        datetime (datetime.datetime): The date of the target observation
        v_bary (float): Barycentric velocity of the target star
        wave (1darray): The wavelength data for the spectrum, geocentric reference frame,
            will be updated by the functions
        bary_wave (1darray): The wavelength data for the spectrum, barycentric reference frame,
            will be updated by the functions
        flux (1darray): The flux data for the spectrum,
            will be updated by the functions
        wave_units (str): The units of the wavelength data
        flux_units (str): The units of the flux data
        raw_wave (1darray): The wavelength data for the spectrum, geocentric reference frame,
            will not be updated by the functions
        raw_bary_wave (1darray): The wavelength data for the spectrum, barycentric reference frame,
            will not be updated by the functions
        raw_flux (1darray): The flux data for the spectrum,
            will not be updated by the functions

    """

    def __init__(self, filename, noDATADIR=False):
        """Filename is relative to the EDIBLES_DATADIR environment variable

        """
        self.filename = DATADIR + filename

        if noDATADIR is True:
            self.filename = filename

        self.loadSpectrum()


    def loadSpectrum(self):

        hdu = fits.open(self.filename)
        self.header = hdu[0].header
        self.target = self.header["OBJECT"]
        self.date = self.header["DATE-OBS"]
        self.datetime = datetime.strptime(self.header["DATE-OBS"], '%Y-%m-%dT%H:%M:%S.%f')
        self.v_bary = self.header["HIERARCH ESO QC VRAD BARYCOR"]

        self.flux = hdu[0].data
        crval1 = self.header["CRVAL1"]
        cdelt1 = self.header["CDELT1"]
        lenwave = len(self.flux)
        grid = np.arange(0, lenwave, 1)
        self.wave = (grid) * cdelt1 + crval1
        self.bary_wave = self.wave + (self.v_bary / cst.c.to("km/s").value) * self.wave

        self.raw_wave = (grid) * cdelt1 + crval1
        self.raw_bary_wave = self.wave + (self.v_bary / cst.c.to("km/s").value) * self.wave
        self.raw_flux = hdu[0].data

        self.wave_units = "AA"
        self.flux_units = "arbitrary"


    def spec_grid(self):

        xmin = np.min(self.raw_wave)
        xmax = np.max(self.raw_wave)
        cent = np.mean([xmin, xmax])
        R = 80000
        oversample = 2

        spacing = cent / R / oversample
        grid = np.arange(3000, 10500, spacing)

        return grid


    def getSpectrum(self, xmin=None, xmax=None):
        """Function to update the wavelength region held in an EdiblesSpectrum object.

        Args:
            xmin (float): Minimum wavelength (Optional)
            xmax (float): Maximum wavelength (Optional)

        """

        assert xmin is not None, "xmin is not defined"
        assert xmax is not None, "xmax is not defined"
        assert xmin < xmax, "xmin must be less than xmax"

        t_idx = np.where(np.logical_and(self.wave > xmin, self.wave < xmax))
        self.wave = self.raw_wave[t_idx]
        self.flux = self.raw_flux[t_idx]

        b_idx = np.where(np.logical_and(self.bary_wave > xmin, self.bary_wave < xmax))
        self.bary_wave = self.bary_wave[b_idx]
        self.bary_flux = self.raw_flux[b_idx]

        grid = self.spec_grid()
        interp_idx = np.where(np.logical_and(grid > xmin, grid < xmax))
        self.grid = grid[interp_idx]

        f = interp1d(self.raw_wave, self.raw_flux)
        i_flux = f(self.grid)
        self.interp_flux = i_flux

        bf = interp1d(self.raw_bary_wave, self.raw_flux)
        b_flux = bf(self.grid)
        self.interp_bary_flux = b_flux




if __name__ == "__main__":
    filename = "/HD170740/RED_860/HD170740_w860_redl_20140915_O12.fits"
    sp = EdiblesSpectrum(filename)
    print(sp.target)
    print(sp.datetime.date())
    print("Barycentric Velocity is", sp.v_bary)
    plt.plot(sp.wave, sp.flux, label="Geocentric")
    plt.show()

    sp.getSpectrum(xmin=7660, xmax=7680)

    plt.plot(sp.wave, sp.flux, label="Geocentric Subset")
    plt.plot(sp.bary_wave, sp.bary_flux, label="Barycentric Subset")
    plt.legend()
    plt.show()

    plt.plot(sp.grid, sp.interp_flux, label='Geocentric Interpolation')
    plt.plot(sp.grid, sp.interp_bary_flux, label='Barycentric Interpolation')
    plt.legend()
    plt.show()
