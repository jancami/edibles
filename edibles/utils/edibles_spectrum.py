import numpy as np
from astropy.io import fits
import astropy.constants as cst
import astropy.units as u
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from datetime import datetime
from specutils.utils.wcs_utils import vac_to_air

from edibles import DATADIR
from edibles import PYTHONDIR
from edibles.utils.functions import make_grid


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
        raw_wave (1darray): The wavelength data for the spectrum, geocentric reference frame,
            will not be updated by the functions
        raw_bary_wave (1darray): The wavelength data for the spectrum, barycentric reference frame,
            will not be updated by the functions
        raw_flux (1darray): The flux data for the spectrum,
            will not be updated by the functions
        raw_grid (1darray): A grid covering the entire spectral range used for interpolation
        raw_sky_wave (1darray): Telluric transmission data covering the entire spectral range
        raw_sky_flux (1darray): Telluric transmission data covering the entire spectral range
        wave (1darray): The wavelength data for the spectrum, geocentric reference frame,
            will be updated by the functions
        bary_wave (1darray): The wavelength data for the spectrum, barycentric reference frame,
            will be updated by the functions
        flux (1darray): The flux data for the spectrum,
            will be updated by the functions
        xmin (float): minimum wavelength boundary of data subset - input to getSpectrum
        xmax (float): maximum wavelength boundary of data subset - input to getSpectrum
        sky_wave (1darray): Telluric transmission data - created and updated by getSpectrum
        sky_flux (1darray): Telluric transmission data - created and updated by getSpectrum
        grid (1darray): Interpolation grid - created by _interpolate
        interp_flux (1darray): Interpolated geocentric flux - created by _interpolate
        interp_bary_flux (1darray): Interpolated barycentric flux - created by _interpolate
        wave_units (str): The units of the wavelength data
        flux_units (str): The units of the flux data

    """

    def __init__(self, filename, noDATADIR=False):
        """Filename is relative to the EDIBLES_DATADIR environment variable

        """
        self.filename = DATADIR + filename

        if noDATADIR is True:
            self.filename = filename

        self._loadSpectrum()
        self._spec_grid()
        self._sky_transmission()


    def _loadSpectrum(self):

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
        self.raw_bary_wave = self.raw_wave + (self.v_bary / cst.c.to("km/s").value) * self.raw_wave
        self.raw_flux = hdu[0].data

        self.wave_units = "AA"
        self.flux_units = "arbitrary"


    def _spec_grid(self):
        '''Creates a grid used for interpolation.

        '''
        grid = make_grid(3000, 10500, resolution=80000, oversample=2)
        self.raw_grid = grid


    def _sky_transmission(self):
        '''A function that adds the telluric transmission data to the EdiblesSpectrum model.

        '''
        filename = PYTHONDIR + "/edibles/data/auxillary_data/sky_transmission/transmission.dat"
        sky_transmission = np.loadtxt(filename)

        vac_wave = sky_transmission[:, 0] * 10
        sky_wave = vac_to_air(vac_wave * u.AA, method='Ciddor1996').value

        sky_flux = sky_transmission[:, 1]

        self.raw_sky_wave = sky_wave
        self.raw_sky_flux = sky_flux


    def getSpectrum(self, xmin, xmax):
        """Function to update the wavelength region held in an EdiblesSpectrum object.

        Args:
            xmin (float): Minimum wavelength
            xmax (float): Maximum wavelength

        """
        assert xmin < xmax, "xmin must be less than xmax"
        assert xmin > np.min(self.raw_wave), "xmin outside bounds"
        assert xmax < np.max(self.raw_wave), "xmax outside bounds"

        self.xmin = xmin
        self.xmax = xmax

        # Geocentric data
        t_idx = np.where(np.logical_and(self.raw_wave > xmin, self.raw_wave < xmax))
        self.wave = self.raw_wave[t_idx]
        self.flux = self.raw_flux[t_idx]

        # Barycentric data
        b_idx = np.where(np.logical_and(self.raw_bary_wave > xmin, self.raw_bary_wave < xmax))
        self.bary_wave = self.raw_bary_wave[b_idx]
        self.bary_flux = self.raw_flux[b_idx]

        # Sky transmission data
        sky_idx = np.where(np.logical_and(self.raw_sky_wave > xmin, self.raw_sky_wave < xmax))
        self.sky_wave = self.raw_sky_wave[sky_idx]
        self.sky_flux = self.raw_sky_flux[sky_idx]

        self._interpolate(initial=True)


    def _interpolate(self, initial=False):
        '''Interpolation function used in shift().

        Warning:
            This function has not been tested for use on its own.
            It is implemented in getSpectrum and shift.

        '''
        xmin = np.max([np.min(self.wave), np.min(self.bary_wave)])
        xmax = np.min([np.max(self.wave), np.max(self.bary_wave)])

        grid_idx = np.where(np.logical_and(self.raw_grid > xmin,
                                           self.raw_grid < xmax))
        self.grid = self.raw_grid[grid_idx]

        if initial:
            # Interpolate geocentric flux data
            f = interp1d(self.raw_wave, self.raw_flux)
            self.interp_flux = f(self.grid)

            # Interpolate barycentric flux data
            bf = interp1d(self.raw_bary_wave, self.raw_flux)
            self.interp_bary_flux = bf(self.grid)

        else:
            # Interpolate geocentric flux data
            func = interp1d(self.wave, self.flux)
            self.interp_flux = func(self.grid)

            # Interpolate barycentric flux data
            bfunc = interp1d(self.bary_wave, self.bary_flux)
            self.interp_bary_flux = bfunc(self.grid)


    def shift(self, shift, zoom_xmin, zoom_xmax):
        '''Shift the geocentric and update the barycentric wavelength data.
Reinterpolates after shifting.

        Args:
            shift (float or 1darray): Amount to shift wavelength grid by.
If shift is an array, it must be the same length as the wavelength grid.
            zoom_xmin (float): New minimum wavelength (must be > old xmin)
            zoom_xmax (float): New maximum wavelength (must be < old xmax)

        '''

        # Add shift to wavelength data
        new_wave = self.wave + shift

        # Input checking
        assert (zoom_xmin > self.xmin) and (
            zoom_xmin > np.min(self.wave)
        ), 'zoom_xmin must be greater than ' + str(np.max([self.xmin, np.min(new_wave)]))
        assert (zoom_xmax < self.xmax) and (
            zoom_xmax < np.max(self.wave)
        ), 'zoom_xmax must be less than ' + str(np.min([self.xmax, np.max(new_wave)]))
        if isinstance(shift, np.ndarray):
            assert len(self.wave) == len(shift), 'Length of shift not equal to length of wave'

        self.wave = new_wave
        self.bary_wave = new_wave + (self.v_bary / cst.c.to("km/s").value) * new_wave

        # Cutoff 'edges' of data
        b_idx = np.where(np.logical_and(self.bary_wave > zoom_xmin, self.bary_wave < zoom_xmax))
        self.bary_wave = self.bary_wave[b_idx]
        self.bary_flux = self.flux[b_idx]

        t_idx = np.where(np.logical_and(self.wave > zoom_xmin, self.wave < zoom_xmax))
        self.wave = self.wave[t_idx]
        self.flux = self.flux[t_idx]

        # Sky transmission data
        sky_idx = np.where(np.logical_and(self.raw_sky_wave > zoom_xmin,
                                          self.raw_sky_wave < zoom_xmax))
        self.sky_wave = self.raw_sky_wave[sky_idx]
        self.sky_flux = self.raw_sky_flux[sky_idx]

        self._interpolate()

        self.xmin = zoom_xmin
        self.xmax = zoom_xmax


if __name__ == "__main__":
    filename = "/HD170740/RED_860/HD170740_w860_redl_20140915_O12.fits"
    sp = EdiblesSpectrum(filename)
    print(sp.target)
    print(sp.datetime.date())
    print("Barycentric Velocity is", sp.v_bary)
    plt.plot(sp.wave, sp.flux, label="Geocentric")
    plt.title('Entire Order')
    plt.xlabel(r'Wavelength ($\AA$)')
    plt.ylabel('Flux')
    plt.legend()
    plt.show()

    sp.getSpectrum(xmin=7660, xmax=7680)
    plt.plot(sp.wave, sp.flux, label="Geocentric Subset")
    plt.plot(sp.bary_wave, sp.bary_flux, label="Barycentric Subset")
    plt.plot(sp.grid, sp.interp_flux, label='Geocentric Interpolation')
    plt.plot(sp.grid, sp.interp_bary_flux, label='Barycentric Interpolation')
    plt.plot(sp.sky_wave, sp.sky_flux, 'k')
    plt.title('Data and Interpolations')
    plt.xlabel(r'Wavelength ($\AA$)')
    plt.ylabel('Flux')
    plt.legend()
    plt.show()

    plt.plot(sp.wave, sp.flux, label="Geocentric")
    plt.plot(sp.bary_wave, sp.bary_flux, label="Barycentric")
    shift = 0.05
    sp.shift(shift=shift, zoom_xmin=7661, zoom_xmax=7679)
    plt.plot(sp.wave, sp.flux, label='Shifted Geocentric')
    plt.plot(sp.bary_wave, sp.bary_flux, label="Shifted Barycentric")
    plt.title('Data and Shifted Data')
    plt.xlabel(r'Wavelength ($\AA$)')
    plt.ylabel('Flux')
    plt.legend()
    plt.show()

    plt.plot(sp.wave, sp.flux, label='Shifted Geocentric')
    plt.plot(sp.bary_wave, sp.bary_flux, label="Shifted Barycentric")
    plt.plot(sp.grid, sp.interp_flux, label='Shifted Geocentric Interpolation')
    plt.plot(sp.grid, sp.interp_bary_flux, label='Shifted Barycentric Interpolation')
    plt.title('Shifted Data and Shifted Interpolations')
    plt.xlabel(r'Wavelength ($\AA$)')
    plt.ylabel('Flux')
    plt.legend()
    plt.show()
