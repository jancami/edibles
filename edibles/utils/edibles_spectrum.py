import os
import glob
import numpy as np
from astropy.io import fits
import astropy.constants as cst
import astropy.units as u
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d
from datetime import datetime
from specutils.utils.wcs_utils import vac_to_air

from pathlib import Path
from edibles import EDIBLES_PYTHONDIR
from edibles import DATADIR, DATARELEASE


from edibles.utils.functions import make_grid


class EdiblesSpectrum:
    """
    This class takes a spectrum file from EDIBLES,
    reads the header and data, and creates a DataFrame.
    The class will also contain a set of methods to operate on the data.

    Args:
        filename (str): Name of the file, starting with the target
        fully_featured (bool): If true, EdiblesSpectrum generates the gky transmission and
            corrected spectrum
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
        corrected_wave (1darray): The wavelength data for the telluric corrected spectrum,
            created by corrected_spectrum
        flux_initial (1darray): The initial flux data for the telluric corrected spectrum,
            created by corrected_spectrum
        flux_corrO2 (1array): The O2 corrected flux data for the telluric corrected spectrum,
            created by corrected_spectrum
        flux_corrO2_H2O (1array): O2 and H2O corrected flux for the telluric corrected spectrum,
            created by corrected_spectrum
        wave_units (str): The units of the wavelength data
        flux_units (str): The units of the flux data
        continuum_filename (astropy.io.fits.header.Header): Name of file with continuum points

    """

    def __init__(self, filename, fully_featured=False, noDATADIR=False):
        """Filename is relative to the EDIBLES_DATADIR environment variable

        """
        if filename.startswith('/'):
            filename = filename[1:]
        self.filename = Path(DATADIR) / filename
        self.fully_featured = fully_featured

        if noDATADIR is True:
            self.filename = filename

        self._loadSpectrum()
        self._spec_grid()
        if self.fully_featured:
            self._sky_transmission()
            self._corrected_spectrum()

    def _loadSpectrum(self):
        with fits.open(self.filename) as hdulist:
            self.header = hdulist[0].header
            self.target = self.header["OBJECT"]
            self.date = self.header["DATE-OBS"]
            self.datetime = datetime.strptime(self.header["DATE-OBS"], '%Y-%m-%dT%H:%M:%S.%f')
            self.v_bary = self.header["HIERARCH ESO QC VRAD BARYCOR"]

            if DATARELEASE == 'DR5':
                self.wave = hdulist[1].data['WAVE']
                self.flux = hdulist[1].data['FLUX']
                self.raw_wave = np.copy(self.wave)
                self.raw_flux = np.copy(self.flux)

            else:
                self.flux = hdulist[0].data
                crval1 = self.header["CRVAL1"]
                cdelt1 = self.header["CDELT1"]
                lenwave = len(self.flux)
                grid = np.arange(0, lenwave, 1)
                self.wave = (grid) * cdelt1 + crval1
                self.raw_wave = (grid) * cdelt1 + crval1
                self.raw_flux = hdulist[0].data
                
            self.raw_bary_wave = self.raw_wave + \
                (self.v_bary / cst.c.to("km/s").value) * \
                self.raw_wave

            self.bary_wave = self.wave + (self.v_bary / cst.c.to("km/s").value) * self.wave
            self.wave_units = "AA"
            self.flux_units = "arbitrary"

            csv_file = str(self.filename).replace(".fits", ".csv").replace(
                "/DR4/data/", "/DR4/continuum/").replace(r"\DR4\data", r"\DR4\continuum")

            if os.path.isfile(csv_file):
                self.continuum_filename = csv_file

    def _spec_grid(self):
        '''Creates a grid used for interpolation.

        '''
        grid = make_grid(3000, 10500, resolution=80000, oversample=2)
        self.raw_grid = grid

    def _sky_transmission(self):
        '''A function that adds the telluric transmission data to the EdiblesSpectrum model.

        '''
        filename = EDIBLES_PYTHONDIR / 'data/auxiliary_data/sky_transmission/transmission.dat'
        #print(filename)
        sky_transmission = np.loadtxt(filename)

        vac_wave = sky_transmission[:, 0] * 10
        sky_wave = vac_to_air(vac_wave * u.AA, method='Ciddor1996').value

        sky_flux = sky_transmission[:, 1]

        self.raw_sky_wave = sky_wave
        self.raw_sky_flux = sky_flux

    def _corrected_spectrum(self):
        '''A function that adds the telluric corrected spectrum data to the EdiblesSpectrum model.

        '''
        #print(self.datetime.date())

        stripped_date = str(self.datetime.date()).replace('-', '')

        search_path = EDIBLES_PYTHONDIR / 'data/telluric_corrected_data'      
        search_string = self.target + "*" + stripped_date + "*.ascii"
        #print(search_path)
        #print(search_string)
        filename = list(search_path.glob(search_string))
        #print(filename)
        #filename = glob.glob(
        #    PYTHONDIR + "/data/telluric_corrected_data/" +
        #    self.target + "*" + stripped_date + "*.ascii"
        #)

        if len(filename) != 0:
            filename = filename[0]
            data = pd.read_csv(
                filename,
                sep=" |:",
                header=0,
                names=["wave", "init", "O2", "H2O"],
                engine="python"
            )

            self.corrected_wave = data["wave"].to_numpy()
            self.flux_initial = data["init"].to_numpy()
            self.flux_corrO2 = data["O2"].to_numpy()
            self.flux_corrO2_h2O = data["H2O"].to_numpy()

        else:
            print('no corrected spectra available')

    def getSpectrum(self, xmin, xmax):
        """Function to update the wavelength region held in an EdiblesSpectrum object.

        Args:
            xmin (float): Minimum wavelength
            xmax (float): Maximum wavelength

        """
        assert xmin < xmax, "xmin must be less than xmax"
        assert xmin > np.nanmin(self.raw_wave), "xmin outside bounds"
        assert xmax < np.nanmax(self.raw_wave), "xmax outside bounds"

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

        try:
            # Sky transmission data
            sky_idx = np.where(np.logical_and(self.raw_sky_wave > xmin, self.raw_sky_wave < xmax))
            self.sky_wave = self.raw_sky_wave[sky_idx]
            self.sky_flux = self.raw_sky_flux[sky_idx]
        except AttributeError:
            pass

        self._interpolate(initial=True)

    def _interpolate(self, initial=False):
        '''Interpolation function used in shift().

        Warning:
            This function has not been tested for use on its own.
            It is implemented in getSpectrum and shift.

        '''

        # xmin = np.max([np.min(self.wave), np.min(self.bary_wave)])
        # xmax = np.min([np.max(self.wave), np.max(self.bary_wave)])

        grid_idx = np.where(np.logical_and(self.raw_grid > self.xmin,
                                           self.raw_grid < self.xmax))
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
            bfunc = interp1d(self.bary_wave, self.flux)
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
        self.wave = np.add(self.wave, shift)

        # Input checking
        assert (zoom_xmin > self.xmin) and (
            zoom_xmin > np.min(self.wave)
        ), 'zoom_xmin must be greater than ' + str(np.max([self.xmin, np.min(self.wave)]))
        assert (zoom_xmax < self.xmax) and (
            zoom_xmax < np.max(self.wave)
        ), 'zoom_xmax must be less than ' + str(np.min([self.xmax, np.max(self.wave)]))
        if isinstance(shift, np.ndarray):
            assert len(self.wave) == len(shift), 'Length of shift not equal to length of wave'

        self.xmin = zoom_xmin
        self.xmax = zoom_xmax

        self.bary_wave = self.wave + (self.v_bary / cst.c.to("km/s").value) * self.wave

        self._interpolate()

        # Zoom
        b_idx = np.where(np.logical_and(self.bary_wave > zoom_xmin, self.bary_wave < zoom_xmax))
        self.bary_wave = self.bary_wave[b_idx]
        self.bary_flux = self.flux[b_idx]

        t_idx = np.where(np.logical_and(self.wave > zoom_xmin, self.wave < zoom_xmax))
        self.wave = self.wave[t_idx]
        self.flux = self.flux[t_idx]

        try:
            sky_idx = np.where(np.logical_and(self.raw_sky_wave > zoom_xmin,
                                              self.raw_sky_wave < zoom_xmax))
            self.sky_wave = self.raw_sky_wave[sky_idx]
            self.sky_flux = self.raw_sky_flux[sky_idx]
        except AttributeError:
            pass

    def get_extra_data(self):
        '''Get sky transmission and atmosphere corrected data about EdiblesSpectrum if
            fully_featured is False

        '''

        assert not self.fully_featured, 'EdiblesSpectrum is already fully_featured!'

        self._sky_transmission()
        self._corrected_spectrum()
        self.fully_featured = True


def measure_snr(wave, flux, block_size=1.0, do_plot=False):
    """
    Estimate SNR of given spectral data
    :param wave: wavelength grid
    :type wave: ndarray
    :param flux: flux
    :type flux: ndarray
    :param do_plot: if set, make SNR plot
    :type  do_plot: bool

    :return: SNR, LAM, SNR and median flux of each of the 1A block
    :rtype: list
    """
    # split in blocks of given size
    xmin = wave[0]
    xmax = xmin + block_size
    SNR, LAM = [], []
    while xmin < wave[-1]:
        flux_block = flux[np.where((wave > xmin) & (wave < xmax))]
        if len(flux_block) == 1:
            break
        if (np.nanmean(flux_block) > 0.0):
            sigma_block = np.nanmean(flux_block) / np.nanstd(flux_block)
            SNR.append(sigma_block)
            LAM.append(xmin + (xmax - xmin) / 2.0)
        xmin = xmax.copy()
        xmax = xmin + block_size
    if (do_plot == True):
        plt.plot(LAM, SNR)
        plt.plot(LAM, np.convolve(SNR, np.ones(10) / 10, mode='same'))
        plt.show()
    return SNR, LAM

if __name__ == "__main__":
    # filename = "/HD170740/RED_860/HD170740_w860_redl_20140915_U.fits"

    filename = "/HD23466/BLUE_346/HD23466_w346_blue_20180731_O11.fits"
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




    # sp.get_extra_data()




    sp.getSpectrum(xmin=3275, xmax=3305)
    plt.plot(sp.wave, sp.flux, label="Geocentric Subset")
    plt.plot(sp.bary_wave, sp.bary_flux, label="Barycentric Subset")
    plt.plot(sp.grid, sp.interp_flux, label='Geocentric Interpolation')
    plt.plot(sp.grid, sp.interp_bary_flux, label='Barycentric Interpolation')
    if sp.fully_featured:
        plt.plot(sp.sky_wave, sp.sky_flux, 'k')
    plt.title('Data and Interpolations')
    plt.xlabel(r'Wavelength ($\AA$)')
    plt.ylabel('Flux')
    plt.legend()
    plt.show()

    plt.plot(sp.wave, sp.flux, label="Geocentric")
    plt.plot(sp.bary_wave, sp.bary_flux, label="Barycentric")
    shift = 0.05
    sp.shift(shift=shift, zoom_xmin=3285, zoom_xmax=3295)
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



    sp.get_extra_data()




    if sp.fully_featured:
        plt.plot(sp.corrected_wave, sp.flux_initial, 'r', label='Initial Flux')
        plt.plot(sp.corrected_wave, sp.flux_corrO2, 'g--', label='Corrected Flux O2')
        plt.plot(sp.corrected_wave, sp.flux_corrO2_h2O, 'b--', label='Corrected Flux H2O and O2')
        plt.legend(fontsize='small')
        plt.xlabel(r'Wavelength ($\AA$)')
        plt.ylabel('Flux')
        plt.title('Data and Telluric Data')
        plt.show()
