import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from edibles_settings import *

class EdiblesSpectrum:
# This object will contain a spectrum from EDIBLES, and a set of methods to operate on the data. 

	def load_spectrum (self):
		# Assume the file is a DR3 product here. 
		hdu = fits.open(self.filename)
		self.header = hdu[0].header
		self.flux = hdu[0].data
		self.flux_units="arbitrary"
		crval1 = hdu[0].header["CRVAL1"]
		cdelt1 = hdu[0].header["CDELT1"]
		nwave = len(self.flux)
		grid = np.arange(0, nwave, 1)
		self.wave = (grid) * cdelt1 + crval1
		self.wave_units = "AA"
		self.reference_frame = "geocentric"
		self.v_bary = hdu[0].header["HIERARCH ESO QC VRAD BARYCOR"]

	def __init__(self, filename):
		self.filename = filename
		self.load_spectrum()

	def GetSpectrum (self):
		return self.wave, self.flux

    def getSpectrumRange(self, xmin, xmax):

        assert xmin < xmax, 'xmin must be less than xmax'


        # create data subset
        min_idx = (np.abs(self.wave - xmin)).argmin()
        max_idx = (np.abs(self.wave - xmax)).argmin()
        wave_range = self.wave[min_idx:max_idx]
        flux_range = self.flux[min_idx:max_idx]

        data = (wave_range, flux_range)

        return data

if __name__ == '__main__':
    sp = EdiblesSpectrum(datadir+"/HD170740/RED_860/HD170740_w860_n20_20140916_L.fits")
    print("Barycentric Velocity is", sp.v_bary)
    wave,flux = sp.getSpectrum()
    plt.plot(wave, flux)
    axes = plt.gca()
    axes.set_xlim([7660,7705])
    axes.set_ylim([0,160])
    plt.vlines((7667.021,7701.093), 0, 160, linestyles='dashed', colors='r')
    plt.show()

    wave_range, flux_range = sp.getSpectrumRange(7660,7705)
    plt.plot(wave_range, flux_range)
    plt.show()

