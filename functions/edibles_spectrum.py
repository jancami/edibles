import numpy as np
from astropy.io import fits
import astropy.constants as cst
import matplotlib.pyplot as plt
from edibles.edibles_settings import *

class EdiblesSpectrum:
# This object will contain a spectrum from EDIBLES, and a set of methods to operate on the data. 

    def loadSpectrum (self):
        # Assume the file is a DR3 product here. 
        hdu = fits.open(self.filename)
        self.header = hdu[0].header
        self.date = self.header["DATE-OBS"]
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
        self.bary_wave = self.wave + (self.v_bary/cst.c.to('km/s').value)*self.wave

    def __init__(self, filename):
        """
        Filename is relative to the DR3 directory
        """
        self.filename = datadir + filename
        self.loadSpectrum()

    def getSpectrum(self, xmin=None, xmax=None):

        if (xmin is not None) and (xmax is not None):

            assert xmin < xmax, 'xmin must be less than xmax'
            idx = (self.wave > xmin) * (self.wave < xmax)

            return self.wave[np.where(idx)], self.flux[np.where(idx)]
        return self.wave, self.flux


if __name__ == '__main__':
    filename = '/HD170740/RED_860/HD170740_w860_n20_20140916_L.fits'
    sp = EdiblesSpectrum(filename)
    print("Barycentric Velocity is", sp.v_bary)
    plt.plot(sp.wave, sp.flux, label='Geocentric')
    plt.plot(sp.bary_wave, sp.flux, label='Barycentric')
    axes = plt.gca()
    axes.set_xlim([7660,7705])
    axes.set_ylim([0,160])
    plt.vlines((7667.021,7701.093), 0, 160, linestyles='dashed', colors='r')
    plt.legend()
    plt.show()

