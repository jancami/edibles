import numpy as np
from astropy.io import fits
import astropy.constants as cst
import matplotlib.pyplot as plt
from edibles.edibles import DATADIR


class EdiblesSpectrum:
    # This object will contain a spectrum from EDIBLES,
    # and a set of methods to operate on the data.

    def loadSpectrum(self):
        # Assume the file is a DR3 product here.
        hdu = fits.open(self.filename)
        self.header = hdu[0].header
        self.target = self.header["OBJECT"]
        self.date = self.header["DATE-OBS"]
        self.flux = hdu[0].data
        self.flux_units = "arbitrary"
        crval1 = self.header["CRVAL1"]
        cdelt1 = self.header["CDELT1"]
        nwave = len(self.flux)
        grid = np.arange(0, nwave, 1)
        self.wave = (grid) * cdelt1 + crval1
        self.wave_units = "AA"
        self.reference_frame = "geocentric"
        self.v_bary = self.header["HIERARCH ESO QC VRAD BARYCOR"]
        self.bary_wave = self.wave + (self.v_bary/cst.c.to('km/s').value)*self.wave

    def __init__(self, filename):
        """
        Filename is relative to the DR3 directory
        """
        self.filename = DATADIR + filename
        self.loadSpectrum()

    def getSpectrum(self, xmin=None, xmax=None, bary=False):

        if bary is True:
            if (xmin is not None) and (xmax is not None):
                assert xmin < xmax, 'xmin must be less than xmax'
                idx = (self.bary_wave > xmin) * (self.bary_wave < xmax)

                return self.bary_wave[np.where(idx)], self.flux[np.where(idx)]

        else:
            if (xmin is not None) and (xmax is not None):
                assert xmin < xmax, 'xmin must be less than xmax'
                idx = (self.wave > xmin) * (self.wave < xmax)

                return self.wave[np.where(idx)], self.flux[np.where(idx)]

        return self.wave, self.flux


if __name__ == '__main__':
    filename = '/HD170740/RED_860/HD170740_w860_redl_20140915_O12.fits'
    sp = EdiblesSpectrum(filename)
    print("Barycentric Velocity is", sp.v_bary)
    print(sp.target)
    plt.plot(sp.wave, sp.flux, label='Geocentric')

    bary_data = sp.getSpectrum(xmin=7660, xmax=7705, bary=True)

    plt.plot(bary_data[0], bary_data[1], label='Barycentric')
    axes = plt.gca()
    # axes.set_xlim([7660, 7705])
    # axes.set_ylim([0, 160])
    # plt.vlines((7667.021, 7701.093), 0, 160, linestyles='dashed', colors='r')
    plt.legend()
    plt.show()
