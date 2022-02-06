## This script is an example on measuring DIB5780 in HD73882

# import statements
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Qt5Agg', force=True)
import os

# EDIBLES related
from edibles.utils.edibles_oracle import EdiblesOracle  # file selection
from edibles.utils.edibles_spectrum import EdiblesSpectrum  # data I/O
from edibles.utils.ContinuumFitter import ContinuumFitter # data normalization
from edibles.utils.ISLineFitter import measure_snr


# define a class that will measure DIB5780 from normalized data
class DIB5780():
    def __init__(self, wave, flux, snr):
        self.wave = wave
        self.flux = flux
        self.snr = snr
        self.cont = 1 - 1 / snr

        # properties that I'm interested in
        self.wave_idx = np.where(flux == np.min(flux))[0][0]
        self.CD = 1 - self.flux[self.wave_idx]
        self.HalfMax = 1 - 0.5 * self.CD
        self.wavelength = wave[self.wave_idx]

        self.EW = 0
        self.err = None
        self.FWHM_boundaries = [None, self.wavelength, None]
        self.Int_boundaries = [None, self.wavelength, None]

        self.__integration(1)  # RHS
        self.__integration(-1)  # LHS

        self.EW = self.EW * 1000  # convert to mA
        self.err = 1000 * 1.064 * (self.FWHM_boundaries[-1] - self.FWHM_boundaries[0]) / self.snr

    def __integration(self, direction):
        # direction = 1 --> right hand side
        # direction = -1 --> left hand side
        idx_now = self.wave_idx + direction
        while True:
            # hit continuum if next five points are within 1 / 1/SNR
            if (self.flux[idx_now+direction*4:idx_now-direction:-direction] >= self.cont).all():
                self.Int_boundaries[direction+1] = self.wave[idx_now]
                break

            # integration
            dx = np.abs(self.wave[idx_now] - self.wave[idx_now - direction])
            dy = 2 - self.flux[idx_now] - self.flux[idx_now - direction]
            self.EW = self.EW + 0.5 * dx * dy

            # FWHM, only write once
            if self.flux[idx_now] > self.HalfMax and self.FWHM_boundaries[direction+1] is None:
                self.FWHM_boundaries[direction+1] = self.wave[idx_now]

            # Move on
            idx_now = idx_now + direction



# set upt rootpath
RootPath = "/Users/haoyufan/DonOutput"

# file selection and import data
fits_log = EdiblesOracle()
file_all = fits_log.getFilteredObsList(object=["HD 73882"],
                                       OrdersOnly=True,
                                       Wave=5780)

# Two orders available, keep both for cross reference
for filename in file_all:
    order = filename.split("_")[-1].split(".")[0]
    if not os.path.exists(os.path.join(RootPath, order)):
        os.mkdir(os.path.join(RootPath, order))

    sp = EdiblesSpectrum(filename)
    wave, flux = sp.bary_wave, sp.flux
    idx = (wave > 5775) & (wave < 5788)
    wave, flux = wave[idx], flux[idx]
    SNR = np.max(measure_snr(wave, flux))

    # Fitting Continuum using ContinuumFitter
    normalizer = ContinuumFitter(wave=wave, flux=flux)
    cont, anchor = normalizer.SplineManualRegion(n_anchors=4, n_regions=99)

    # Graph output 1, un normalized data
    plt.plot(wave, flux)
    plt.plot(wave, cont(wave))
    plt.fill_between(wave,
                     cont(wave) * (1 + 1/SNR),
                     cont(wave) * (1 - 1/SNR),
                     color="C1",
                     alpha=0.2)
    plt.scatter(anchor.T[0], anchor.T[1], marker="x", s=80, color="r")
    plt.grid()
    plt.xlabel(filename)
    plt.ylabel("Raw Data, Spline Continuum, and Uncertainties")
    plt.savefig(os.path.join(RootPath, order, order + "_Normalization_1.png"))
    plt.close()

    # Normalization and Grapth output 2
    flux = flux / cont(wave)
    SNR = np.max(measure_snr(wave, flux))

    plt.plot(wave, flux, label="SNR = %.2f" % SNR)
    plt.plot(wave, np.ones_like(wave))
    plt.fill_between(wave,
                     np.ones_like(wave) * (1 + 1/SNR),
                     np.ones_like(wave) * (1 - 1/SNR),
                     color="C1", alpha=0.2)
    plt.scatter(anchor.T[0], np.ones_like(anchor.T[0]), marker="x", s=80, color="r")
    plt.grid()
    plt.xlabel(filename)
    plt.ylabel("Normalized Data and SNR")
    plt.legend(loc="lower right")
    plt.savefig(os.path.join(RootPath, order, order + "_Normalization_2.png"))
    plt.close()

    # Get DIB measurement and output result
    measurer = DIB5780(wave, flux, SNR)
    plt.plot(wave, flux, label="EW = %.2f pm %.2f mAA" % (measurer.EW, measurer.err))
    plt.plot(wave, np.ones_like(wave))
    plt.fill_between(wave,
                     np.ones_like(wave) * (1 + 1 / SNR),
                     np.ones_like(wave) * (1 - 1 / SNR),
                     color="C1", alpha=0.2)
    plt.plot([measurer.Int_boundaries[0]] * 2, [0.995, 1.005], color="r")
    plt.plot([measurer.Int_boundaries[-1]] * 2, [0.995, 1.005], color="r")
    plt.plot([measurer.wavelength] * 2,
             [measurer.HalfMax + 0.5 * measurer.CD, measurer.HalfMax - 0.5 * measurer.CD],
             color="0.5")
    plt.plot(measurer.FWHM_boundaries[0::2],
             [measurer.HalfMax]*2,
             color="0.5",
             label="FWHM = %.2f AA" % (measurer.FWHM_boundaries[-1] - measurer.FWHM_boundaries[0]))
    plt.legend(loc="lower right")
    plt.xlabel(filename)
    plt.ylabel("Measurement")
    plt.grid()
    plt.savefig(os.path.join(RootPath, order, order + "_Measurement.png"))
    plt.close()

