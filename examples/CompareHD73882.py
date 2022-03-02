import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Qt5Agg', force=True)
import os
import astropy.constants as cst

# EDIBLES related
from edibles.utils.edibles_oracle import EdiblesOracle  # file selection
from edibles.utils.edibles_spectrum import EdiblesSpectrum  # data I/O
from edibles.utils.ContinuumFitter import ContinuumFitter # data normalization
from edibles.utils.ISLineFitter import measure_snr

# 5780 measurer
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



# Target info
all_targets = [
    ["HD 73882", 20.3, "HD73882, O8.5 IV, Ebv=0.68, fH2=0.65"],
    ["HD 54662", 29.9, "HD54662, O7 V, Ebv=0.32, fH2=0.08"],
    ["HD 93222", -8.2, "HD93222, O7 V, Ebv=0.34, fH2=0.05"],
    ["HD 149038", -7.4, "HD149038, O9.7 Iab, Ebv=0.31, fH2=0.36"],
    #["HD 149757", -15.1, "HD149757, O9.2 IV, Ebv=0.32, fH2=0.63"],
]
order = "O3"
scale = False
fits_log = EdiblesOracle()


for i, target in enumerate(all_targets):
    all_files = fits_log.getFilteredObsList(object=[target[0]],
                                            OrdersOnly=True,
                                            Wave=5780)
    fits_filename = [item for item in all_files if order in item][0]

    sp = EdiblesSpectrum(fits_filename)
    wave, flux = sp.bary_wave, sp.flux
    wave = wave * (1 - target[1] / cst.c.to("km/s").value)
    idx = (wave > 5775) & (wave < 5788)
    wave, flux = wave[idx], flux[idx]

    normalizer = ContinuumFitter(wave=wave, flux=flux)
    cont, anchor = normalizer.SplineManualRegion(n_anchors=4, n_regions=99)
    flux = flux / cont(wave)
    SNR = np.max(measure_snr(wave, flux))

    measurer = DIB5780(wave, flux, SNR)

    if scale:
        if i == 0:
            EW2match = measurer.EW
            scaler = 1.0
        else:
            scaler = EW2match / measurer.EW
        flux = (flux - 1.0) * scaler + 1.0

    plt.plot(wave, flux,
             label=target[-1] + ", EW=%.2f mAA" % measurer.EW)

plt.grid()
plt.legend(loc="lower right")
plt.xlabel("Wavelength (AA)")
plt.ylabel("Normalized Flux")
plt.show()
