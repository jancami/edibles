from scipy.signal import find_peaks
from heapq import nsmallest

import numpy as np


def largest_peak_wavelength(wave, flux, n=2):
    """
    Returns the wavelengths of the lowest n peaks of the equation
    """
    peaks, _ = find_peaks(-flux)  # indices of peaks
    peak_flux = nsmallest(n, flux[peaks])  # smallest two flux values at peaks
    peak_wavelength = [wave[np.where(flux == peak)][0] for peak in peak_flux]  # corresponding wavelengths
    return peak_wavelength


def all_prominent_peak_wavelength(wave, flux, prominence):
    yrange = max(flux) - min(flux)
    peaks, _ = find_peaks(-flux, prominence=prominence * yrange)
    return [wave[peak] for peak in peaks]
