import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from scipy.signal import find_peaks
from heapq import nsmallest
from edibles.edibles_spectrum import EdiblesSpectrum
from edibles.edibles_settings import *
from edibles.functions.find_f_known import AtomicLines


def wavelength_at_peaks(wave, flux, n=2):
    """
    Returns the wavelengths of the lowest n peaks of the equation
    """
    peaks, _ = find_peaks(-flux)  # indices of peaks
    peak_flux = nsmallest(2, flux[peaks])  # smallest two flux values at peaks
    peak_wavelength = [wave[np.where(flux == peak)][0] for peak in peak_flux]  # corresponding wavelengths
    return peak_wavelength


def filter_dataframe(df, star, time, lab_wavelength):
    df.DateObs = df.DateObs.apply(parse_time)  # format date
    df2 = df[(df.DateObs == time) & (df.Object == star)]  # get correct date

    # get only in range
    wave1 = lab_wavelength[0]
    filt = (df2.WaveMax > wave1) & (df2.WaveMin < wave1)
    for wave in lab_wavelength[1:]:
        filt = filt | ((df2.WaveMax > wave) & (df2.WaveMin < wave))
    df2 = df2[filt]

    return df2[~df2.Filename.str.contains('O')]  # get only the largest range


def parse_time(timestamp):
    """
    :param timestamp: ex. 2014-10-29T07:01:33.557
    :return: 20141029
    """
    return timestamp.split('T')[0].replace('-', '')


def velocity_space(star, time, lab_wavelength):
    """

    :param star: star name
    :param time:
    :param lab_wavelength: Array of wavelengths
    :return:
    """
    c = 2.9979 * 10 ** 5  # km/s

    df = pd.read_csv('./data/DR3_ObsLog.csv')
    df = filter_dataframe(df, star, time, lab_wavelength)
    # pd.set_option('display.max_colwidth', -1)
    # print(df.to_string())
    df = df.reset_index(drop=True)

    for index, row in df.iterrows():
        # get wavelength and flux values from file
        xrange = [row.WaveMin, row.WaveMax]
        wavelengths = sorted([wavelength for wavelength in lab_wavelength if xrange[0] < wavelength < xrange[1]])
        wave_dist = 2 * (wavelengths[-1] - wavelengths[0])  # arbitrary parameter to limit range of spectrum
        sp = EdiblesSpectrum(datadir + row.Filename)
        wave, flux = sp.getSpectrum(wavelengths[0] - wave_dist, wavelengths[-1] + wave_dist)

        plt.figure()
        plt.plot(wave, flux)
        plt.vlines(wavelengths, min(flux), max(flux))

        # TODO: doesn't allow for size 1 AND size > 2
        # TODO: doesn't allow for non-adjacent peaks
        peak_wavelength = sorted(wavelength_at_peaks(wave, flux))
        # bisector = sum(peak_wavelength) / len(peak_wavelength)
        peak_distance = abs(peak_wavelength[1] - peak_wavelength[0])

        for peakn in range(len(wavelengths)):
            def transform(wavelength):
                return (wavelengths[peakn] - wavelength) / wavelengths[peakn] * c

            v = transform(wave)

            # TODO: normalize flux values before plotting
            plt.figure('Velocity Space')
            plt.plot(v, flux)
            av_v = transform(peak_wavelength[peakn])
            approx_v_range = peak_distance / wavelengths[peakn] * c
            plt.xlim([av_v - 0.5 * approx_v_range, av_v + 0.5 * approx_v_range])
            plt.xlabel('Velocity (km/s)')
            plt.ylabel('Flux')

    plt.show()


if __name__ == '__main__':
    star = 'HD170740'
    time = '20140916'
    na_observed_wavelength = AtomicLines().getAllLabWavelength('Na I')

    velocity_space(star, time, na_observed_wavelength)
