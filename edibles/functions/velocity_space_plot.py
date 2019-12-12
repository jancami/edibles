import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from scipy.signal import find_peaks
from heapq import nsmallest
from edibles.edibles.functions.edibles_spectrum import EdiblesSpectrum
from edibles.edibles import DATADIR
from edibles.edibles.functions.atomic_line_tool import AtomicLines
from edibles.edibles.functions.file_search import FilterDR
from astropy.constants import c
import bisect

from edibles.edibles.functions.functions import peak_wavelength_largest

def spliceToRange(l, minimum, maximum):
    def find_gt(a, x):
        'Find leftmost value greater than x'
        return bisect.bisect_right(a, x)

    def find_lt(a, x):
        'Find rightmost value less than x'
        return bisect.bisect_left(a, x)

    return find_gt(l, minimum), find_lt(l, maximum)


def velocity_space(star, time, ion):
    """

    :param star: star name
    :param time:
    :param lab_wavelength: Array of wavelengths
    """

    lab_wavelength = AtomicLines().getAllLabWavelength(ion)

    data = FilterDR().filterAll(star=star, date=time, wavelength=lab_wavelength, order=[4, 11])
    print(data)
    df = data.getDataFrame()

    for index, row in df.iterrows():
        # get wavelength and flux values from file
        xrange = [row.WaveMin, row.WaveMax]
        wavelengths = sorted([wavelength for wavelength in lab_wavelength if xrange[0] < wavelength < xrange[1]])
        wave_dist = 2 * (wavelengths[-1] - wavelengths[0])  # arbitrary parameter to limit range of spectrum
        sp = EdiblesSpectrum(row.Filename)
        wave, flux = sp.getSpectrum(wavelengths[0] - wave_dist, wavelengths[-1] + wave_dist)

        plt.figure('Figure ' + str(index))
        plt.plot(wave, flux)
        plt.title(star + ' at ' + time[:4] + '-' + time[4:6] + '-' + time[6:])
        plt.xlabel('Wavelength (AA)')
        plt.ylabel('Flux')
        ax = plt.gca()
        plot_lines(ax, ion, wavelengths, flux)

        # TODO: doesn't allow for size 1 AND size > 2
        # TODO: doesn't allow for non-adjacent peaks
        peak_wavelength = peak_wavelength_largest(wave, flux)
        # bisector = sum(peak_wavelength) / len(peak_wavelength)
        peak_distance = abs(peak_wavelength[1] - peak_wavelength[0])

        for peakn in range(len(wavelengths)):
            def transform(wavelength):
                return (wavelengths[peakn] - wavelength) / wavelengths[peakn] * c.to('km/s').value

            v = transform(wave)

            av_v = transform(peak_wavelength[peakn])
            approx_v_range = peak_distance / wavelengths[peakn] * c.to('km/s').value

            v_min, v_max = spliceToRange(v, av_v - 0.5 * approx_v_range, av_v + 0.5 * approx_v_range)
            median_flux = flux[(v_max - v_min) // 2]

            plt.figure('Velocity Space')
            plt.plot(v, flux / median_flux, label='$\lambda$=' + str(wavelengths[peakn]))
            plt.xlim([av_v - 0.5 * approx_v_range, av_v + 0.5 * approx_v_range])

        plt.title(star + ' at ' + time[:4] + '-' + time[4:6] + '-' + time[6:])
        plt.xlabel('Velocity (km/s)')
        plt.ylabel('Flux')
        plt.legend()



def plot_lines(ax, ion, wavelengths, flux):
    """

    :type ax: matplotlib.axis
    """
    ymax = max(flux)
    ax.vlines(wavelengths, 0.95 * min(flux), 0.98 * ymax)
    for wavelength in wavelengths:
        ax.text(wavelength, ymax, ion, horizontalalignment='center')


if __name__ == '__main__':
    star = 'HD170740'
    ion = 'Na I'

    filter = FilterDR().filterAll(star=star, wavelength=AtomicLines().getAllLabWavelength(ion))
    time_list = filter.getDates()
    print(time_list)

    velocity_space(star, time_list[0], ion)
    plt.show()