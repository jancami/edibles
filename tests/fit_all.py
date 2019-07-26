from edibles.functions.voigt import voigtAbsorptionLine, voigtOpticalDepth, voigtNorm
from edibles.fit.models.create_model import createLine, createCont
import astropy.constants as cst

from edibles.fit.fit import fit
from edibles.functions.file_search import FilterDR
from edibles.functions.edibles_spectrum import EdiblesSpectrum
from edibles.functions.peak_wavelength import largest_peak_wavelength, all_prominent_peak_wavelength

import matplotlib.pyplot as plt
import numpy as np


def line_parameters(peak_wavelength):
    return (peak_wavelength, 2.34333, 0.00328429, 21.8628 * 0.00406 / (
            cst.m_e.to('g').value * (cst.c.to('cm/s').value) ** 2 / (
            np.pi * (cst.e.esu.value) ** 2 * 1e8 * peak_wavelength)))


if __name__ == '__main__':
    star = 'HD170740'
    wav_range = (7652.5,7668)
    wavelength = sum(wav_range) / 2
    date = FilterDR().filterAll(star=star, wavelength=wavelength).getDates()[0]
    print(FilterDR().filterAll(star=star, wavelength=wavelength).getDates())

    # for date in FilterDR().filterAll(star=star, wavelength=wavelength).getDates()[:1]:
    filename = FilterDR().filterAll(star=star, date=date, wavelength=wavelength).getAllFileNames()[0]
    wave, flux = EdiblesSpectrum(filename).getSpectrum(*wav_range)

    # peak_wavelength = largest_peak_wavelength(wave, flux, n=12)
    peak_wavelength = all_prominent_peak_wavelength(wave, flux, 0.05)
    print(peak_wavelength)

    n_points = 2 * len(peak_wavelength) + 1
    model = createCont((wave, flux), n_points)
    lines = [createLine("MODEL", *line_parameters(wav)) for wav in peak_wavelength]
    for line in lines:
        model *= line

    fit_model = fit(star, (wave, flux), model)
