from edibles.fit.models.create_model import createLine, createCont
import astropy.constants as cst

from edibles.fit.fit import fit
from edibles.edibles.functions.file_search import FilterDR
from edibles.edibles.functions.edibles_spectrum import EdiblesSpectrum
from edibles.edibles.functions.peak_wavelength import largest_peak_wavelength, all_prominent_peak_wavelength
from astropy.constants import c

import matplotlib.pyplot as plt
import numpy as np


def line_parameters(peak_wavelength):
    return (peak_wavelength, 2.34333, 0.00328429,
            6.61966e+13 * 0.00406 * 1e8 * np.pi * (cst.e.esu.value) ** 2 * (1e-8 * peak_wavelength) ** 2 / (
                    cst.m_e.to('g').value * (cst.c.to('cm/s').value) ** 2))


def get_parameters(fit_model):
    parts = fit_model.parts
    parameter_data = []

    while True:
        absorption_line = parts[1]
        parameter_data.append({
            'lam_0': absorption_line.lam_0.val,
            'b': absorption_line.b.val,
            'd': absorption_line.d.val,
            # gamma value doesn't align with nist data
            'gamma': 4 * np.pi * c.to('AA/s').value * absorption_line.d.val /
                     absorption_line.lam_0.val ** 2,
            'tau_0': absorption_line.tau_0.val
        })
        try:
            parts = parts[0].parts
        except AttributeError:
            break
    return parameter_data


if __name__ == '__main__':
    star = 'HD170740'
    wav_range = (3301, 3303.5)
    wavelength = sum(wav_range) / 2
    date = FilterDR().filterAll(star=star, wavelength=wavelength).getDates()[0]
    print(FilterDR().filterAll(star=star, wavelength=wavelength).getDates())

    # for date in FilterDR().filterAll(star=star, wavelength=wavelength).getDates()[:1]:
    filename = FilterDR().filterAll(star=star, date=date, wavelength=wavelength).getAllFileNames()[0]
    wave, flux = EdiblesSpectrum(filename).getSpectrum(*wav_range)

    # peak_wavelength = largest_peak_wavelength(wave, flux, n=12)
    peak_wavelength = all_prominent_peak_wavelength(wave, flux, 0.2)
    print(peak_wavelength)

    n_points = 2 * len(peak_wavelength) + 1
    model = createCont((wave, flux), n_points)
    lines = [createLine("peak_{:.2f}".format(wav), *line_parameters(wav)) for wav in peak_wavelength]
    for line in lines:
        model *= line

    fit_model = fit(star, (wave, flux), model)

    # calculate gamma for all peaks
    data = get_parameters(fit_model)
    print(data)
