from edibles.utils.edibles_oracle import EdiblesOracle
from edibles.projects.edr5_integration import dr5_io, util_functions
from edibles import DATADIR
import matplotlib.pyplot as plt
from edibles.utils.voigt_profile import voigt_optical_depth
from lmfit import Model
import numpy as np
from importlib.resources import files
import pandas as pd
from typing import Callable

atomic_line_file = files('edibles') / 'data/auxiliary_data/line_catalogs/edibles_linelist_atoms.csv'
atomic_line_list = pd.read_csv(atomic_line_file)
print(atomic_line_list)

obs_log = pd.read_csv(files('edibles') / 'data/DR5_ObsLog.csv')

# KI
elem_inds = [12, 15]
elem_row = atomic_line_list.loc[elem_inds]
order1 = 10
order2 = 13

# NaI
elem_inds = [20, 21, 22, 23]
elem_row = atomic_line_list.loc[elem_inds]
orders = [12, 12, 5, 5]

print(elem_row)

elem_wave = elem_row['WavelengthAir']
my_f = elem_row['OscillatorStrength']
my_gamma = elem_row['Gamma']

print(my_f)


def voigt_slope(x, cont, slope, lambda0, b, n, f, gamma, v_rad):
    return (cont + slope * x) * np.exp(-voigt_optical_depth(x, lambda0=lambda0, b=b, N=n, f=f, gamma=gamma, v_rad=v_rad))


wave_ranges = [[4043, 4045], [7697, 7701]]

def make_multi_comp_voigt(n_components: int, wave_ranges: list) -> Callable:
    for comp, wave_range in zip(n_components, wave_ranges):


def two_comp_voigt(x, b, n, v_rad, cont1, slope1, lambda01, f1, gamma1, cont2, slope2, lambda02, f2, gamma2):
    # separate windows
    x1 = x[x<=wave_ranges[0][1]]
    x2 = x[x>=wave_ranges[1][0]]
    y1 = voigt_slope(x1, cont1, slope1, lambda01, b, n, f1, gamma1, v_rad)
    y2 = voigt_slope(x2, cont2, slope2, lambda02, b, n, f2, gamma2, v_rad)
    return np.concatenate((y1, y2))
    

vmodel = Model(two_comp_voigt)
params = vmodel.make_params()

params['lambda01'].value = elem_wave[elem_inds[0]]
params['lambda01'].vary = False
params['lambda02'].value = elem_wave[elem_inds[1]]
params['lambda02'].vary = False
params['f1'].value = my_f[elem_inds[0]]
params['f1'].vary = False
params['gamma1'].value = my_gamma[elem_inds[0]]
params['gamma1'].vary = False
params['f2'].value = my_f[elem_inds[1]]
params['f2'].vary = False
params['gamma2'].value = my_gamma[elem_inds[1]]
params['gamma2'].vary = False
params['n'].min = 0
params['b'].min = 0
params['n'].value = 1e11
params['b'].value = 0.1
params['slope1'].value = 0
params['slope2'].value = 0
params['v_rad'].value = 0
params['v_rad'].min = -20
params['v_rad'].max = 20


# obs_times = obs_log['DateObs'].unique()

# for obs_time in obs_times:
#     filtered_obs = obs_log.loc[obs_log['DateObs'] == obs_time]
#     filtered_obs = filtered_obs[filtered_obs['Order'] == 'ALL']
#     print(filtered_obs)

# old single cloud sight lines
scl_old = ['HD 23180', 'HD 24398', 'HD 144470', 'HD 147165', 'HD 147683', 'HD 149757', 'HD 166937', 'HD 170740',
           'HD 184915', 'HD 185418', 'HD 185859', 'HD 203532']

for star_name in scl_old:
    file_lists = []
    for wave_range in wave_ranges:
        pythia = EdiblesOracle()
        file_list = pythia.getFilteredObsList(object=[star_name], MergedOnly=True, Wave=np.mean(wave_range))
        file_lists.append(file_list)

    file_lists = np.array(file_lists).T

    for file1, file2 in file_lists:
        spec1 = dr5_io.read_combined_spec(DATADIR / file1, bary_corr=True)
        spec2 = dr5_io.read_combined_spec(DATADIR / file2, bary_corr=True)

        spec1 = spec1[:, spec1[4]==order1]
        spec1 = util_functions.crop_spectrum(spec1, *wave_ranges[0])

        spec2 = spec2[:, spec2[4]==order2]
        spec2 = util_functions.crop_spectrum(spec2, *wave_ranges[1])

        fit_wave = np.concatenate((spec1[0], spec2[0]))
        fit_flux = np.concatenate((spec1[1], spec2[6]))

        fit_spec = np.array([fit_wave, fit_flux])

        plt.plot(fit_spec[0], fit_spec[1])
        plt.show()


        params['cont1'].value = np.nanmedian(spec1[1])
        params['cont2'].value = np.nanmedian(spec2[6])
        result = vmodel.fit(fit_flux, params, x=fit_spec[0])

        # print(fit_spec)
        # plt.plot(fit_spec[0], fit_flux)
        # plt.plot(fit_spec[0], result.best_fit)
        plt.plot(fit_flux)
        plt.plot(result.best_fit)

        plt.show()
