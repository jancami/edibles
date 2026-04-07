from edibles.utils.edibles_oracle import EdiblesOracle
from edibles.projects.edr5_integration import dr5_io, util_functions
from edibles import DATADIR
import matplotlib.pyplot as plt
from edibles.utils.voigt_profile import voigt_optical_depth
from lmfit import Model
import numpy as np
from importlib.resources import files
import pandas as pd

pythia = EdiblesOracle()
file_list = pythia.getFilteredObsList(object=['HD 186841'], MergedOnly=True, Wave=6707)

atomic_line_file = files('edibles') / 'data/auxiliary_data/line_catalogs/edibles_linelist_atoms.csv'
atomic_line_list = pd.read_csv(atomic_line_file)
print(atomic_line_list)


print(file_list)

elem_ind = 15
elem_row = atomic_line_list.loc[elem_ind]

print(elem_row)

elem_wave = elem_row['WavelengthAir']
my_f = elem_row['OscillatorStrength']
my_gamma = elem_row['Gamma']

order = 13

wave_range = [7696, 7720]

def voigt_slope(x, cont, slope, lambda0, b, n, v_rad):
    return (cont + slope * x) * np.exp(-voigt_optical_depth(x, lambda0=lambda0, b=b, N=n, f=my_f, gamma=my_gamma, v_rad=v_rad))


vmodel = Model(voigt_slope)
params = vmodel.make_params()

params['lambda0'].value = elem_wave
params['n'].min = 0
params['b'].min = 0
params['n'].value = 1e13
params['b'].value = 0.1
params['slope'].value = 0
params['v_rad'].value = 0


for file in file_list:
    spec = dr5_io.read_combined_spec(DATADIR / file)

    spec = spec[:, spec[4]==order]
    spec = util_functions.crop_spectrum(spec, *wave_range)

    fit_flux = spec[6]

    params['cont'].value = np.nanmedian(fit_flux)
    result = vmodel.fit(fit_flux, params, x=spec[0])

    print(spec)
    plt.plot(spec[0], fit_flux)
    # plt.plot(spec[0], spec[4])
    plt.plot(spec[0], result.best_fit)

    plt.show()
