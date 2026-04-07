from edibles.utils.edibles_oracle import EdiblesOracle
from edibles.projects.edr5_integration import dr5_io, util_functions, transformations
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

rv_range = 30
wave_ranges = [[4043, 4045], [7697, 7701]]
wave_ranges = np.array([transformations.doppler_shift_wl(atomic_line_list.loc[elem_inds, 'WavelengthAir'], -rv_range), 
                transformations.doppler_shift_wl(atomic_line_list.loc[elem_inds, 'WavelengthAir'], rv_range)]).T

print('wave_ranges', wave_ranges)

def make_multi_comp_voigt(n_components: int, wave_ranges: list) -> Callable:
    """
    Generate a multi-component Voigt profile fitting function compatible with lmfit.Model.
    
    Parameters
    ----------
    n_components : int
        Number of spectral components/windows.
    wave_ranges : list of tuples
        List of (min, max) wavelength ranges for each component window.
    
    Returns
    -------
    Callable
        A function with an explicit signature:
        f(x, b, n, v_rad, cont1, slope1, lambda01, f1, gamma1, cont2, ...)
    """
    if len(wave_ranges) != n_components:
        raise ValueError(f"wave_ranges must have {n_components} entries, got {len(wave_ranges)}")

    # Build explicit parameter name list
    comp_param_names = []
    for i in range(1, n_components + 1):
        comp_param_names += [f"cont{i}", f"slope{i}", f"lambda0{i}", f"f{i}", f"gamma{i}"]

    all_params = ["x", "b", "n", "v_rad"] + comp_param_names
    signature_str = ", ".join(all_params)
    func_name = f"voigt_{n_components}comp"

    # Build the function body
    body_lines = [f"def {func_name}({signature_str}):"]
    body_lines.append("    segments = []")

    params_per_comp = 5  # cont, slope, lambda0, f, gamma
    for i in range(n_components):
        idx = i + 1
        cont    = f"cont{idx}"
        slope   = f"slope{idx}"
        lambda0 = f"lambda0{idx}"
        f_      = f"f{idx}"
        gamma   = f"gamma{idx}"

        lo, hi = wave_ranges[i]

        if i == 0:
            body_lines.append(f"    x{idx} = x[x <= {hi}]")
        elif i == n_components - 1:
            body_lines.append(f"    x{idx} = x[x >= {lo}]")
        else:
            body_lines.append(f"    x{idx} = x[(x >= {lo}) & (x <= {hi})]")

        body_lines.append(
            f"    y{idx} = voigt_slope(x{idx}, {cont}, {slope}, {lambda0}, b, n, {f_}, {gamma}, v_rad)"
        )
        body_lines.append(f"    segments.append(y{idx})")

    body_lines.append("    return np.concatenate(segments)")

    func_code = "\n".join(body_lines)

    # Execute in a namespace that includes dependencies
    namespace = {"np": np, "voigt_slope": voigt_slope}
    exec(func_code, namespace)
    func = namespace[func_name]

    func.__doc__ = (
        f"Auto-generated {n_components}-component Voigt profile function.\n\n"
        f"Parameters: {signature_str}"
    )

    return func

generated_fit_function = make_multi_comp_voigt(len(elem_inds), wave_ranges)
    

vmodel = Model(generated_fit_function)
params = vmodel.make_params()

# Fixed atomic parameters — generalized over all components
for i, ind in enumerate(elem_inds, start=1):
    params[f'lambda0{i}'].set(value=elem_wave[ind],  vary=False)
    params[f'f{i}'].set(    value=my_f[ind],         vary=False)
    params[f'gamma{i}'].set(value=my_gamma[ind],     vary=False)
    params[f'slope{i}'].set(value=0)

# Shared parameters
params['b'].set(    value=0.1,  min=0)
params['n'].set(    value=1e11, min=0)
params['v_rad'].set(value=0,    min=-20, max=20)


# obs_times = obs_log['DateObs'].unique()

# for obs_time in obs_times:
#     filtered_obs = obs_log.loc[obs_log['DateObs'] == obs_time]
#     filtered_obs = filtered_obs[filtered_obs['Order'] == 'ALL']
#     print(filtered_obs)

# old single cloud sight lines
scl_old = ['HD 23180', 'HD 24398', 'HD 144470', 'HD 147165', 'HD 147683', 'HD 149757', 'HD 166937', 'HD 170740',
           'HD 184915', 'HD 185418', 'HD 185859', 'HD 203532']

for star_name in scl_old:

    # Gather file lists for each wavelength window
    file_lists = []
    for wave_range in wave_ranges:
        pythia    = EdiblesOracle()
        file_list = pythia.getFilteredObsList(object=[star_name], MergedOnly=True, Wave=np.mean(wave_range))
        file_lists.append(file_list)

    file_lists = np.array(file_lists).T  # shape: (n_observations, n_components)

    for obs_files in file_lists:  # one row = one observation epoch
        # Load and crop each component spectrum
        spectra = []
        for i, (fpath, order, wave_range) in enumerate(zip(obs_files, orders, wave_ranges)):
            spec = dr5_io.read_combined_spec(DATADIR / fpath, bary_corr=True)
            spec = spec[:, spec[4] == order]
            spec = util_functions.crop_spectrum(spec, *wave_range)
            spectra.append(spec)

        # Stitch into a single spectrum for fitting
        fit_wave = np.concatenate([s[0] for s in spectra])
        fit_flux = np.concatenate([s[1] for s in spectra])

        plt.plot(fit_wave, fit_flux, label=star_name)
        plt.show()

        # Set continuum initial guesses per component
        for i, spec in enumerate(spectra, start=1):
            params[f'cont{i}'].set(value=np.nanmedian(spec[1]))

        result = vmodel.fit(fit_flux, params, x=fit_wave)

        plt.plot(fit_flux,        label='data')
        plt.plot(result.best_fit, label='fit')
        plt.legend()
        plt.show()