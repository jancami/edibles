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

# NaI
elem_inds = [20, 21, 22, 23]
elem_df = atomic_line_list.loc[elem_inds]

print(elem_df)

range_list = [[3301, 3304], [5888, 5900]]

n_comp = 3

comp_list = []
fit_df = pd.DataFrame()

# Define cloud components
for c_comp in range(n_comp):
    i_df = elem_df.copy()
    i_df.loc[:, 'c_comp'] = c_comp
    fit_df = pd.concat((fit_df, i_df), ignore_index=True)

# Define wave range for each line
for i, row in fit_df.iterrows():
    w_min, w_max = [(low, high) for low, high in range_list if low <= row['WavelengthAir'] <= high][0]
    fit_df.loc[i, 'w_min'] = w_min
    fit_df.loc[i, 'w_max'] = w_max


print(fit_df)


def cont_slope(x, cont, slope):
    return (cont + slope * x)

def add_voigt(x, y, lambda0, b, n, f, gamma, v_rad):
    return y * np.exp(-voigt_optical_depth(x, lambda0=lambda0, b=b, N=n, f=f, gamma=gamma, v_rad=v_rad))


def voigt_slope(x, cont, slope, lambda0, b, n, f, gamma, v_rad):
    return (cont + slope * x) * np.exp(-voigt_optical_depth(x, lambda0=lambda0, b=b, N=n, f=f, gamma=gamma, v_rad=v_rad))


def make_multi_comp_voigt(df_list: pd.DataFrame) -> Callable:
    range_df = df_list[['w_min', 'w_max']].drop_duplicates().reset_index(drop=True)
    c_comps = df_list['c_comp'].drop_duplicates().reset_index(drop=True)

    comp_param_names = []
    for i, w_range in range_df.iterrows():
        comp_param_names += [f'cont_{i}', f'slope_{i}']
    for c in c_comps:
        comp_param_names += [f'b_{c}', f'n_{c}', f'v_rad_{c}']
    for j, row in df_list.iterrows():
        comp_param_names += [f'lambda0_{j}', f'f_{j}', f'gamma_{j}']

    all_params = ['x'] + comp_param_names
    signature_str = ", ".join(all_params)
    func_name   = f"voigt_n_comp"
    
    body_lines  = [f"def {func_name}({signature_str}):"]
    body_lines.append("    segments = []")

    for i, w_range in range_df.iterrows():
        body_lines.append(f'    x{i} = x[(x >= {w_range["w_min"]}) & (x <= {w_range["w_max"]})]')
        body_lines.append(f'    y{i} = cont_slope(x{i}, cont_{i}, slope_{i})')
        sub_df = df_list.loc[(df_list['w_min'] == w_range["w_min"]) & (df_list['w_max'] == w_range["w_max"])]
        for j, row in sub_df.iterrows():
            c = row['c_comp']
            body_lines.append(f'    y{i} = add_voigt(x{i}, y{i}, lambda0_{j}, b_{c}, n_{c}, f_{j}, gamma_{j}, v_rad_{c})')
        
        body_lines.append(f'    segments.append(y{i})')
    body_lines.append("    return np.concatenate(segments)")

    
    print("\n".join(body_lines))
    namespace = {"np": np, "cont_slope": cont_slope, "add_voigt": add_voigt}
    exec("\n".join(body_lines), namespace)
    func = namespace[func_name]
    func.__doc__ = f"Auto-generated n-component Voigt profile.\nParameters: {signature_str}"
    return func

generated_fit_function = make_multi_comp_voigt(fit_df)
vmodel = Model(generated_fit_function)
params = vmodel.make_params()


for i, w_range in enumerate(range_list):
    params[f'slope_{i}'].set(value=0)

# Fixed atomic parameters — generalized over all components
for i, row in fit_df.iterrows():
    params[f'lambda0_{i}'].set(value=row['WavelengthAir'], vary=False)
    params[f'f_{i}'].set(value=row['OscillatorStrength'], vary=False)
    params[f'gamma_{i}'].set(value=row['Gamma'], vary=False)

c_comps = fit_df['c_comp'].drop_duplicates().reset_index(drop=True)

for c in c_comps:
    # Shared parameters
    params[f'b_{c}'].set(    value=0.1,  min=0)
    params[f'n_{c}'].set(    value=1e11, min=0)
    params[f'v_rad_{c}'].set(value=0,    min=-20, max=20)


# old single cloud sight lines
scl_old = ['HD 23180', 'HD 24398', 'HD 144470', 'HD 147165', 'HD 147683', 'HD 149757', 'HD 166937', 'HD 170740',
           'HD 184915', 'HD 185418', 'HD 185859', 'HD 203532']

for star_name in scl_old[4:]:
    file_lists = []
    for wave_range in range_list:
        pythia = EdiblesOracle()
        file_list = pythia.getFilteredObsList(object=[star_name], MergedOnly=True, Wave=np.mean(wave_range))
        file_lists.append(file_list)

    file_lists = np.array(file_lists).T

    for obs_files in file_lists:  # one row = one observation epoch
        # Load and crop each component spectrum
        spectra = []
        for i, (fpath, wave_range) in enumerate(zip(obs_files, range_list)):
            spec = dr5_io.read_combined_spec(DATADIR / fpath, bary_corr=True)
            spec = util_functions.crop_spectrum(spec, *wave_range)
            my_order = np.nanmedian(spec[4])
            spec = spec[:, spec[4]==my_order]
            print(spec.shape)
            if len(spec) > 5:
                spec[1] = spec[6]
            spectra.append(spec[:5])

            params[f'cont_{i}'].value = np.nanmedian(spec[1])

        fit_spec = np.concatenate(spectra, axis=1)

        # plt.plot(fit_spec[0], fit_spec[1])
        # plt.show()

        result = vmodel.fit(fit_spec[1], params, x=fit_spec[0])

        plt.plot(fit_spec[1])
        plt.plot(result.best_fit)
        plt.show()
