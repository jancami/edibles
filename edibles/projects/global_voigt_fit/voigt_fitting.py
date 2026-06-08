from edibles.utils.edibles_oracle import EdiblesOracle
from edibles.projects.edr5_integration import dr5_io, util_functions
from edibles import DATADIR
import matplotlib.pyplot as plt
from edibles.utils.voigt_profile import voigt_optical_depth
from lmfit import Model
from lmfit.model import ModelResult
import numpy as np
from importlib.resources import files
import pandas as pd
from typing import Callable
from pprint import pprint
from PyAstronomy import pyasl


def add_voigt(x: np.array, y: np.array, lambda0: float, b: float, n: float, f: float, gamma: float, v_rad: float) -> np.array:
    """
    Adds a voigt absorption line of an element or molecule to an optical depth array.

    Parameters
    ----------
    x : np.array
        Wavelength array.
    y : np.array
        Optical depth array.
    lambda0 : float
        Central wavelength of voigt line.
    b : float
        Gaussian broadening parameter of voigt line.
    n : float
        Column desity of the species.
    f : float
        Oscillator strength of the transition.
    gamma : float
        Natural line width of the transition (Lorentzian component).
    v_rad : float
        Radial velocity.

    Returns
    -------
    np.array
        New optical depth array with additional voigt component.
    """
    return y + voigt_optical_depth(x, lambda0=lambda0, b=b, N=n, f=f, gamma=gamma, v_rad=v_rad)

def make_multi_comp_voigt(input_df: pd.DataFrame) -> Callable:
    """
    Makes a multi component voigt function based on an input DataFrame.

    Parameters
    ----------
    input_df : pd.DataFrame
        DataFrame including species names, wavelengths, values and cloud component numbers for v_rad, b and n. 

    Returns
    -------
    Callable
        Multi voigt absorption function.
    """
    # Get individual wavelength ranges
    range_df = input_df[['w_min', 'w_max']].drop_duplicates().reset_index(drop=True).sort_values(by=['w_min'])
    # Get individual velocity, Gaussian width and column density components.
    c_comps = input_df[['v_comp', 'b_comp', 'n_comp', 'Species']].drop_duplicates().reset_index(drop=True)

    # initialize list of parameters for the generated function
    comp_param_names = []
    # initialize list of v_rad, b and n components which are already included so they dont get included twice.
    incl_vrad = []
    incl_b = []
    incl_n = []
    for _, row in c_comps.iterrows():
        v_comp = row['v_comp']  # Radial velocity component
        b_comp = row['b_comp']  # Gaussian width component
        n_comp = row['n_comp']  # Column density component
        if v_comp not in incl_vrad:
            comp_param_names += [f'v_rad_{v_comp}']  # Add column density parameter of component v_rad 
        if b_comp not in incl_b:
            comp_param_names += [f'b_{b_comp}']  # Add b parameter of component b
        if n_comp not in incl_n:
            comp_param_names += [f'n_{n_comp}']  # Add column density parameter of component n 

            # add patameters to list of included components, so they dont get included twice.
            incl_vrad.append(v_comp)
            incl_b.append(b_comp)
            incl_n.append(n_comp)

    # Add atomic/molecular data for transitions
    for j, row in input_df.iterrows():
        comp_param_names += [f'lambda0_{j}', f'f_{j}', f'gamma_{j}']

    comp_param_names += ['v_rad_corr']  # Add radial velocity correction parameter

    all_params = ['x'] + comp_param_names  # Add wavelength as parameter
    signature_str = ", ".join(all_params)  # Join list to string
    func_name   = f"voigt_n_comp"  # Define function name
    
    # Start writing lines of function definition
    body_lines  = [f"def {func_name}({signature_str}):"]
    body_lines.append("    segments = []")  # List of flux output

    # Iterate trough wavelength windows
    for i, w_range in range_df.iterrows():
        # define instrumental resolution 
        # TODO: refine
        mean_wl = np.mean(w_range)  # mean wavelength of wavelength window
        if mean_wl < 5000:
            inst_res = 80000
        else:
            inst_res = 110000
        # Take spectrum within window
        body_lines.append(f'    x{i} = x[(x >= {w_range["w_min"]}) & (x <= {w_range["w_max"]})]')
        # Initalize optical depth of 0
        body_lines.append(f'    y{i} = np.zeros(len(x{i}))')
        # Get lines which are in the wavelength window
        sub_df = input_df.loc[(input_df['w_min'] == w_range["w_min"]) & (input_df['w_max'] == w_range["w_max"])]
        # Add voigt component for each line in the window
        for j, row in sub_df.iterrows():
            v_comp = row['v_comp']  # Radial velocity component
            b_comp = row['b_comp']  # Gaussian width component
            n_comp = row['n_comp']  # Column density component
            # Add component
            if row['WavelengthAir'] < 3303:
                body_lines.append(f'    y{i} = add_voigt(x{i}, y{i}, lambda0_{j}, b_{b_comp}, n_{n_comp}, f_{j}, gamma_{j}, v_rad_{v_comp} + v_rad_corr)')
            else:
                body_lines.append(f'    y{i} = add_voigt(x{i}, y{i}, lambda0_{j}, b_{b_comp}, n_{n_comp}, f_{j}, gamma_{j}, v_rad_{v_comp})')
        # convert to transmission spectrum
        body_lines.append(f'    y{i} = np.exp(-y{i})')
        # Add instrumental broadening
        body_lines.append(f'    y{i} = pyasl.instrBroadGaussFast(x{i}, y{i}, {inst_res}, edgeHandling="firstlast")')
        
        # Append flux to output flux array
        body_lines.append(f'    segments.append(y{i})')
    body_lines.append("    return np.concatenate(segments)")  # Return flux

    
    print("\n".join(body_lines))  # Join lines to one string
    namespace = {"np": np, "add_voigt": add_voigt, "pyasl": pyasl}  # Define used functions
    # Turn string into function
    exec("\n".join(body_lines), namespace)
    func = namespace[func_name]
    func.__doc__ = f"Auto-generated n-component Voigt profile.\nParameters: {signature_str}"
    return func  # Return function


def voigt_fit_wrapper(fit_df: pd.DataFrame, fit_spec: np.array, fit=True) -> ModelResult:
    """
    Generates a voigt model from a fitting DataFrame and an input spectrum.

    Parameters
    ----------
    fit_df : pd.DataFrame
        DataFrame including species names, wavelengths, values and cloud component numbers for v_rad, b and n. 
    fit_spec : np.array
        Fitted spectrum. [wave, flux, error]

    Returns
    -------
    ModelResult
        Fitting result of multi component voigt model.
    """
    # generate fitting function
    generated_fit_function = make_multi_comp_voigt(fit_df)
    # Make lmfit model
    vmodel = Model(generated_fit_function)
    # generate parameters
    params = vmodel.make_params()

    params['v_rad_corr'].set(value=0, min=-3, max=3)  # Set radial velocity correction parameter

    # Fixed atomic parameters — generalized over all components. Fixing them like this does not significantly decrease the fitting performance.
    for i, row in fit_df.iterrows():
        params[f'lambda0_{i}'].set(value=row['WavelengthAir'], vary=False)
        params[f'f_{i}'].set(value=row['OscillatorStrength'], vary=False)
        params[f'gamma_{i}'].set(value=row['Gamma'], vary=False)

    # extract the doppler and b components
    c_comps = fit_df[['v_comp', 'b_comp', 'n_comp', 'Species', 'v_rad_init', 'b_init', 'b_min', 'b_max', 'v_rad_min', 'v_rad_max', 'n_init']].drop_duplicates().reset_index(drop=True)

    # set initial values and bounds of b values and v_rad
    for _, row in c_comps.iterrows():
        v_comp = row['v_comp']
        b_comp = row['b_comp']
        n_comp = row['n_comp']
        # Shared parameters
        params[f'b_{b_comp}'].set(value=row['b_init'],  min=row['b_min'], max=row['b_max'], vary=True)
        params[f'v_rad_{v_comp}'].set(value=row[f'v_rad_init'],    min=row['v_rad_min'], max=row['v_rad_max'])
        params[f'n_{n_comp}'].set(value=row['n_init'], min=0)

    # convert weights for individual windows to array of length fit_spec
    weights = np.ones(len(fit_spec[0]))
    range_df = fit_df[['w_min', 'w_max']].drop_duplicates().reset_index(drop=True).sort_values(by=['w_min'])

    for i, wave_range in range_df.iterrows():
        weights[(fit_spec[0] >= wave_range['w_min']) & (fit_spec[0] <= wave_range['w_max'])] = fit_df.loc[fit_df['w_min'] == wave_range['w_min'], 'weight'].values[0]

    if fit:
        result = vmodel.fit(fit_spec[1], params, x=fit_spec[0], weights=weights/fit_spec[2]**2)  # Fit with combined weights of S/N and user-defined weights

        return result
    else:
        guess = vmodel.eval(params, x=fit_spec[0])
        return params, guess



def main():
    atomic_line_file = files('edibles') / 'data/auxiliary_data/line_catalogs/edibles_linelist_atoms.csv'
    atomic_line_list = pd.read_csv(atomic_line_file)
    print(atomic_line_list)

    obs_log = pd.read_csv(files('edibles') / 'data/DR5_ObsLog.csv')

    # NaI + KI
    elem_inds = [12, 15, 20, 21, 22, 23]
    elem_df = atomic_line_list.loc[elem_inds]
    range_list = [[3301, 3304], [4043, 4045], [5888, 5900], [7697, 7701]]

    # # NaI
    # elem_inds = [20, 21, 22, 23]
    # elem_df = atomic_line_list.loc[elem_inds]
    # range_list = [[3301, 3304], [5888, 5900]]

    # # LiI
    # elem_inds = [16, 17, 18, 19]
    # elem_df = atomic_line_list.loc[elem_inds]
    # range_list = [[6706, 6709.5]]

    n_comp = 1


    # old single cloud sight lines
    scl_old = ['HD 23180', 'HD 24398', 'HD 144470', 'HD 147165', 'HD 147683', 'HD 149757', 'HD 166937', 'HD 170740',
            'HD 184915', 'HD 185418', 'HD 185859', 'HD 203532']

    scl_rv = {'HD 23180': 13.3, 'HD 24398': 13.8, 'HD 144470': -10.1, 'HD 147165': -6.4, 'HD 147683': -0.8, 'HD 149757': -13.9,
            'HD 166937': -6.4, 'HD 170740': -10.1, 'HD 184915': -12.0, 'HD 185418': -10.1, 'HD 185859': -8.2, 'HD 203532': 14.2}

    for star_name in scl_old[0:]:
        fit_df = pd.DataFrame()
        # Define cloud components and initial values, limits
        for v_comp in range(n_comp):
            i_df = elem_df.copy()
            i_df.loc[:, 'v_comp'] = v_comp
            i_df.loc[:, f'v_rad_init'] = scl_rv[star_name]
            i_df.loc[:, f'v_rad_min'] = -20
            i_df.loc[:, f'v_rad_max'] = 20
            i_df.loc[:, 'b_comp'] = v_comp
            i_df.loc[:, f'b_init'] = 10
            i_df.loc[:, f'b_min'] = 0
            i_df.loc[:, f'b_max'] = 20
            fit_df = pd.concat((fit_df, i_df), ignore_index=True)

        # Define wave range for each line
        for i, row in fit_df.iterrows():
            w_min, w_max = [(low, high) for low, high in range_list if low <= row['WavelengthAir'] <= high][0]
            fit_df.loc[i, 'w_min'] = w_min
            fit_df.loc[i, 'w_max'] = w_max

        print(fit_df)

        file_lists = []
        for wave_range in range_list:
            pythia = EdiblesOracle()
            file_list = pythia.getFilteredObsList(object=[star_name], MergedOnly=True, Wave=np.mean(wave_range))
            file_lists.append(file_list[:2])

        file_lists = np.array(file_lists).T

        for obs_files in file_lists:  # one row = one observation epoch
            # Load and crop each component spectrum
            spectra = []
            for i, (fpath, wave_range) in enumerate(zip(obs_files, range_list)):
                spec = dr5_io.read_combined_spec(DATADIR / fpath, bary_corr=True)
                spec = util_functions.crop_spectrum(spec, *wave_range)
                my_order = np.nanmedian(spec[4])
                spec = spec[:, spec[4]==my_order]
                if len(spec) > 5:
                    if not np.isnan(spec[6]).all():
                        spec[1] = spec[6]

                x, y = pyasl.equidistantInterpolation(spec[0], spec[1], '2x')
                spectra.append(np.array([x, y]))
            fit_spec = np.concatenate(spectra, axis=1)

            plt.figure(figsize=(20, 10))
            plt.plot(fit_spec[0], fit_spec[1])
            plt.show()

            result = voigt_fit_wrapper(fit_df, fit_spec)

            
            pprint(result.best_values)

            plt.figure(figsize=(20, 10))
            plt.plot(fit_spec[1])
            plt.plot(result.best_fit)
            plt.show()

if __name__ == '__main__':
    main()