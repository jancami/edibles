"""All the fitting machinery — builds the lmfit Parameters and runs the fit.

Called from main_run.astrovoigtfit_run. The flow is:

    astro_simultaneous_fit          <- the one you call from outside
        _generate_smart_knots       puts knots in sensible places
        lmfit.Model(continuum_voigt_wrapper).fit
            continuum_voigt_wrapper = spline(knots) * Voigt_fit_wrapper(...)
            Voigt_fit_wrapper       repacks lmfit's flat params and hands them
                                    to main.master_function (the actual model)

Only astro_simultaneous_fit is meant to be called from outside this file. The
*_wrapper functions are just model callables that lmfit drives, and
_generate_smart_knots is a helper.
"""

import numpy as np
from scipy.interpolate import CubicSpline, make_interp_spline
from lmfit import Parameters, Model

from edibles.projects.GUI.main import master_function


def Voigt_fit_wrapper(**params_list):
    """Voigt model that lmfit drives.

    lmfit passes everything as a flat dict of scalar params, so we have to
    unpack the per-species / per-transition / per-component pieces ourselves
    and stitch them back into arrays. Param naming convention:
        n_species, v_resolution, n_step, wavegrid              (globals)
        n_trans_i, n_component_i                               (per species i)
        lambda_i_j, f_i_j, gamma_i_j                           (transition j of i)
        v_rad_i_k, b_i_k, N_i_k                                (component k of i)
    Once unpacked we just call master_function and let it do the actual work.
    """
    n_species = params_list['n_species']
    wavegrid = params_list['wavegrid']
    v_resolution = params_list['v_resolution']
    n_step = params_list['n_step']

    species_data = {}
    for species_idx in range(int(n_species)):
        n_trans = int(float(params_list[f'n_trans_{species_idx}']))
        n_component = int(params_list[f'n_component_{species_idx}'])

        all_lambda = np.empty(n_trans)
        all_f = np.empty(n_trans)
        all_gamma = np.empty(n_trans)
        for i in range(n_trans):
            all_lambda[i] = params_list[f'lambda_{species_idx}_{i}']
            all_f[i] = params_list[f'f_{species_idx}_{i}']
            all_gamma[i] = params_list[f'gamma_{species_idx}_{i}']

        all_v_rad = np.empty(n_component)
        all_b = np.empty(n_component)
        all_N = np.empty(n_component)
        for i in range(n_component):
            all_v_rad[i] = params_list[f'v_rad_{species_idx}_{i}']
            all_b[i] = params_list[f'b_{species_idx}_{i}']
            all_N[i] = params_list[f'N_{species_idx}_{i}']

        species_data[species_idx] = {
            'lambda': all_lambda, 'f': all_f, 'gamma': all_gamma,
            'v_rad': all_v_rad, 'b': all_b, 'N': all_N,
        }

    # Translate to the suffix-based kwargs that master_function expects
    master_kwargs = {}
    for species_idx, data in species_data.items():
        if species_idx == 0:
            suffix = '1st'
        elif species_idx == 1:
            suffix = '2nd'
        elif species_idx == 2:
            suffix = '3rd'
        else:
            suffix = f'{species_idx + 1}th'

        master_kwargs[f'lambda_{suffix}'] = data['lambda']
        master_kwargs[f'f_{suffix}'] = data['f']
        master_kwargs[f'gamma_{suffix}'] = data['gamma']
        master_kwargs[f'b_{suffix}'] = data['b']
        master_kwargs[f'N_{suffix}'] = data['N']
        master_kwargs[f'v_rad_{suffix}'] = data['v_rad']

    return master_function(wavegrid, v_resolution=v_resolution,
                           n_step=n_step, **master_kwargs)


def continuum_voigt_wrapper(**params_list):
    """The actual model lmfit fits: continuum * absorption.

    Continuum is a spline through (knot_x, knot_y) — knot_x is fixed, knot_y
    is what gets optimized. Absorption is whatever Voigt_fit_wrapper returns.
    Multiplying them gives the observed (un-normalized) flux.
    """
    absorption_model = Voigt_fit_wrapper(**params_list)

    wavegrid = params_list['wavegrid']
    knot_x_array = params_list['knot_x_array']
    n_knots = len(knot_x_array)

    spline_order = params_list.get('spline_order', 3)
    if hasattr(spline_order, 'shape') and spline_order.shape == ():
        spline_order = int(spline_order)

    knot_y_values = [params_list[f'knot_y_{i}'] for i in range(n_knots)]

    if spline_order == 3:
        spline = CubicSpline(knot_x_array, knot_y_values)
    else:
        spline = make_interp_spline(knot_x_array, knot_y_values, k=spline_order)

    return spline(wavegrid) * absorption_model


def _generate_smart_knots(wavegrid, species_params, n_knots,
                          avoidance_width=0.5, absorption_ranges=None):
    """Pick spline knot positions that don't sit on top of absorption lines.

    Two paths:
      - If you pass absorption_ranges, knots inside those ranges get dropped
        and we add anchor knots just outside each range.
      - Otherwise we figure out where the lines actually land (rest wavelength
        Doppler-shifted by each component's v_rad guess) and nudge any knot
        that ends up too close.
    """
    w_min, w_max = np.min(wavegrid), np.max(wavegrid)
    initial_knots = np.linspace(w_min, w_max, n_knots)

    if absorption_ranges is not None:
        if isinstance(absorption_ranges, tuple):
            absorption_ranges = [absorption_ranges]

        final_knots = []
        for k in initial_knots:
            is_bad = any(start <= k <= end for (start, end) in absorption_ranges)
            if not is_bad:
                final_knots.append(k)

        for (start, end) in absorption_ranges:
            if start - avoidance_width >= w_min:
                final_knots.append(start - avoidance_width)
            if end + avoidance_width <= w_max:
                final_knots.append(end + avoidance_width)

        return np.sort(np.unique(final_knots))

    c_light = 299792.458  # km/s
    shifted_line_centers = []
    for species_idx in species_params:
        sp = species_params[species_idx]
        lambdas = sp.get('lambda', [])
        lambdas = [lambdas] if np.isscalar(lambdas) else np.asarray(lambdas).flatten()
        v_rads = sp.get('v_rad', [])
        v_rads = [v_rads] if np.isscalar(v_rads) else np.asarray(v_rads).flatten()

        for lam in lambdas:
            for v in v_rads:
                shifted_lam = lam * (1 + v / c_light)
                if w_min - 2.0 < shifted_lam < w_max + 2.0:
                    shifted_line_centers.append(shifted_lam)

    shifted_line_centers = np.array(shifted_line_centers)

    final_knots = []
    for k in initial_knots:
        if len(shifted_line_centers) > 0:
            dists = np.abs(shifted_line_centers - k)
            if np.any(dists < avoidance_width):
                new_k = k + avoidance_width
                if new_k > w_max:
                    new_k = k - avoidance_width
                k = new_k
        final_knots.append(k)

    return np.sort(np.unique(final_knots))

def astro_simultaneous_fit(wavegrid, ydata, species_params,
                           n_knots=10, knots_x_array=None,
                           v_resolution=0.0, n_step=25, std_dev=0.002,
                           avoidance_width=0, absorption_ranges=None,
                           spline_order=3):
    """Fit the spline continuum and the Voigt absorption together in one shot.

    Doing them in one fit (rather than normalizing first and then fitting the
    Voigt profile) avoids the error-propagation problem you get from the
    sequential approach.

    wavegrid, ydata        : the observed wavelength grid and the raw (un-normalized) flux
    species_params         : same per-species dict used everywhere else
                             (lambda / f / gamma / b / N / v_rad arrays; species
                             after the first can also set tie_b / tie_v_rad to
                             share those parameters with species 0)
    n_knots / knots_x_array: pick a count (we'll auto-place them) OR pass exact
                             knot positions
    spline_order           : 1 = linear, 2 = quadratic, 3 = cubic
    """
    if knots_x_array is None:
        knots_x_array = _generate_smart_knots(
            wavegrid, species_params, n_knots, avoidance_width, absorption_ranges
        )
    else:
        knots_x_array = np.asarray(knots_x_array)
    n_knots = len(knots_x_array)

    n_species = len(species_params)
    for species_idx in species_params:
        for key in ['lambda', 'f', 'gamma', 'b', 'N', 'v_rad']:
            species_params[species_idx][key] = np.asarray(species_params[species_idx][key])

    params = Parameters()
    params.add('n_species', value=n_species, vary=False)
    params.add('v_resolution', value=v_resolution, vary=False)
    params.add('n_step', value=n_step, vary=False)

    for species_idx in range(n_species):
        species_data = species_params[species_idx]

        n_trans = species_data['lambda'].size
        params.add(f'n_trans_{species_idx}', value=n_trans, vary=False)
        for i in range(n_trans):
            params.add(f'lambda_{species_idx}_{i}', value=species_data['lambda'][i], vary=False)
            params.add(f'f_{species_idx}_{i}', value=species_data['f'][i], vary=False)
            params.add(f'gamma_{species_idx}_{i}', value=species_data['gamma'][i], vary=False)

        n_component = species_data['v_rad'].size
        params.add(f'n_component_{species_idx}', value=n_component, vary=False)

        tie_b = species_data.get('tie_b', False) if species_idx > 0 else False
        tie_v_rad = species_data.get('tie_v_rad', False) if species_idx > 0 else False

        for i in range(n_component):
            if tie_b:
                params.add(f'b_{species_idx}_{i}', expr=f'b_0_{i}')
            else:
                params.add(f'b_{species_idx}_{i}', value=species_data['b'][i],
                           min=0.2, max=5.5, vary=True)

            params.add(f'N_{species_idx}_{i}', value=species_data['N'][i], min=0, vary=True)

            if tie_v_rad:
                params.add(f'v_rad_{species_idx}_{i}', expr=f'v_rad_0_{i}')
            else:
                params.add(f'v_rad_{species_idx}_{i}', value=species_data['v_rad'][i],
                        min=-30, max=30, vary=True)

    if len(ydata) == len(wavegrid):
        initial_knot_y = np.interp(knots_x_array, wavegrid, ydata)
    else:
        initial_knot_y = np.ones(n_knots)

    for i in range(n_knots):
        params.add(f'knot_y_{i}', value=initial_knot_y[i], min=0, vary=True)

    model = Model(continuum_voigt_wrapper,
                  independent_vars=['wavegrid', 'knot_x_array', 'spline_order'])
    return model.fit(ydata, params, wavegrid=wavegrid,
                     knot_x_array=knots_x_array,
                     spline_order=spline_order,
                     weights=1 / std_dev)
