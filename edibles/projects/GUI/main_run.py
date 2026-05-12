import re
import numpy as np
from scipy.interpolate import CubicSpline, make_interp_spline

from astrovoigtfit import astro_simultaneous_fit


def get_species_data(species_file='species.txt'):
    """Parses species.txt to return available species and their wavelength ranges.

    Returns
    -------
    dict
        ``{species_name: {'wrange': [min, max], 'line': line_value}}``
    """
    species_data = {}
    with open(species_file) as f:
        lines = f.readlines()

    for line in lines[1:]:
        parts = line.split()
        if len(parts) < 6:
            continue

        name = parts[0]
        wrange_str = parts[-1].strip('[]')
        try:
            wrange = [float(x) for x in wrange_str.split(',')]
            line_val = parts[1]
            species_data[name] = {'wrange': wrange, 'line': line_val}
        except ValueError:
            pass

    return species_data


def _parse_bracketed_list(s):
    inner = s.strip().strip('[]')
    if not inner:
        return []
    return [float(x.strip()) for x in inner.split(',') if x.strip()]


def get_species_params(species_file, species_params, molecule):
    """Populate species_params with atomic data (lambda, f, gamma) from species_file."""
    with open(species_file) as f:
        lines = f.readlines()

    for idx, species in enumerate(molecule):
        for line in lines[1:]:
            if not line.strip():
                continue
            parts = line.split()
            if parts[0] != species:
                continue

            bracketed = re.findall(r'\[.*?\]', line)
            if len(bracketed) < 3:
                raise ValueError(f"Line for species {species} does not have expected bracketed fields: {line!r}")

            lambda_vals = _parse_bracketed_list(bracketed[0])
            f_vals = _parse_bracketed_list(bracketed[1])
            gamma_vals = _parse_bracketed_list(bracketed[2])

            # Handle case where lambda is [0] (e.g. K_4044) — fall back to the 'line' column
            if len(lambda_vals) == 1 and lambda_vals[0] == 0:
                try:
                    line_val = float(parts[1])
                    lambda_vals = [line_val]
                    print(f"Warning: Lambda is 0 for {species}, using line value {line_val}")
                except ValueError:
                    pass

            print(f"Species: {species}, Lambda: {lambda_vals}, f: {f_vals}, Gamma: {gamma_vals}")

            reordered = {
                'lambda': lambda_vals,
                'f': f_vals,
                'gamma': gamma_vals,
            }
            for key, value in species_params[idx].items():
                reordered[key] = value

            species_params[idx] = reordered
            break

    return species_params


def astrovoigtfit_run(wave, flux, molecules, species_params,
                      absorption_range, n_knots=10,
                      knots_x_array=None, species_file='species.txt',
                      spline_order=3, std_dev=0.002):
    """Run the joint spline-continuum + Voigt fit on (wave, flux).

    ``std_dev`` is the noise level of the co-added spectrum. It just gets
    forwarded to ``astro_simultaneous_fit``, which uses ``1/std_dev`` as the
    lmfit weight. The GUI's "Fit σ" box is what feeds this — bump it up if
    your co-added spectrum is noisier than the 0.002 default.

    Returns (fitresult, normalized_flux, continuum, residual_std). The last one
    is just ``std(flux - best_fit)`` for diagnostics — don't confuse it with
    the ``std_dev`` you passed in.
    """
    species_param_updated = get_species_params(species_file, species_params, molecules)

    fitresult = astro_simultaneous_fit(
        wavegrid=wave,
        ydata=flux,
        species_params=species_param_updated,
        n_knots=n_knots,
        knots_x_array=knots_x_array,
        v_resolution=3,
        n_step=25,
        std_dev=std_dev,
        absorption_ranges=[absorption_range] if absorption_range else None,
        spline_order=spline_order,
    )

    fitresult.params.pretty_print()
    print("chi-square value ", fitresult.chisqr)
    print("reduced chi-square value ", fitresult.redchi)
    print("FITTING RESULT :", fitresult.success)

    # Pull the continuum back out of the fit result. lmfit only gives us
    # best_fit = continuum * absorption, so we rebuild the spline from the
    # fitted knot_y values to get the continuum on its own. (No re-fitting
    # is happening here — we're just evaluating the curve the fit chose.)
    used_knots_x = fitresult.userkws.get('knot_x_array')
    knot_y_values = []
    i = 0
    while True:
        name = f'knot_y_{i}'
        if name in fitresult.params:
            knot_y_values.append(fitresult.params[name].value)
            i += 1
        else:
            break

    if used_knots_x is not None and len(knot_y_values) == len(used_knots_x):
        if spline_order == 3:
            spline = CubicSpline(used_knots_x, knot_y_values)
        else:
            spline = make_interp_spline(used_knots_x, knot_y_values, k=spline_order)
        continuum = spline(wave)
    else:
        # Shouldn't get here under normal operation — means the fit didn't
        # produce knot_y_i params or the count is off. Just punt with a flat
        # continuum so the GUI can still draw something.
        print("Warning: Could not reconstruct continuum perfectly.")
        continuum = np.ones_like(wave)

    normalized_flux = flux / continuum
    std_dev = np.std(flux - fitresult.best_fit)

    return fitresult, normalized_flux, continuum, std_dev
