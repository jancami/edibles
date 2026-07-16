import numpy as np
from numpy.polynomial import Chebyshev
from edibles.projects.DIBCode_Empirical.plot_final_fit import plot_final_fit
from edibles.projects.DIBCode_Empirical.get_component_list import get_component_list
from edibles.projects.DIBCode_Empirical.fitting_objective import fitting_objective
from edibles.projects.DIBCode_Empirical.plot_alignment_verification import plot_alignment_verification
from edibles.projects.DIBCode_Empirical.calculate_velocity_alignment import calculate_velocity_alignment

def perform_fit(x_obs, y_obs, params, figpath, c0, c1=0, target='Target', fitlimit=0, use_mask=False, model_v_shift=0):
    """
    Fits absorption profiles with a 1st order Chebyshev continuum.
    Args:
        x_obs: Observed wavelength values.
        y_obs: Coadded flux as a NumPy array.
        params: Parameters of the model.
        c0: First-order Chebyeshev polynomial parameter
        c1: Second-order Chebyeshev polynomial parameter
        target: Targeted sightline.
        fitlimit: Extent of the region the model is trying to fit.
        use_mask: Boolean - whether or not the mask is used for the fit
        model_v_shift: Velocity shift of the model.
    Returns:
        best_fit: Fit with the velocity that gives the highest correlation.
    """

    for key in params.keys():
        if key.endswith('_cen'):
            params[key] *= (1 + (model_v_shift / 299792.458))

    x_ref = np.asarray(x_obs).copy()

    x_min, x_max = x_ref.min(), x_ref.max()
    x_norm = 2 * (x_ref - x_min) / (x_max - x_min) - 1
    initial_continuum = c0 + c1 * x_norm

    comp_list, comp_ids = get_component_list(x_ref, params)
    ref_y_profiles = np.sum(comp_list, axis=0)
    model_template = initial_continuum + ref_y_profiles

    # 3. Masking and Alignment
    mask = (model_template < fitlimit) if use_mask else np.ones(len(model_template), dtype=bool)
    if use_mask and np.sum(mask) < 2:
        mask = np.ones(len(model_template), dtype=bool)

    v_grid = np.arange(-40.2, 40.2, 0.2)
    best_v, max_r, correlations, best_y_data = calculate_velocity_alignment(x_ref, y_obs, model_template, v_grid, mask)

    plot_alignment_verification(x_ref, y_obs, best_y_data, model_template, v_grid, correlations, best_v, max_r, mask,
                                use_mask, target, figpath)

    # 4. Optimization Loop
    n = len(best_y_data)
    best_bic, best_fit = np.inf, None

    for order in range(1, 25):
        # We pass initial_continuum to the objective so it knows the baseline
        from scipy.optimize import minimize as scipy_minimize
        res = scipy_minimize(fitting_objective, x0=[1.0],
                             args=(x_ref, best_y_data, initial_continuum, ref_y_profiles, order),
                             bounds=[(0, None)])

        scale_opt, rss = float(res.x[0]), res.fun
        k = (order + 1) + 1
        current_bic = n * np.log(rss / n) + k * np.log(n)

        if current_bic < best_bic:
            best_bic = current_bic
            final_scaled_absorption = scale_opt * ref_y_profiles

            # Fit the polynomial to the residuals (Data - Scaled Absorption)
            # The result 'poly' already includes the baseline offset, no need to add c0/c1 again manually
            poly = Chebyshev.fit(x_ref, best_y_data - final_scaled_absorption, deg=order)
            y_cont = poly(x_ref)

            best_fit = {
                'x_corr': x_ref, 'y_obs': best_y_data, 'y_cont': y_cont,
                'scaled_individual_comps': [scale_opt * p for p in comp_list],
                'comp_ids': comp_ids,
                'y_total_model': y_cont + final_scaled_absorption,
                'residuals': best_y_data - (y_cont + final_scaled_absorption),
                'order': order, 'scale_factor': scale_opt,
                'v_shift': best_v, 'max_correlation': max_r
            }
        else:
            break

    if best_fit:
        print(
            f"Fit Complete: Shift={best_fit['v_shift']:.2f} km/s, Scale={best_fit['scale_factor']:.5f}, Poly Order={best_fit['order']}")
        plot_final_fit(best_fit, target, figpath)

    return best_fit
