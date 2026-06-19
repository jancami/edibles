import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial import Chebyshev
import pandas as pd
from scipy.optimize import minimize as minimize
from lmfit import Parameters, minimize
from edibles.projects.DIBCode_Empirical.calculate_chisqr_manual import calculate_chisqr_manual, 
from edibles.projects.DIBCode_Empirical.calculate_reduced_chi_square import calculate_reduced_chi_square
from edibles.projects.DIBCode_Empirical.format_params_grouped import format_params_grouped
from edibles.projects.DIBCode_Empirical.define_aic_bic import calculate_aic, calculate_bic
from edibles.projects.DIBCode_Empirical.initialize_gaussian import initialize_gaussian
from edibles.projects.DIBCode_Empirical.initialize_lorentzian import initialize_lorentzian
from edibles.projects.DIBCode_Empirical.initialize_chebyshev import initialize_chebyshev
from edibles.projects.DIBCode_Empirical.plot_iteration import plot_iteration
from edibles.projects.DIBCode_Empirical.composite_model import composite_model
def iterative_fit(x, data, uncertainties=None,
                  max_iterations=20,
                  significance_level=0.05,
                  criterion='bic',
                  verbose=True,
                  plot=True,
                  true_params=None):
    """Iterative model building with forward selection for curve fitting.

    Initalizes the routine with a constant baseline where it then
    progresses to the main iterative loop by adding one unique component everytime.
    The routine then evaluates and picks the best candidate and decides wheather
    to accept of reject based on fit.

    Args:
        x: Independent variable data points (wavelengths) as a NumPy array.
        data: Observed dependent variable data points.
        uncertainties: Array-like or None measurement error at each data point
        max_iterations: integer, prevents infinite fitting
        significance_level: float, threshold for statistical significance
        criterion: str ('bic', 'aic', or 'ftest'), how the algorithm decides whether a new component is worth adding
        verbose: bool, print detailed progress
        plot: bool, visualize each iteration
        true_params: dictionary or None, helps evaluate accuracy of recovered parameters


    Returns:
        Best_result: diagnostic information about the final fit, chi-squared, number of parameters, residuals, and success flag
        best_params: the best parameters from the final iteration determined by the routine
        history: the record of the step-by-step model evaluation

    """

    n_data = len(x)
    history = []

    # Verify uncertainties
    if uncertainties is not None:
        print(f"Using uncertainties: min={np.min(uncertainties):.3f}, "
              f"max={np.max(uncertainties):.3f}, mean={np.mean(uncertainties):.3f}")
    else:
        print("WARNING: No uncertainties provided. Chi-square will be unweighted.")

    # ========================================================================
    # STEP 1: Initialize with JUST baseline (no components)
    # ========================================================================

    params = Parameters()
    params.add('c0', value=np.median(data))  # Just constant baseline

    # Fit initial model
    result = minimize(residual, params, args=(x, data, uncertainties), method='leastsq')
    # **CRITICAL FIX: Recalculate chi-square with fitted parameters**
    result.chisqr = calculate_chisqr_manual(result.params, x, data, uncertainties)

    # Verify chi-square calculation
    chisqr_manual = calculate_chisqr_manual(result.params, x, data, uncertainties)
    if verbose:
        print(f"\n{'=' * 70}")
        print(f"INITIAL MODEL: Constant baseline only")
        print(f"  c0 = {result.params['c0'].value:.4f}")
        print(f"  Chi-square (lmfit): {result.chisqr:.2f}")
        print(f"  Chi-square (manual): {chisqr_manual:.2f}")
        print(f"  Reduced χ²: {calculate_reduced_chi_square(result, n_data):.4f}")
        print(f"  DOF: {n_data - result.nvarys}")
        print(f"  AIC: {calculate_aic(result.chisqr, result.nvarys, n_data):.2f}")
        print(f"  BIC: {calculate_bic(result.chisqr, result.nvarys, n_data):.2f}")

    best_result = result
    best_params = result.params.copy()

    history.append({
        'iteration': 0,
        'component': 'Baseline only',
        'n_params': result.nvarys,
        'chisqr': result.chisqr,
        'reduced_chisqr': calculate_reduced_chi_square(result, n_data),
        'aic': calculate_aic(result.chisqr, result.nvarys, n_data),
        'bic': calculate_bic(result.chisqr, result.nvarys, n_data)
    })

    # ========================================================================
    # STEP 2: Iterative improvement
    # ========================================================================

    for iteration in range(1, max_iterations + 1):

        if verbose:
            print(f"\n{'=' * 70}")
            print(f"ITERATION {iteration}: Testing candidate additions...")
            print(f"\nCurrent best model parameters:")
            print(format_params_grouped(best_params))
            print(f"Current best chi²: {best_result.chisqr:.2f}")

        current_residuals = data - composite_model(x, best_params)

        # Count current components
        n_gaussians = sum(1 for name in best_params if name.startswith('g') and name.endswith('_amp'))
        n_lorentzians = sum(1 for name in best_params if name.startswith('l') and name.endswith('_amp'))
        n_cheb = sum(1 for name in best_params if name.startswith('c'))

        # ====================================================================
        # Test all possible additions
        # ====================================================================

        candidates = []

        # --- Candidate 1: Add Gaussian ---
        print(f"\nTrying to add Gaussian {n_gaussians}:")
        test_params = best_params.copy()

        # **FIX: Explicitly set all parameters to vary=True**
        for param in test_params.values():
            param.vary = True

        new_gauss = initialize_gaussian(x, current_residuals, prefix=f'g{n_gaussians}_')
        for name, config in new_gauss.items():
            test_params.add(name, **config)

        test_result = minimize(residual, test_params, args=(x, data, uncertainties), method='leastsq')
        # **CRITICAL FIX: Recalculate chi-square with fitted parameters**
        test_result.chisqr = calculate_chisqr_manual(test_result.params, x, data, uncertainties)

        print(f"  Chi-square AFTER fit: {test_result.chisqr:.2f}")
        print(f"  Fit succeeded: {test_result.success}")
        print(f"  All parameters after refit:")
        print("    " + format_params_grouped(test_result.params).replace('\n', '\n    '))

        candidates.append({
            'name': f'+ Gaussian {n_gaussians}',
            'params': test_result.params,
            'result': test_result,
            'type': 'gaussian'
        })

        # --- Candidate 2: Add Lorentzian ---
        print(f"\nTrying to add Lorentzian {n_lorentzians}:")
        test_params = best_params.copy()

        # **FIX: Explicitly set all parameters to vary=True**
        for param in test_params.values():
            param.vary = True

        new_lorentz = initialize_lorentzian(x, current_residuals, prefix=f'l{n_lorentzians}_')
        for name, config in new_lorentz.items():
            test_params.add(name, **config)

        test_result = minimize(residual, test_params, args=(x, data, uncertainties), method='leastsq')
        # **CRITICAL FIX: Recalculate chi-square with fitted parameters**
        test_result.chisqr = calculate_chisqr_manual(test_result.params, x, data, uncertainties)

        print(f"  Chi-square AFTER fit: {test_result.chisqr:.2f}")
        print(f"  Fit succeeded: {test_result.success}")
        print(f"  All parameters after refit:")
        print("    " + format_params_grouped(test_result.params).replace('\n', '\n    '))

        candidates.append({
            'name': f'+ Lorentzian {n_lorentzians}',
            'params': test_result.params,
            'result': test_result,
            'type': 'lorentzian'
        })

        # --- Candidate 3: Add Chebyshev term ---
        if n_cheb < 6:  # Limit polynomial order
            print(f"\nTrying to add Chebyshev c{n_cheb}:")
            test_params = best_params.copy()

            # **FIX: Explicitly set all parameters to vary=True**
            for param in test_params.values():
                param.vary = True

            new_cheb = initialize_chebyshev(n_cheb, prefix='c')
            for name, config in new_cheb.items():
                test_params.add(name, **config)

            test_result = minimize(residual, test_params, args=(x, data, uncertainties), method='leastsq')
            # **CRITICAL FIX: Recalculate chi-square with fitted parameters**
            test_result.chisqr = calculate_chisqr_manual(test_result.params, x, data, uncertainties)

            print(f"  Chi-square AFTER fit: {test_result.chisqr:.2f}")
            print(f"  Fit succeeded: {test_result.success}")
            print(f"  All parameters after refit:")
            print("    " + format_params_grouped(test_result.params).replace('\n', '\n    '))

            candidates.append({
                'name': f'+ Chebyshev c{n_cheb}',
                'params': test_result.params,
                'result': test_result,
                'type': 'chebyshev'
            })

        # ====================================================================
        # Evaluate candidates
        # ====================================================================

        for cand in candidates:
            cand['aic'] = calculate_aic(cand['result'].chisqr, cand['result'].nvarys, n_data)
            cand['bic'] = calculate_bic(cand['result'].chisqr, cand['result'].nvarys, n_data)
            cand['reduced_chisqr'] = calculate_reduced_chi_square(cand['result'], n_data)

            # F-test vs current best
            F, p_value = f_test(
                best_result.chisqr,
                cand['result'].chisqr,
                n_data - best_result.nvarys,
                n_data - cand['result'].nvarys
            )
            cand['f_statistic'] = F
            cand['p_value'] = p_value

            if verbose:
                print(f"\n  {cand['name']}:")
                print(f"    Reduced χ²: {cand['reduced_chisqr']:.4f}")
                print(f"    BIC: {cand['bic']:.2f}")
                print(f"    AIC: {cand['aic']:.2f}")
                print(f"    p-value: {cand['p_value']:.4e}")

        # ====================================================================
        # Select best candidate
        # ====================================================================

        if criterion == 'aic':
            best_candidate_idx = min(range(len(candidates)), key=lambda i: candidates[i]['aic'])
            best_candidate = candidates[best_candidate_idx]
            improvement = best_result.chisqr - best_candidate['result'].chisqr
            accept = best_candidate['aic'] < calculate_aic(best_result.chisqr, best_result.nvarys, n_data)

        elif criterion == 'bic':
            best_candidate_idx = min(range(len(candidates)), key=lambda i: candidates[i]['bic'])
            best_candidate = candidates[best_candidate_idx]
            improvement = best_result.chisqr - best_candidate['result'].chisqr
            accept = best_candidate['bic'] < calculate_bic(best_result.chisqr, best_result.nvarys, n_data)

        else:  # 'ftest'
            best_candidate_idx = min(range(len(candidates)), key=lambda i: candidates[i]['result'].chisqr)
            best_candidate = candidates[best_candidate_idx]
            improvement = best_result.chisqr - best_candidate['result'].chisqr
            accept = best_candidate['p_value'] < significance_level

        # ====================================================================
        # Visualize
        # ====================================================================

        if plot:
            plot_iteration(x, data, best_params, best_result, candidates, best_candidate_idx,
                           iteration, n_data, uncertainties)

        # ====================================================================
        # Decision: accept or stop
        # ====================================================================

        if accept and improvement > 0:
            best_result = best_candidate['result']
            best_params = best_candidate['params']

            if verbose:
                print(f"\n  ✓ ACCEPTED: {best_candidate['name']}")
                print(f"    χ² improvement: {improvement:.2f}")

            history.append({
                'iteration': iteration,
                'component': best_candidate['name'],
                'n_params': best_result.nvarys,
                'chisqr': best_result.chisqr,
                'reduced_chisqr': best_candidate['reduced_chisqr'],
                'aic': best_candidate['aic'],
                'bic': best_candidate['bic'],
                'p_value': best_candidate.get('p_value', None)
            })

        else:
            if verbose:
                print(f"\n  ✗ STOPPING: No significant improvement")
                print(f"    Best candidate was {best_candidate['name']}")
                print(f"    but did not meet acceptance criterion")
            break

    # ========================================================================
    # STEP 3: Compare with true parameters if provided
    # ========================================================================

    if true_params is not None and verbose:
        print(f"\n{'=' * 70}")
        print("COMPARISON: TRUE vs FITTED PARAMETERS")
        print(f"{'=' * 70}")

        for param_type in ['gaussian', 'lorentzian', 'chebyshev']:
            if param_type == 'gaussian':
                prefix = 'g'
                suffix = ['amp', 'cen', 'wid']
                true_list = true_params.get('gaussians', [])
            elif param_type == 'lorentzian':
                prefix = 'l'
                suffix = ['amp', 'cen', 'wid']
                true_list = true_params.get('lorentzians', [])
            else:
                prefix = 'c'
                suffix = None
                true_list = true_params.get('chebyshev', [])

            if param_type in ['gaussian', 'lorentzian']:
                for i, true_comp in enumerate(true_list):
                    print(f"\n{param_type.upper()} {i}:")
                    for s in suffix:
                        true_val = true_comp.get(s, None)
                        param_name = f'{prefix}{i}_{s}'
                        if param_name in best_params:
                            fit_val = best_params[param_name].value
                            if true_val is not None:
                                diff = fit_val - true_val
                                pct = 100 * diff / true_val if true_val != 0 else 0
                                print(f"  {s:4s}: True={true_val:8.3f}, Fit={fit_val:8.3f}, "
                                      f"Diff={diff:+8.3f} ({pct:+6.1f}%)")
                        else:
                            if true_val is not None:
                                print(f"  {s:4s}: True={true_val:8.3f}, Fit=NOT FOUND")
            else:
                print(f"\nCHEBYSHEV/POLYNOMIAL COEFFICIENTS:")
                for i, true_val in enumerate(true_list):
                    param_name = f'c{i}'
                    if param_name in best_params:
                        fit_val = best_params[param_name].value
                        diff = fit_val - true_val
                        print(f"  c{i}: True={true_val:8.4f}, Fit={fit_val:8.4f}, Diff={diff:+8.4f}")
                    else:
                        print(f"  c{i}: True={true_val:8.4f}, Fit=NOT FOUND")

    if verbose:
        print(f"\n{'=' * 70}")
        print(f"FINAL MODEL:")
        n_g = sum(1 for name in best_params if name.startswith('g') and name.endswith('_amp'))
        n_l = sum(1 for name in best_params if name.startswith('l') and name.endswith('_amp'))
        n_c = sum(1 for name in best_params if name.startswith('c'))
        print(f"  {n_g} Gaussian(s) + {n_l} Lorentzian(s) + Chebyshev order {n_c - 1}")
        print(f"  Total parameters: {best_result.nvarys}")
        print(f"  Reduced χ²: {calculate_reduced_chi_square(best_result, n_data):.4f}")
        print(f"{'=' * 70}\n")

    return best_result, best_params, history
