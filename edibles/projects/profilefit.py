import numpy as np
from lmfit import Parameters, minimize
from numpy.polynomial import chebyshev
from scipy import stats
import matplotlib.pyplot as plt

# ============================================================================
# PART 1: MODEL DEFINITION (stays constant throughout)
# ============================================================================

def composite_model(x, params):
    """
    Flexible model: sum of Gaussians, Lorentzians, and Chebyshev polynomial.
    
    Parameters naming convention:
    - Gaussian i: 'g{i}_amp', 'g{i}_cen', 'g{i}_wid'
    - Lorentzian i: 'l{i}_amp', 'l{i}_cen', 'l{i}_wid'
    - Chebyshev: 'c0', 'c1', 'c2', ...
    
    NOTE: For absorption features, amp should be NEGATIVE.
    """
    y = np.zeros_like(x)
    
    # Add Gaussian components
    i = 0
    while f'g{i}_amp' in params:
        amp = params[f'g{i}_amp'].value
        cen = params[f'g{i}_cen'].value
        wid = params[f'g{i}_wid'].value
        y += amp * np.exp(-((x - cen) / wid)**2)
        i += 1
    
    # Add Lorentzian components
    i = 0
    while f'l{i}_amp' in params:
        amp = params[f'l{i}_amp'].value
        cen = params[f'l{i}_cen'].value
        wid = params[f'l{i}_wid'].value
        y += amp * wid**2 / ((x - cen)**2 + wid**2)
        i += 1
    
    # Add Chebyshev polynomial
    cheb_coeffs = []
    i = 0
    while f'c{i}' in params:
        cheb_coeffs.append(params[f'c{i}'].value)
        i += 1
    
    if cheb_coeffs:
        # Normalize x to [-1, 1] for numerical stability
        x_norm = 2 * (x - x.min()) / (x.max() - x.min()) - 1
        y += chebyshev.chebval(x_norm, cheb_coeffs)
    
    return y


def residual(params, x, data, uncertainties=None):
    """
    Residual function for least-squares minimization.
    
    Returns weighted residuals: (data - model) / uncertainties
    This ensures chi^2 = sum(residuals^2) is properly weighted.
    """
    model = composite_model(x, params)
    if uncertainties is not None:
        return (data - model) / uncertainties
    return data - model


def calculate_chisqr_manual(params, x, data, uncertainties):
    """
    Manually calculate chi-square for verification.
    chi^2 = sum((data - model)^2 / sigma^2)
    """
    model = composite_model(x, params)
    if uncertainties is not None:
        chisqr = np.sum(((data - model) / uncertainties)**2)
    else:
        chisqr = np.sum((data - model)**2)
    return chisqr


# ============================================================================
# PART 2: IMPROVED PARAMETER INITIALIZATION
# ============================================================================

def estimate_feature_width(x, residuals, center_idx):
    """
    Estimate width of a feature by finding where it drops to half maximum.
    """
    peak_val = abs(residuals[center_idx])
    half_max = peak_val / 2.0
    
    # Search left
    left_idx = center_idx
    while left_idx > 0 and abs(residuals[left_idx]) > half_max:
        left_idx -= 1
    
    # Search right
    right_idx = center_idx
    while right_idx < len(residuals) - 1 and abs(residuals[right_idx]) > half_max:
        right_idx += 1
    
    # FWHM estimate
    fwhm = x[right_idx] - x[left_idx]
    
    # Convert to sigma (for our Gaussian parameterization)
    width_guess = fwhm / 1.665 if fwhm > 0 else np.ptp(x) * 0.01
    
    # Sanity check
    x_step = np.median(np.diff(x))
    min_width = 2 * x_step  # At least 2 data points wide
    max_width = np.ptp(x) * 0.2  # At most 20% of range
    
    width_guess = np.clip(width_guess, min_width, max_width)
    
    return width_guess


def initialize_gaussian(x, residuals, prefix='g0_'):
    """
    Smart initialization: place Gaussian at strongest absorption feature.
    
    For ABSORPTION: amplitude should be NEGATIVE.
    """
    # Find location of most negative residual (deepest absorption)
    idx_min = np.argmin(residuals)
    center_guess = x[idx_min]
    amp_guess = residuals[idx_min]  # Should be negative for absorption
    
    # Estimate width from feature shape
    width_guess = estimate_feature_width(x, residuals, idx_min)
    
    print(f"  Initializing {prefix[:-1]}: amp={amp_guess:.2f}, cen={center_guess:.1f}, wid={width_guess:.1f}")
    
    return {
        f'{prefix}amp': {'value': amp_guess, 'min': -np.inf, 'max': 0},  # ABSORPTION: amp <= 0
        f'{prefix}cen': {'value': center_guess, 'min': x.min(), 'max': x.max()},
        f'{prefix}wid': {'value': width_guess, 'min': np.median(np.diff(x)), 'max': np.ptp(x) * 0.3}
    }


def initialize_lorentzian(x, residuals, prefix='l0_'):
    """
    Smart initialization for Lorentzian (absorption).
    """
    idx_min = np.argmin(residuals)
    center_guess = x[idx_min]
    amp_guess = residuals[idx_min]  # Should be negative for absorption
    
    # Estimate width (Lorentzian HWHM)
    width_guess = estimate_feature_width(x, residuals, idx_min)
    
    print(f"  Initializing {prefix[:-1]}: amp={amp_guess:.2f}, cen={center_guess:.1f}, wid={width_guess:.1f}")
    
    return {
        f'{prefix}amp': {'value': amp_guess, 'min': -np.inf, 'max': 0},  # ABSORPTION: amp <= 0
        f'{prefix}cen': {'value': center_guess, 'min': x.min(), 'max': x.max()},
        f'{prefix}wid': {'value': width_guess, 'min': np.median(np.diff(x)), 'max': np.ptp(x) * 0.3}
    }


def initialize_chebyshev(order, prefix='c'):
    """Initialize next Chebyshev coefficient to zero."""
    return {f'{prefix}{order}': {'value': 0.0}}


# ============================================================================
# PART 3: PARAMETER FORMATTING FOR DISPLAY
# ============================================================================

def format_params_grouped(params):
    """
    Format parameters grouped by type for display.
    
    Returns a formatted string with:
    - Continuum (Chebyshev coefficients)
    - Gaussians
    - Lorentzians
    """
    lines = []
    
    # Continuum (Chebyshev)
    cheb_params = []
    i = 0
    while f'c{i}' in params:
        val = params[f'c{i}'].value
        cheb_params.append(f'c{i}={val:.4f}')
        i += 1
    if cheb_params:
        lines.append('Continuum: ' + ', '.join(cheb_params))
    
    # Gaussians
    gauss_list = []
    i = 0
    while f'g{i}_amp' in params:
        amp = params[f'g{i}_amp'].value
        cen = params[f'g{i}_cen'].value
        wid = params[f'g{i}_wid'].value
        gauss_list.append(f'g{i}[{amp:.2f}, {cen:.0f}, {wid:.1f}]')
        i += 1
    if gauss_list:
        lines.append('Gaussians: ' + ', '.join(gauss_list))
    
    # Lorentzians
    lorentz_list = []
    i = 0
    while f'l{i}_amp' in params:
        amp = params[f'l{i}_amp'].value
        cen = params[f'l{i}_cen'].value
        wid = params[f'l{i}_wid'].value
        lorentz_list.append(f'l{i}[{amp:.2f}, {cen:.0f}, {wid:.1f}]')
        i += 1
    if lorentz_list:
        lines.append('Lorentzians: ' + ', '.join(lorentz_list))
    
    return '\n'.join(lines)


def check_vary_status(params):
    """Check and report which parameters are set to vary."""
    vary_true = []
    vary_false = []
    for name, param in params.items():
        if param.vary:
            vary_true.append(name)
        else:
            vary_false.append(name)
    return vary_true, vary_false


# ============================================================================
# PART 4: MODEL COMPARISON
# ============================================================================

def calculate_aic(chi_square, n_params, n_data):
    """Akaike Information Criterion (lower is better)."""
    return chi_square + 2 * n_params


def calculate_bic(chi_square, n_params, n_data):
    """Bayesian Information Criterion (lower is better, penalizes complexity more)."""
    return chi_square + n_params * np.log(n_data)


def calculate_reduced_chi_square(result, n_data):
    """Reduced chi-square statistic."""
    return result.chisqr / (n_data - result.nvarys)


def f_test(chi2_simple, chi2_complex, df_simple, df_complex):
    """
    F-test for nested models.
    
    Returns F-statistic and p-value.
    H0: complex model does not significantly improve fit.
    """
    delta_chi2 = chi2_simple - chi2_complex
    delta_df = df_simple - df_complex
    
    if delta_df <= 0 or delta_chi2 <= 0:
        return 0, 1.0  # No improvement
    
    F = (delta_chi2 / delta_df) / (chi2_complex / df_complex)
    p_value = 1 - stats.f.cdf(F, delta_df, df_complex)
    
    return F, p_value


def get_component_params_text(params, cand_type, cand_idx):
    """Extract and format parameters for display."""
    if cand_type == 'gaussian':
        prefix = f'g{cand_idx}_'
        if f'{prefix}amp' in params:
            amp = params[f'{prefix}amp'].value
            cen = params[f'{prefix}cen'].value
            wid = params[f'{prefix}wid'].value
            return f"Amp={amp:.2f}, Cen={cen:.1f}, Wid={wid:.1f}"
    elif cand_type == 'lorentzian':
        prefix = f'l{cand_idx}_'
        if f'{prefix}amp' in params:
            amp = params[f'{prefix}amp'].value
            cen = params[f'{prefix}cen'].value
            wid = params[f'{prefix}wid'].value
            return f"Amp={amp:.2f}, Cen={cen:.1f}, Wid={wid:.1f}"
    elif cand_type == 'chebyshev':
        if f'c{cand_idx}' in params:
            val = params[f'c{cand_idx}'].value
            return f"c{cand_idx}={val:.4f}"
    return ""


# ============================================================================
# PART 5: VISUALIZATION (RESIZED FOR MACBOOK)
# ============================================================================

def plot_iteration(x, data, best_params, best_result, candidates, best_candidate_idx, 
                   iteration, n_data, uncertainties=None):
    """
    Plot each candidate in a separate subplot with fit statistics.
    Sized appropriately for MacBook screens.
    NOW WITH GROUPED PARAMETER DISPLAY AND CHI-SQUARE FOR CURRENT BEST!
    """
    n_candidates = len(candidates)
    
    # Smaller figure size for MacBook (max height ~10 inches)
    fig_height = min(2.5 * (n_candidates + 1), 10)
    fig, axes = plt.subplots(n_candidates + 1, 2, figsize=(12, fig_height))
    fig.suptitle(f'Iteration {iteration}', fontsize=14, fontweight='bold', y=0.998)
    
    # If only one row of subplots, make it 2D
    if n_candidates == 0:
        axes = axes.reshape(1, 2)
    
    # ========================================================================
    # Top row: Current best model
    # ========================================================================
    ax_best = axes[0, 0]
    ax_resid_best = axes[0, 1]
    
    if uncertainties is not None:
        ax_best.errorbar(x, data, yerr=uncertainties, fmt='k.', alpha=0.3, 
                        label='Data', capsize=0, elinewidth=0.5, markersize=2)
    else:
        ax_best.plot(x, data, 'k.', alpha=0.3, markersize=2, label='Data')
    
    best_model = composite_model(x, best_params)
    ax_best.plot(x, best_model, 'b-', lw=1.5, label='Current best')
    
    # Count components
    n_g = sum(1 for name in best_params if name.startswith('g') and name.endswith('_amp'))
    n_l = sum(1 for name in best_params if name.startswith('l') and name.endswith('_amp'))
    n_c = sum(1 for name in best_params if name.startswith('c'))
    
    ax_best.set_ylabel('Flux', fontsize=9)
    ax_best.set_title(f'CURRENT BEST: {n_g}G + {n_l}L + Cheb({n_c-1})', 
                     fontsize=10, fontweight='bold', color='blue')
    ax_best.legend(loc='lower right', fontsize=8)  # **CHANGED: moved to lower right**
    ax_best.grid(True, alpha=0.3)
    ax_best.tick_params(labelsize=8)
    
    # Add fit statistics for current best (LEFT SIDE)
    best_reduced_chisqr = calculate_reduced_chi_square(best_result, n_data)
    best_aic = calculate_aic(best_result.chisqr, best_result.nvarys, n_data)
    best_bic = calculate_bic(best_result.chisqr, best_result.nvarys, n_data)
    
    stats_text = (f"χ²ᵣ={best_reduced_chisqr:.3f}\n"
                 f"BIC={best_bic:.1f}\n"
                 f"AIC={best_aic:.1f}")
    ax_best.text(0.02, 0.98, stats_text, transform=ax_best.transAxes,
                fontsize=7, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.7))
    
    # Add parameter text box for current best (RIGHT SIDE)
    param_text = format_params_grouped(best_params)
    ax_best.text(0.98, 0.97, param_text, transform=ax_best.transAxes,
                fontsize=6, verticalalignment='top', horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.7),
                family='monospace')
    
    # Residuals for current best
    residuals_best = data - best_model
    ax_resid_best.plot(x, residuals_best, 'k.', alpha=0.4, markersize=2)
    ax_resid_best.axhline(0, color='b', linestyle='--', lw=1)
    if uncertainties is not None:
        ax_resid_best.fill_between(x, -uncertainties, uncertainties, 
                                   color='gray', alpha=0.2)
    ax_resid_best.set_ylabel('Residuals', fontsize=9)
    ax_resid_best.set_title('Current Best Residuals', fontsize=9)
    ax_resid_best.grid(True, alpha=0.3)
    ax_resid_best.tick_params(labelsize=8)
    
    # ========================================================================
    # Subsequent rows: Each candidate
    # ========================================================================
    for idx, cand in enumerate(candidates):
        ax_data = axes[idx + 1, 0]
        ax_resid = axes[idx + 1, 1]
        
        is_selected = (idx == best_candidate_idx)
        color = 'green' if is_selected else 'red'
        marker = '✓' if is_selected else '✗'
        
        # Plot data and model
        if uncertainties is not None:
            ax_data.errorbar(x, data, yerr=uncertainties, fmt='k.', alpha=0.3,
                           capsize=0, elinewidth=0.5, markersize=2)
        else:
            ax_data.plot(x, data, 'k.', alpha=0.3, markersize=2)
        
        cand_model = composite_model(x, cand['params'])
        ax_data.plot(x, cand_model, color=color, lw=2, alpha=0.9)
        
        # Get component parameters for the NEW component
        n_g_current = sum(1 for name in best_params if name.startswith('g') and name.endswith('_amp'))
        n_l_current = sum(1 for name in best_params if name.startswith('l') and name.endswith('_amp'))
        n_c_current = sum(1 for name in best_params if name.startswith('c'))
        
        if cand['type'] == 'gaussian':
            param_text_short = get_component_params_text(cand['params'], 'gaussian', n_g_current)
        elif cand['type'] == 'lorentzian':
            param_text_short = get_component_params_text(cand['params'], 'lorentzian', n_l_current)
        else:  # chebyshev
            param_text_short = get_component_params_text(cand['params'], 'chebyshev', n_c_current)
        
        # Title with fit statistics
        title = f"{marker} {cand['name']}\n{param_text_short}"
        ax_data.set_title(title, fontsize=9, fontweight='bold', color=color)
        ax_data.set_ylabel('Flux', fontsize=9)
        ax_data.grid(True, alpha=0.3)
        ax_data.tick_params(labelsize=8)
        
        # Add text box with fit statistics (left side)
        stats_text = (f"χ²ᵣ={cand['reduced_chisqr']:.3f}\n"
                     f"BIC={cand['bic']:.1f}\n"
                     f"AIC={cand['aic']:.1f}\n"
                     f"p={cand['p_value']:.1e}")
        ax_data.text(0.02, 0.98, stats_text, transform=ax_data.transAxes,
                    fontsize=7, verticalalignment='top',
                    bbox=dict(boxstyle='round', facecolor=color, alpha=0.2))
        
        # Add ALL parameters text box (right side)
        param_text_full = format_params_grouped(cand['params'])
        ax_data.text(0.98, 0.97, param_text_full, transform=ax_data.transAxes,
                    fontsize=6, verticalalignment='top', horizontalalignment='right',
                    bbox=dict(boxstyle='round', facecolor='lightyellow' if not is_selected else 'lightgreen', 
                             alpha=0.7),
                    family='monospace')
        
        # Plot residuals
        cand_residuals = data - cand_model
        ax_resid.plot(x, cand_residuals, 'k.', alpha=0.4, markersize=2)
        ax_resid.axhline(0, color=color, linestyle='--', lw=1)
        if uncertainties is not None:
            ax_resid.fill_between(x, -uncertainties, uncertainties,
                                 color='gray', alpha=0.2)
        ax_resid.set_ylabel('Residuals', fontsize=9)
        ax_resid.set_title(f'Residuals: {cand["name"]}', fontsize=9, color=color)
        ax_resid.grid(True, alpha=0.3)
        ax_resid.tick_params(labelsize=8)
    
    # Set x-labels on bottom row only
    axes[-1, 0].set_xlabel('Wavelength', fontsize=9)
    axes[-1, 1].set_xlabel('Wavelength', fontsize=9)
    
    plt.tight_layout()
    plt.savefig(f'iteration_{iteration:02d}.png', dpi=120, bbox_inches='tight')
    plt.show(block=False)
    plt.draw()
    
    # Wait for user input
    print("\n" + "="*70)
    print("Press ENTER to continue to next iteration...")
    print("="*70)
    input()
    plt.close(fig)


# ============================================================================
# PART 6: ITERATIVE ALGORITHM (WITH FIX FOR CHI-SQUARE!)
# ============================================================================

def iterative_fit(x, data, uncertainties=None, 
                  max_iterations=20, 
                  significance_level=0.05,
                  criterion='bic',
                  verbose=True,
                  plot=True,
                  true_params=None):
    """
    Iterative model building with forward selection.
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
        print(f"\n{'='*70}")
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
            print(f"\n{'='*70}")
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
        print(f"\n{'='*70}")
        print("COMPARISON: TRUE vs FITTED PARAMETERS")
        print(f"{'='*70}")
        
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
        print(f"\n{'='*70}")
        print(f"FINAL MODEL:")
        n_g = sum(1 for name in best_params if name.startswith('g') and name.endswith('_amp'))
        n_l = sum(1 for name in best_params if name.startswith('l') and name.endswith('_amp'))
        n_c = sum(1 for name in best_params if name.startswith('c'))
        print(f"  {n_g} Gaussian(s) + {n_l} Lorentzian(s) + Chebyshev order {n_c-1}")
        print(f"  Total parameters: {best_result.nvarys}")
        print(f"  Reduced χ²: {calculate_reduced_chi_square(best_result, n_data):.4f}")
        print(f"{'='*70}\n")
    
    return best_result, best_params, history


# ============================================================================
# PART 7: USAGE EXAMPLE WITH ABSORPTION FEATURES
# ============================================================================

if __name__ == "__main__":
    # Generate synthetic data for testing
    np.random.seed(42)
    x = np.linspace(5000, 6000, 500)
    
    # Polynomial continuum with 3rd order term (weak)
    continuum = 1.0 + 0.0001 * (x - 5500) - 2e-8 * (x - 5500)**2 + 5e-13 * (x - 5500)**3
    
    # Define 4 ABSORPTION features (now including two blended Gaussians)
    absorption1 = -0.15 * np.exp(-((x - 5250) / 25)**2)  # Gaussian absorption
    
    # Two blended Gaussians at ~5500, separated by ~half width (9 Angstroms)
    absorption2a = -0.12 * np.exp(-((x - 5475) / 18)**2)  # First Gaussian (weaker)
    absorption2b = -0.18 * np.exp(-((x - 5505) / 18)**2)  # Second Gaussian (stronger)
    
    absorption3 = -0.12 * 20**2 / ((x - 5750)**2 + 20**2)  # Lorentzian absorption
    
    # Define true parameters for comparison
    true_params = {
        'gaussians': [
            {'amp': -0.15, 'cen': 5250, 'wid': 25},
            {'amp': -0.12, 'cen': 5475, 'wid': 18},  # Blended pair
            {'amp': -0.18, 'cen': 5505, 'wid': 18}   # Blended pair
        ],
        'lorentzians': [
            {'amp': -0.12, 'cen': 5750, 'wid': 20}
        ],
        'chebyshev': []
    }
    
    # Build true model
    true_model = continuum + absorption1 + absorption2a + absorption2b + absorption3
    
    # Add noise
    noise_level = 0.005
    noise = np.random.normal(0, noise_level, len(x))
    data = true_model + noise
    uncertainties = np.ones_like(data) * noise_level
    
    print(f"Generated synthetic absorption spectrum:")
    print(f"  Wavelength range: {x.min():.1f} - {x.max():.1f}")
    print(f"  Number of points: {len(x)}")
    print(f"  Noise level: {noise_level:.4f} ({100*noise_level:.2f}%)")
    print(f"  True features: 3 Gaussian absorption lines (2 blended) + 1 Lorentzian")
    print(f"  Continuum: 3rd-order polynomial (will be fitted with Chebyshev)")
    
    # Run iterative fitting
    result, params, history = iterative_fit(
        x, data, uncertainties,
        max_iterations=15,
        significance_level=0.01,
        criterion='bic',
        verbose=True,
        plot=True,
        true_params=true_params
    )
    
    # Print fitting history
    import pandas as pd
    df = pd.DataFrame(history)
    print("\nFitting History:")
    print(df.to_string(index=False))
    
    # Final summary plot
    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
    
    axes[0].errorbar(x, data, yerr=uncertainties, fmt='k.', alpha=0.3, 
                     label='Data', capsize=0, markersize=2)
    axes[0].plot(x, true_model, 'b--', lw=2, label='True model', alpha=0.7)
    axes[0].plot(x, composite_model(x, params), 'r-', lw=2, label='Fitted model')
    axes[0].set_ylabel('Normalized Flux', fontsize=11)
    axes[0].legend(fontsize=10, loc='best')
    axes[0].set_title('ABSORPTION SPECTRUM: FINAL MODEL COMPARISON', fontsize=12, fontweight='bold')
    axes[0].grid(True, alpha=0.3)
    
    axes[1].plot(x, continuum, 'b--', lw=2, alpha=0.7, label='True continuum')
    
    fitted_continuum = np.zeros_like(x)
    i = 0
    x_norm = 2 * (x - x.min()) / (x.max() - x.min()) - 1
    while f'c{i}' in params:
        coeff_array = [0] * (i+1)
        coeff_array[i] = params[f'c{i}'].value
        fitted_continuum += chebyshev.chebval(x_norm, coeff_array)
        i += 1
    
    axes[1].plot(x, fitted_continuum, 'r-', lw=2, label='Fitted continuum (Chebyshev)')
    axes[1].set_ylabel('Continuum Flux', fontsize=11)
    axes[1].legend(fontsize=10, loc='best')
    axes[1].set_title('Continuum Comparison', fontsize=11)
    axes[1].grid(True, alpha=0.3)
    
    residuals = data - composite_model(x, params)
    axes[2].plot(x, residuals, 'k.', alpha=0.4, markersize=2)
    axes[2].axhline(0, color='r', linestyle='--', lw=2)
    axes[2].fill_between(x, -uncertainties, uncertainties, color='gray', alpha=0.2, label='±1σ')
    axes[2].set_xlabel('Wavelength (Å)', fontsize=11)
    axes[2].set_ylabel('Residuals', fontsize=11)
    axes[2].legend(fontsize=9)
    axes[2].grid(True, alpha=0.3)
    
    rms = np.std(residuals)
    axes[2].text(0.02, 0.95, f'RMS = {rms:.5f}\nExpected = {noise_level:.5f}',
                transform=axes[2].transAxes, verticalalignment='top', fontsize=9,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    plt.savefig('final_absorption_comparison.png', dpi=120)
    plt.show()
