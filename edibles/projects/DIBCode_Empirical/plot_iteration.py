import matplotlib.pyplot as plt
def plot_iteration(x, data, best_params, best_result, candidates, best_candidate_idx,
                   iteration, n_data, uncertainties=None):
    """Plots each candidate fitting option in a separate subplot.

    Renders a vertical grid of plots comparing the current running best model
    fit against a list of newly computed component candidates. Each row displays
    the data and model fit on the left, alongside the generated residuals on
    the right.

    Args:
        x: Independent variable data points (wavelengths) as a NumPy array.
        data: Observed dependent variable data points.
        best_params: The current authoritative lmfit Parameters before the loop.
        best_result: Minimizer result object tied to the best_params state.
        candidates: A list of dictionaries, where each dict represents a
          potential addition and contains keys: 'params', 'name', 'type',
          'reduced_chisqr', 'bic', 'aic', and 'p_value'.
        best_candidate_idx: List index of the winning candidate in this cycle.
        iteration: The numerical sequence loop index for titling and file saves.
        n_data: Number of data points used in evaluating degrees of freedom.
        uncertainties: Optional array of uncertainties (sigma) for the data
          points.
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
    ax_best.set_title(f'CURRENT BEST: {n_g}G + {n_l}L + Cheb({n_c - 1})',
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
    print("\n" + "=" * 70)
    print("Press ENTER to continue to next iteration...")
    print("=" * 70)
    input()
    plt.close(fig)
