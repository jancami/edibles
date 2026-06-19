import numpy as np
import matplotlib.pyplot as plt
def plot_alignment_verification(x_ref, y_obs, best_y_data, model_template, v_grid, correlations, best_v, max_r, mask,
                                use_mask, target):
    """
    Plots the cross-correlation results and shifted spectra
    Args:
        x_ref: Initial independent variable data points (wavelengths) as a NumPy array.
        y_obs: Observed flux values from the sightline.
        best_y_data: The best profile for the DIB fit.
        model_template: Sum of the observed continuum and the initial model fit.
        v_grid: List of possible velocity values to run correlations over.
        correlations: List of correlation values for each velocity.
        best_v: Velocity that yields the best correlation
        max_r: Maximum correlation value.
        use_mask: Boolean - whether or not the mask is used for the fit
        target: Target sightline being observed.
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), gridspec_kw={'height_ratios': [2, 1]})

    ax1.plot(x_ref, model_template, color='black', lw=2, label='Reference Model', zorder=3)
    ax1.plot(x_ref, y_obs, color='gray', alpha=0.3, linestyle='--', label='Original Data (0 km/s)', zorder=1)
    ax1.plot(x_ref, best_y_data, color='royalblue', alpha=0.8, label=f'Shifted Data ({best_v:.2f} km/s)', zorder=2)

    if use_mask:
        ax1.fill_between(x_ref, 0, 1, where=mask, color='salmon', alpha=0.2,
                         label='Velocity Mask (y)', transform=ax1.get_xaxis_transform())

    ax1.set_title(f"Velocity Alignment Verification: {target}")
    ax1.set_ylabel("Normalized Intensity")
    ax1.legend(loc='upper right', fontsize='small')
    ax1.set_ylim(np.min(y_obs) * 0.95, np.max(y_obs) * 1.05)

    ax2.plot(v_grid, correlations, color='forestgreen', lw=1.5)
    ax2.axvline(best_v, color='red', linestyle=':', label=f'Peak: {best_v:.2f} km/s (r={max_r:.4f})')
    ax2.set_title(f"Correlation Coefficient vs. Velocity Shift {target}")
    ax2.set_xlabel("Velocity Shift (km/s)")
    ax2.set_ylabel("Correlation (r)")
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    plt.tight_layout()
    plt.show()

