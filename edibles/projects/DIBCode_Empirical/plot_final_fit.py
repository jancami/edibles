import matplotlib.pyplot as plt
def plot_final_fit(fd, target, figpath):
    """
    Plots the final continuum fit, components, and residuals
    Args:
        fd: Dictionary of components for the best fit
        target: Target sightline being observed
    Returns:
        best_fit: A list of best fit parameters
        - Plot of final continuum fit
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 9), sharex=True, gridspec_kw={'height_ratios': [3, 1]})

    ax1.plot(fd['x_corr'], fd['y_obs'], color='silver', alpha=0.7, label='Corrected Data', zorder=1)
    ax1.plot(fd['x_corr'], fd['y_cont'], color='red', lw=2, label=f'Chebyshev Continuum (Order {fd["order"]})',
             zorder=2)
    ax1.plot(fd['x_corr'], fd['y_total_model'], color='black', lw=1.5, label='Total Fit (Poly + Scaled Comps)',
             zorder=5)

    for i, profile in enumerate(fd['scaled_individual_comps']):
        ax1.plot(fd['x_corr'], fd['y_cont'] + profile, label=f"Comp: {fd['comp_ids'][i]}", linestyle='--', alpha=0.8,
                 zorder=4)

    ax1.set_title(f"Spectral Fit: BIC Optimization & Component Visualization for {target}")
    ax1.set_ylabel("Intensity")
    ax1.legend(loc='lower left', fontsize='small', ncol=2)

    ax2.plot(fd['x_corr'], fd['residuals'], color='gray', lw=0.8)
    ax2.axhline(0, color='red', linestyle='--', alpha=0.6)
    ax2.set_ylabel("Residuals")
    ax2.set_xlabel("Wavelength (Velocity Corrected)")

    plt.tight_layout()
    plt.savefig(f"{figpath}/Component_Visualization_{target}")
    plt.show()
