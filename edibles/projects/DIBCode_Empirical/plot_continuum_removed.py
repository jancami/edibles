import numpy as np
import matplotlib.pyplot as plt
def plot_continuum_removed(spectra_data,figpath):
    """
      Divides out the Chebyshev continuum.
      Plots normalized and scaled spectra between 0.5 and 1
      Args:
          spectra_data: data of the spectra received from plotall()
          figpath: path for the final plotted figure to go into as a .png file
      Returns: Plot of the normalized and scaled spectra.
      """
    fig1, (ax1_obs, ax1_mod) = plt.subplots(1, 2, figsize=(16, 6), sharex=True, sharey=True)
    fig2, (ax2_obs, ax2_mod) = plt.subplots(1, 2, figsize=(16, 6), sharex=True, sharey=True)

    colors = plt.cm.viridis(np.linspace(0, 1, len(spectra_data)))

    for i, spec in enumerate(spectra_data):
        c = colors[i]
        fd = spec['fit_dict']
        y_obs_norm = fd['y_obs'] / fd['y_cont']
        y_model_norm = fd['y_total_model'] / fd['y_cont']
        m_min = np.min(y_model_norm)

        def scale_norm(arr, model_min):
            """
            Normalized the values according to the model output.
            Args:
                arr: Series of values in an array
                model_min: Minimum value of the model function
            Returns:
                Array values rescaled so that 1 is the maximum of the model and 0.5 is the minimum.
            """
            if model_min >= 1.0: return arr
            return 1.0 - ((1.0 - arr) / (1.0 - model_min)) * 0.5

        y_obs_scaled = scale_norm(y_obs_norm, m_min)
        y_model_scaled = scale_norm(y_model_norm, m_min)

        # Normalized
        ax1_obs.plot(spec['x'], y_obs_norm, color=c, alpha=0.3, lw=0.7)
        ax1_mod.plot(spec['x'], y_model_norm, color=c, alpha=0.9, label=spec['label'])

        # Normalized and Scaled (1.0 to 0.5)
        ax2_obs.plot(spec['x'], y_obs_scaled, color=c, alpha=0.3, lw=0.7)
        ax2_mod.plot(spec['x'], y_model_scaled, color=c, alpha=0.9, label=spec['label'])

    fig1.suptitle('Continuum Removed (Baseline = 1.0)', fontsize=16)
    ax1_obs.set_title('Normalized OBSERVED')
    ax1_mod.set_title('Normalized MODELS')
    fig2.suptitle('Continuum Removed & Scaled (Baseline 1.0, Min 0.5)', fontsize=16)
    ax2_obs.set_title('Scaled OBSERVED')
    ax2_mod.set_title('Scaled MODELS')
    for ax in [ax1_obs, ax1_mod, ax2_obs, ax2_mod]:
        ax.axhline(1.0, color='red', linestyle='--', alpha=0.5)
        ax.set_xlabel('Wavelength (Å)')
        ax.set_ylabel('Intensity')

    #ax1_mod.legend(loc='lower left', fontsize='xx-small', ncol=3)
    #ax2_mod.legend(loc='lower left', fontsize='xx-small', ncol=3)

    fig1.tight_layout()
    fig2.tight_layout()
    plt.savefig(f"{figpath}/Continuum_Removed_Scaled.png")
    plt.show()
