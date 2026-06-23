import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from edibles.projects.DIBCode_Empirical.perform_fit import perform_fit
from edibles.projects.DIBCode_Empirical.importdata import importdata
def plotall(targets, minrange, maxrange, csv_name, c0, c1, model_v_shift, data_piece, ContinuumMin, ContinuumMax):
    '''Plots the model and data residuals for all targets.
    Args:
        targets: List of sightlines.
        minrange: Lowest wavelength of the model plot.
        maxrange: Highest wavelength of the model plot.
        csv_name: Location of the csv in the folder.
        c0: First-order Chebyeshev polynomial parameter.
        c1: Second-order Chebyeshev polynomial parameter.
        model_v_shift: Velocity shift of the model.
        data_piece: String containing the path and order of our target, i.e. 564nm_redu_O9
        ContinuumMin: Minimum wavelength of the continuum range we use for our signal-to-noise measurement.
        ContinuumMax: Maximum wavelength of the continuum range we use for our signal-to-noise measurement.
    Returns:
        Plot of all model components and data residuals for all target sightlines.
        '''
    spectra_data = []

    for i in targets:
        recalibrated_wavelength, coadded_flux, coadd_SNR, target, common_wave_full, coadded_flux_full = importdata(target=i, minrange=minrange,
                                                                              maxrange=maxrange, data_piece=data_piece,
                                                                              ContinuumMin=ContinuumMin,
                                                                              ContinuumMax=ContinuumMax)

        #plt.ioff()
        df = pd.read_csv(csv_name)
        components_params = dict(zip(df.iloc[:, 0], df.iloc[:, 1]))

        fd = perform_fit(recalibrated_wavelength, coadded_flux, components_params,
                         c0=c0, c1=c1,
                         target=target, model_v_shift=model_v_shift)

        if fd:
            y_m, y_o = fd['y_total_model'], fd['y_obs']

            def scale_to_range(arr):
                """
                Rescales the model flux values.
                Arg:
                    arr: Series of values in an array.
                Returns:
                    Array values rescaled into the range 0.5 to 1
                """
                a_min, a_max = np.min(arr), np.max(arr)
                return ((arr - a_min) / (a_max - a_min)) * 0.5 + 0.5

            spectra_data.append({
                'label': target,
                'x': fd['x_corr'],
                'm_raw': y_m,
                'o_raw': y_o,
                'resid': y_o - y_m,  # Calculate residual (Observed - Model)
                'm_scaled': scale_to_range(y_m),
                'o_scaled': scale_to_range(y_o),
                'fit_dict': fd
            })
        plt.close('all')

    fig = plt.figure(figsize=(16, 14))
    gs = fig.add_gridspec(3, 2)

    ax_m_raw = fig.add_subplot(gs[0, 0])
    ax_o_raw = fig.add_subplot(gs[0, 1], sharex=ax_m_raw)
    ax_m_scaled = fig.add_subplot(gs[1, 0], sharex=ax_m_raw)
    ax_o_scaled = fig.add_subplot(gs[1, 1], sharex=ax_m_raw)
    # This axis spans both columns for a wide comparison
    ax_resid = fig.add_subplot(gs[2, :], sharex=ax_m_raw)

    colors = plt.cm.viridis(np.linspace(0, 1, len(spectra_data)))

    for i, spec in enumerate(spectra_data):
        c = colors[i]
        ax_m_raw.plot(spec['x'], spec['m_raw'], color=c, alpha=0.8, label=spec['label'])
        ax_o_raw.plot(spec['x'], spec['o_raw'], color=c, alpha=0.5, lw=0.5)
        ax_m_scaled.plot(spec['x'], spec['m_scaled'], color=c, alpha=0.8)
        ax_o_scaled.plot(spec['x'], spec['o_scaled'], color=c, alpha=0.5, lw=0.5)

        ax_resid.plot(spec['x'], spec['resid'], color=c, alpha=0.7, lw=0.8)

    ax_m_raw.set_title('Unstretched MODELS')
    ax_o_raw.set_title('Unstretched OBSERVED DATA')
    ax_m_scaled.set_title('MODELS (0.5 to 1.0)')
    ax_o_scaled.set_title('OBSERVED DATA (0.5 to 1.0)')
    ax_resid.set_title('RESIDUALS (Observed - Model)')

    for ax in [ax_m_raw, ax_o_raw, ax_m_scaled, ax_o_scaled, ax_resid]:
        ax.grid(True, alpha=0.2)
        ax.set_ylabel('Intensity')

    ax_resid.axhline(0, color='black', linestyle='--', alpha=0.5)  # Zero line for residuals
    ax_resid.set_xlabel('Wavelength (Å)')
    #ax_m_raw.legend(loc='upper right', fontsize='xx-small', ncol=3)

    plt.tight_layout()
    plt.show()

    return spectra_data
