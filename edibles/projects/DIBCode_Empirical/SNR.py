import numpy as np
def SNR(wave, flux, min_wave, max_wave):
    """Calculates signal-to-noise ratio in the selected continuum region
    Args:
        wave:  Array of wavelength values as a NumPy array.
        flux: Array of flux values.
        min_wave: Minimum wavelength in our final model.
        max_wave: Maximum wavelength in our final model.
    Returns:
        Signal-to-noise raito value.
        """
    mask = (wave >= min_wave) & (wave <= max_wave)
    wave_masked = wave[mask]
    flux_masked = flux[mask]

    coeffs = np.polyfit(wave_masked, flux_masked, 1)
    fit = np.polyval(coeffs, wave_masked)

    return 1 / np.std(flux_masked / fit)
