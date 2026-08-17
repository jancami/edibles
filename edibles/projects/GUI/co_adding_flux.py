"""Co-adds a list of spectra into one.

Just numpy — no file I/O, no edibles dependency. You give it a list of
(wave, flux) pairs (the GUI gets these from flux_wave_find.wave_flux_data),
it interpolates everything onto a common wavelength grid (taken from
wave_list[ref_index] by default) and returns one averaged spectrum plus its
uncertainty. The "Co-add Selected" button in the GUI calls this right before
handing the result to the fit.
"""

import numpy as np


def coadd_spectra(wave_list, flux_list, sigma_list=None, ref_index=0):
    """
    Co-add multiple spectra with (optional) constant uncertainties per spectrum.

    Parameters
    ----------
    wave_list : list of 1D arrays
        [wave0, wave1, wave2, ...] one array per spectrum.
    flux_list : list of 1D arrays
        [flux0, flux1, flux2, ...] matching wave_list.
    sigma_list : list or 1D array or float, optional
        If given: one constant sigma per spectrum, e.g. [0.002, 0.003, 0.0015, ...].
        If a single float is given, the same sigma is used for all spectra.
    ref_index : int, optional
        Index of the spectrum whose wavelength grid will be used as the common grid.

    Returns
    -------
    wave_common : 1D array
        Common wavelength grid (same as wave_list[ref_index]).
    flux_coadd : 1D array
        Co-added flux on wave_common.
    sigma_coadd : float or 1D array
        Uncertainty of the co-added spectrum.
        - If sigma_list is provided and constant per spectrum: returns a single float.
        - If sigma_list is None: returns per-pixel std / sqrt(N).
    """

    # Convert to numpy arrays for safety
    wave_list = [np.array(w) for w in wave_list]
    flux_list = [np.array(f) for f in flux_list]

    n_spec = len(wave_list)
    if len(flux_list) != n_spec:
        raise ValueError("wave_list and flux_list must have the same length")

    # Common wavelength grid = chosen reference spectrum
    wave_common = wave_list[ref_index]

    # Interpolate all fluxes onto the common grid
    interp_fluxes = []
    for w, f in zip(wave_list, flux_list):
        f_interp = np.interp(wave_common, w, f)
        interp_fluxes.append(f_interp)

    interp_fluxes = np.array(interp_fluxes)  # shape: (n_spec, n_pixels)

    # If no sigmas given: simple average
    if sigma_list is None:
        flux_coadd = np.mean(interp_fluxes, axis=0)
        # estimate error from scatter: std / sqrt(N)
        sigma_coadd = np.std(interp_fluxes, axis=0, ddof=1) / np.sqrt(n_spec)
        return wave_common, flux_coadd, sigma_coadd

    # Handle sigma_list (constant per spectrum)
    sigma_list = np.array(sigma_list, dtype=float)
    if sigma_list.size == 1:
        # broadcast single sigma to all spectra
        sigma_list = np.full(n_spec, sigma_list.item())

    if sigma_list.size != n_spec:
        raise ValueError("sigma_list must have length 1 or match number of spectra")

    # Weights = 1 / sigma^2 (constants)
    weights = 1.0 / (sigma_list**2)          # shape: (n_spec,)
    Wsum = np.sum(weights)

    # Weighted sum over spectra: sum_i w_i * F_i(lambda)
    # reshape weights to (n_spec, 1) to broadcast over pixels
    flux_coadd = np.sum(interp_fluxes * weights[:, None], axis=0) / Wsum

    # Co-added sigma is the same for all pixels since sigmas are constant per spectrum
    sigma_coadd = 1.0 / np.sqrt(Wsum)

    return wave_common, flux_coadd, sigma_coadd
