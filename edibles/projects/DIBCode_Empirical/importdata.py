import numpy as np

from edibles.utils.edibles_oracle import EdiblesOracle
from edibles.utils.edibles_spectrum import EdiblesSpectrum
from edibles.projects.DIBCode_Empirical.SNR import SNR
def importdata(target, minrange, maxrange, data_piece, ContinuumMin, ContinuumMax):
    """
    Imports the data from the EDIBLES dataset.
    Coadds the data so that it can be used in all other parts of the software.
    Args:
        target: Target sightline being observed
        minrange: Minimum wavelength in our final model.
        maxrange: Maximum wavelength in our final model.
        data_piece: String containing the path and order of our target, i.e. 564nm_redu_O9
        ContinuumMin: Minimum wavelength of the continuum range we use for our signal-to-noise measurement.
        ContinuumMax: Maximum wavelength of the continuum range we use for our signal-to-noise measurement.
    Returns:
        common_wave[mask]: The filtered version of our wavelength range, (only the region for the model fit)
        coadded_flux[mask]: The filtered version of our coadded specta, (only the region for the model fit)
        coadd_SNR: The signal-to-noise ratio of the continuum window
        target: The target sightline for our observations.
    """
    oracle = EdiblesOracle()
    all_specs = list(oracle.getObsListByTarget(target=target, MergedOnly=False))

    target_specs = [f for f in all_specs if data_piece in f]

    print(f"Found {len(target_specs)} matching files.")

    if len(target_specs) == 0:
        raise RuntimeError("No matching files found.")

    DR_interp_flux_list = []
    DR_SNR_list = []

    common_wave = np.linspace(minrange - 20, maxrange + 20, 1000)

    for fname in target_specs:
        print("Loading:", fname)
        edibles_spec = EdiblesSpectrum(fname)

        wave = edibles_spec.bary_wave
        flux = edibles_spec.flux

        norm_flux = flux / np.median(flux)

        snr = SNR(wave, norm_flux, ContinuumMin, ContinuumMax)

        interp_flux = np.interp(common_wave, wave, norm_flux)

        DR_interp_flux_list.append(interp_flux)
        DR_SNR_list.append(snr)

    weights = np.array(DR_SNR_list)
    weights /= np.sum(weights)

    coadded_flux = np.sum([w * f for w, f in zip(weights, DR_interp_flux_list)], axis=0)

    coadd_SNR = SNR(common_wave, coadded_flux, ContinuumMin, ContinuumMax)

    print("Coadded SNR:", coadd_SNR)
    mask = (common_wave >= minrange) & (common_wave <= maxrange)

    return common_wave[mask], coadded_flux[mask], coadd_SNR, target, common_wave, coadded_flux
