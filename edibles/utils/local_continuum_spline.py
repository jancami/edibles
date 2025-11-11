import numpy as np
import matplotlib.pyplot as plt


# Function to fit the local continuum and normalize the spectrum
def local_continuum(wave, flux, positions=None, windows=1.0, spline_order=1, silent=True):
    """A function that will fit a local continuum spline to a "spectrum" using a
    list of anchor points. Each anchor point has a "continuum" window (or the
    same single value for all anchor points). A spline (order s) is fit to the
    "continuum" data points. The continumm is then created on the input
    "wavelength" grid. The input data tuple (wave, flux) is then normalised
    giving (wave, normalised_flux).

    Args:
        data (tuple): In the form (wave, flux)
        positions (list): List of spline anchor points
        windows (list): List of window sizes around each anchor point (in Angstrom)
        spline_order (int): Order (s) of spline fit
        silent (bool): If true, no plots will generate. Default: True

    Returns:
        tuple: tuple containing:

            tuple: Continuum flux in the form (flux) on the same spectral grid as (wave)

            tuple: Normalised input spectrum In the form (normalised_flux)

    """
    if len(windows) == 1:
        win = np.ones(len(positions)) * windows[0]
    elif len(windows) != len(positions):
        raise ValueError("windows needs to be length 1 or same length as positions (number of anchor points)")
    else:
        win = windows

    idx_continuum = []
    for anchor, window in zip(positions, win):
        idx_c = np.where((wave > anchor - window / 2.0) & (wave < anchor + window / 2.0))[0]  # Extract the array from the tuple
        idx_continuum.extend(idx_c)

    idx_continuum = np.array(idx_continuum)
    cont_wave = wave[idx_continuum]
    cont_flux = flux[idx_continuum]

    c1 = np.polyfit(cont_wave, cont_flux, spline_order)
    p1 = np.poly1d(c1)
    continuum = p1(cont_wave)

    flux_std = np.std(cont_flux-continuum)

    if not silent:
        plt.plot(wave, flux, label='Original Flux')
        plt.plot(cont_wave, continuum, label='Fitted Continuum')
        plt.legend()
        plt.show()

    return continuum, flux_std