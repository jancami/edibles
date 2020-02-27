import numpy as np
import matplotlib.pyplot as plt


def local_continuum(data, positions=None, windows=1.0, spline_order=1, silent=True):
    """A function that will fit a local continuum spline to a "spectrum" using a
    list of anchor points. Each anchor point has a "continuum" window (or the
    same single value for all anchor points). A spline (order s) is fit to the
    "continuum" data points. The continumm is then created on the input
    "wavelength" grid. The input data tuple (wave, flux) is then normalised
    giving (wave, normalised_flux).

    INPUT:
    data:         [tuple]             In the form (wave, flux)
    positions:    []                  List of spline anchor points
    windows:      []                  List of window sizes around each anchor point (in Angstrom)
    spline_order: [int - 1]           Order (s) of spline fit
    silent:       [bool - True]       If true, no plots will generate

    OUTPUT:
    continuum:    [tuple]             Continuum flux in the form (flux) on the same
                                      spectral grid as (wave)
    normalised:   [tuple]             Normalised input spectrum In the form (normalised_flux)
    ##c1:         [?]                 polyfit_result - can be used by np.poly1d(polyfit_result)
                                      to create continuum again: continuum = p1(wave_grid)
    """

    wave, flux = data

    if len(windows) == 1:
        # replicate windows to same length as positions; creat n-element array
        # with single value in each element
        win = np.ones(len(positions))
        win = win * windows[0]

    # error if now 'windows' list not has same length as 'positions' list
    if len(win) != len(positions):
        print(
            "error -- windows needs to be length 1 or same length as \
            positions (nr achnor points)"
        )

    # use list of anchor points and window sizes to
    # sub-select the continuum points
    idx_continuum = []
    for anchor, windows in zip(positions, win):
        idx_c = np.where((wave > anchor - win / 2.0) & (wave < anchor + win / 2.0))
        idx_continuum = idx_continuum + idx_c

    cont_wave = wave[idx_continuum]
    cont_flux = flux[idx_continuum]

    c1 = np.polyfit(cont_wave, cont_flux, spline_order)
    p1 = np.poly1d(c1)
    continuum = p1(wave)

    normalised_flux = flux / continuum

    if silent is False:
        plt.plot(wave, flux)
        plt.plot(wave, continuum)
        plt.show()

    return normalised_flux, continuum  # , c1
