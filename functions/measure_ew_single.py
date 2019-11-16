import numpy as np
from scipy.interpolate import CubicSpline
# import matplotlib.pyplot as plt

def measure_ew_single(data, borders=None, silent=True, anchors=None, windows=1.0, spline_order=1):
    # should anchors, windows, spline_order be in **kargs arugment?

    from edibles import local_continuum_spline

    '''A function that will calculate the Equivalent Width (EW) using min,max integration borders.

    A spline can be fit again to a list of continuum anchor points. If anchors==None we take continuum = unity (1).

    INPUT:
    data:         [tuple]             In the form (wave, flux)
    anchors:      []                  List of spline anchor points
    windows:      []                  List of window sizes around each anchor point (in Angstrom)
    borders:      [real,real]         2-element list with the minimum and maximum integration borders: (x1,x2)
    spline_order: [int - 1]           Order (s) of spline fit
    silent:       [bool - True]       If true, no plots will generate

    OUTPUT:
    ew:           []                  equivalent width (direct integration)
    ew_sigma:     []                  equivalent width uncertainty. SNR is computed from the anchor windows
    '''

    wave, flux = data

    # check that x1 < x2, and that x1, x2 are within the data x-range
    if (x2 >= x1):
        print("error borders: x1 needs to be smaller than x2")

    if ( (x1 < wave[0]) | (x1 > wave[-1]) ) :
        print("border x1 falls outisde input range")

    if ( (x2 < wave[0]) | (x2 > wave[-1]) ) :
        print("border x2 falls outisde input range")

    if anchors not None:
        # if list of anchor positions is provided derive a new local cotinuum
        # find the local continuum (default is straight line)
        normalised, continuum = local_continuum_spline(data, positions=anchors, windows=windows, spline_order=spline_order )

    if anchors is None:
        # if no list of anchor positions is provided, assume continuum is 1.0
        continuum = wave*0.0+1.0

    # now integrate between continuum and spectrum from x1 to x2, divide by (x2-x1)
    idx = np.where( ( wave > x1) * (wave < x2) )
    integrated_flux = np.integrate(continuum[idx]-flux[idx])
    ew = integrated_flux / np.abs(x2-x1)

    # compute the SNR from anchor windows
    if len(windows) == 1:
        # replicate windows to same length as positions; creat n-element array with single value in each element
        win = np.ones(len(positions))
        win = win*windows[0]

    # error if now 'windows' list not has same length as 'positions' list
    if len(win) != len(positions):
        print("error -- windows needs to be length 1 or same length as positions (nr achnor points)")

    # use list of anchor points and window sizes to sub-select the continuum points
    idx_continuum = []
    for anchor,windows in zip(positions,win):
        idx_c = np.where( (wave > anchor-win/2.) & (wave < anchor+win/2.) )
        idx_continuum = idx_continuum + idx_c

    cont_wave = wave[idx_continuum]
    cont_flux = flux[idx_continuum]

    # compute the stddev from all continuum windows
    sigma = np.nanstd(cont_flux)

    DELTA_lambda = borders[1]-borders[0]
    delta_lambda = DELTA_lambda/len(flux[idx])

    ew_sigma = UncertaintyEW(delta_lambda=delta_lambda, DELTA_lambda=DELTA_lambda, SNR=(1.0/sigma))

    if silent is False:
        print("Equivalent width ", ew, " +- ", ew_sigma*1000.0, " mA")

    return ew, ew_sigma


def UncertaintyEW(delta_lambda=1, DELTA_lambda=1, SNR=1):
    """ sigma_EW^2 = sigma_noise^2 + sigma_continuum^2
        From Vos et al. 2011, A&A 533, A129 (Appendix A), based on Chalabaev & Maillard (1983)
        Follows Vollmann & Eversberg (2006) in weak line limit:
           sigma_W = sqrt(2) / SNR * sqrt(delta_lambda * Delta_lambda)
        Is upper limit for strong narrow lines such as KI.

        DELTA_lambda = integration range of line and continuum (~equal)
        delta_lambda = spectral dispersion in \AA/pixel.
        SNR = signal-to-noise ratio per pixel = 1/stdev
    """

    ew_sigma =  np.sqrt(2.0)/SNR * np.sqrt(DELTA_lambda * delta_lambda) #

    return ew_sigma
