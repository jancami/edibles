import numpy as np
import matplotlib.pylab as plt
import astropy
from astropy.stats import sigma_clip
from scipy import stats
from scipy.optimize import leastsq
import astropy.io.ascii as ascii


# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%
def cont_histo(flux, rms_noise):
    """
    Create histogram distribution of the flux data
    and select a narrower range around the maximum
    of the histogram distribution

    Parameters
    ----------
    flux : np.ndarray
        One-dimension array of flux values
    rms_noise : float
        The estimated RMS noise level of the data

    Returns
    -------
    all_bins : np.ndarray
        One-dimension array with the value of bins of the histogram
    all_hist : np.ndarray
        One-dimension array with the value of the position of the bins
    sel_bins : np.ndarray
        One-dimension array with the value of bins of the histogram
        for the selected bins around the maximum
    sel_hist : np.ndarray
        One-dimension array with the value of position of the bins
        for the selected bins around the maximum
    sel_flux : np.ndarray
        One-dimension array of the flux values selected
        around the maximum of the histogram
    """


    # creating a general histogram of the flux data
    # main variables are:
    #   all_hist     - counts in each bin of the histogram
    #   all_bins     - location of the bins (fluxes)
    #   all_number_* - index of the array
    # number_bins = int((np.amax(flux)-np.amin(flux))/(1*rms_noise))
    number_bins = 20
    all_hist, all_bin_edges = np.histogram(flux, number_bins)
    all_bins = all_bin_edges[0:len(all_bin_edges)-1]
    all_bins = [x + (all_bins[1]-all_bins[0])/2. for x in all_bins]
    all_number_max_array = (np.where(all_hist == all_hist.max())[0])
    all_number_max = all_number_max_array[0]
    all_bins_max = (all_bin_edges[all_number_max] + (all_bins[1]-all_bins[0])/2.)

    # Gaussian fit around the maximum of the distribution
    # determining the range to fit the Gaussian function
    all_number_left  = (np.where(((all_hist == 0) & (all_bins <= all_bins_max)) | (all_bins == all_bins[0]))[0]).max()
    all_number_right = (np.where(((all_hist == 0) & (all_bins >= all_bins_max)) | (all_bins == all_bins[number_bins-1]))[0]).min()
    all_number_total = abs(all_number_right-all_number_max)+abs(all_number_left-all_number_max)
    emission_absorption_ratio = abs(all_number_right-all_number_max)*1.0/(all_number_total*1.0)
    if (emission_absorption_ratio >= 0.66):
        lower_all_bins = all_bins_max - 8. * (all_bins[1]-all_bins[0])
        upper_all_bins = all_bins_max + 4. * (all_bins[1]-all_bins[0])
    if (emission_absorption_ratio <= 0.33):
        lower_all_bins = all_bins_max - 4. * (all_bins[1]-all_bins[0])
        upper_all_bins = all_bins_max + 8. * (all_bins[1]-all_bins[0])
    if ((emission_absorption_ratio > 0.33) and (emission_absorption_ratio < 0.66)):
        lower_all_bins = all_bins_max - 5. * (all_bins[1]-all_bins[0])
        upper_all_bins = all_bins_max + 5. * (all_bins[1]-all_bins[0])
    sel_bins_array = np.where((all_bins >= lower_all_bins) & (all_bins <= upper_all_bins))[0]
    if (len(sel_bins_array) < 3):
        sel_bins_array = [sel_bins_array[0]-2, sel_bins_array[0]-1, sel_bins_array[0], sel_bins_array[0]+1, sel_bins_array[0]+2]
        lower_all_bins = all_bins[sel_bins_array[0]]
        upper_all_bins = all_bins[sel_bins_array[len(sel_bins_array)-1]]
    sel_bins = all_bins[sel_bins_array[0]:sel_bins_array[len(sel_bins_array)-1]+1]
    sel_hist = all_hist[sel_bins_array[0]:sel_bins_array[len(sel_bins_array)-1]+1]
    sel_flux = flux[(flux >= lower_all_bins) & (flux <= upper_all_bins)]

    return all_bins, all_hist, sel_bins, sel_hist, sel_flux



# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%
def c_Gaussian(flux, rms_noise):
    """
    Perform Gaussian fit to the distribution of variable flux, and determine
    the center and width of the Gaussian. Similarly, perform the Gaussian
    fit to a selected range of the distribution around the maximum of
    the histogram, and determine the center and width of the new Gaussian

    Parameters
    ----------
    flux : np.ndarray
        One-dimension array of flux values
    rms_noise : float
        The estimated RMS noise level of the data

    Returns
    -------
    Gaussian_flux : float
    Gaussian_noise : float
        The measured continuum flux and estimated 1-sigma noise as the
        center and width of the Gaussian fit to the histogram distribution
        The estimated 1-sigma per-channel noise around that measurement
    GaussNw_flux : float
    GaussNw_noise : float
        The measured continuum flux and estimated 1-sigma noise as the
        center and width of the Gaussian fit to a selected range around
        the maximum of the distribution
    """

    fitfunc = lambda p, x: p[0]*np.exp(-0.5*((x-p[1])/p[2])**2.)
    errfunc = lambda p, x, y: (y - fitfunc(p, x))

    all_bins, all_hist, sel_bins, sel_hist, sel_flux = cont_histo(flux, rms_noise)

    meansel_flux = np.mean(sel_flux)
    meansel_sigma = np.std(sel_flux)

    init = [all_hist.max(), meansel_flux, meansel_sigma]
    out = leastsq(errfunc, init, args=(all_bins, all_hist))
    c = out[0]
    Gaussian_flux = c[1]
    Gaussian_noise = c[2]

    init = [all_hist.max(), meansel_flux, meansel_sigma]
    out = leastsq(errfunc, init, args=(sel_bins, sel_hist))
    d = out[0]
    GaussNw_flux = d[1]
    GaussNw_noise = d[2]

    return Gaussian_flux, Gaussian_noise, GaussNw_flux, GaussNw_noise



# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%
def find_flat_region(wave, xcont, ycont):
    '''
    This code find the attached points in a list
    '''
    diff = np.mean(np.diff(wave))
    epsilon = 0.001

    xneighbour = []
    yneighbour = []
    xtemp = []
    ytemp = []
    for loop in range(1,len(xcont)):
        d_temp = abs(xcont[loop] - xcont[loop-1])
        # find the last set
        if (loop == len(xcont)-1) & (len(xtemp) >= 5):
            xneighbour.append(xtemp)
            yneighbour.append(ytemp)
        else:
            # find the middle sets
            if d_temp <= diff+epsilon:
                xtemp.append(xcont[loop])
                ytemp.append(ycont[loop])
            else:
                if len(xtemp) >= 5:
                    xneighbour.append(xtemp)
                    yneighbour.append(ytemp)
                xtemp = []
                ytemp = []

    return xneighbour, yneighbour






# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%
def continuum_sigma(wave, flux, plot=None):

    # derivative the data
    dydx = np.diff(flux)/np.diff(wave)
    wave1 = wave[1:]
    flux1 = dydx

    # wave1 = wave
    # flux1 = flux


    # find the continuum and white noise
    Gaussian_flux, Gaussian_noise, GaussNw_flux, GaussNw_noise = c_Gaussian(flux1, 1)
    if GaussNw_flux > 1.1: GaussNw_flux = np.mean(flux1)

    # print GaussNw_flux

    # find all flux values in this range
    idx = [i for i,x in enumerate(flux1) if (x <= GaussNw_flux + GaussNw_noise) & (x >= GaussNw_flux - GaussNw_noise)]
    xcont = wave1[idx]
    ycont = flux1[idx]


    # plt.plot(wave1, flux1)
    # plt.show()

    # find flat region
    xn, yn = find_flat_region(wave, xcont, ycont)

    # return the longest flat region
    xn_len = []
    for loop_a in range(len(xn)): xn_len.append(len(xn[loop_a]))
    idx_max = [i for i, e in enumerate(xn_len) if e == len(max(xn,key=len))]

    # find the flat region in main flux
    idx = [i for i, e in enumerate(wave) if e in set(xn[idx_max[0]])]
    wave_flat = wave[idx]
    flux_flat = flux[idx]


    # plot
    if plot is True:
        plt.plot(wave, flux, color='black')
        plt.plot(wave_flat, flux_flat, marker='+', linestyle='none', color='red')
        plt.show()

    return np.std(flux_flat)
