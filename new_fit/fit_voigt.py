import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.signal import find_peaks

import edibles.fit.avoigt as fit
import edibles.fit.make_grid as mg
from edibles.new_fit.cont_model import Cont1D
from edibles.new_fit.v_model import Voigt1D
from edibles.new_fit.astro_v_model import AstroVoigt1D
from edibles.functions.continuum_guess import generate_continuum
from edibles.functions.line_merger import line_merger

from sherpa.data import Data1D
from sherpa.plot import DataPlot
from sherpa.plot import ModelPlot
from sherpa.stats import LeastSq
from sherpa.optmethods import LevMar
from sherpa.fit import Fit
from sherpa.plot import FitPlot


def fitter(wave_subset, flux_subset, peak_cutoff=0.5, n_points=5, b_eff=3.47, Gamma=6.0e07, scaling=1.0):
    '''
    Function that fits a spline continum and multiple voigt profile peaks simultaneously. 



    INPUT:
        file:         [str]    /path/to/file
        xmin:         [float]  cropped dataset minimum
        xmax:         [float]  cropped dataset maximum
        peak_cutoff:  [float]  threshold to find peaks, be wary with small values!
        n_points      [int]    number of continuum points, must be 2 < n_points < 8
        alpha:        [float]  Gaussian HWHM component
        gamma:        [float]  Lorentzian HWHM component
        scaling:      [float]  Height of the peak scaling factor


    OUTPUT:
        outputs plots and values

        currently returns None

    USAGE: 
        fitter(file, xmin=xmin, xmax=xmax, peak_cutoff=peak_cutoff, n_points=n_points, 
                alpha=alpha, gamma=gamma, scaling=scaling)

    '''

    n_piece = n_points - 1


    # hdu = fits.open(file)

    # spec_flux = hdu[0].data
    # crval1 = hdu[0].header["CRVAL1"]
    # cdelt1 = hdu[0].header["CDELT1"]
    # nwave = len(spec_flux)
    # wave = np.arange(0, nwave, 1)
    # spec_wave = (wave) * cdelt1 + crval1

    # plt.plot(spec_wave, spec_flux)
    # plt.show()

    # # create data subset
    # min_idx = (np.abs(spec_wave - xmin)).argmin()
    # max_idx = (np.abs(spec_wave - xmax)).argmin()
    # wave_subset = spec_wave[min_idx:max_idx]
    # flux_subset = spec_flux[min_idx:max_idx]

    # =========================================================================

    # CREATE MODEL
    model = 0


    # create initial continuum guess spline points
    y_spline, y_points= generate_continuum((wave_subset, flux_subset), 
                                            delta_v=1000, n_piece=n_piece)

    # create cont object and define parameters
    cont = Cont1D()

    # always at least 2 points / 1 piece
    if n_points >= 1:
        cont.y1            = y_points[0]
        cont.y1.frozen     = False
    if n_points >= 2:
        cont.y2            = y_points[1]
        cont.y2.frozen     = False
    if n_points >= 3:
        cont.y3            = y_points[2]
        cont.y3.frozen     = False
    if n_points >= 4:
        cont.y4            = y_points[3]
        cont.y4.frozen     = False
    if n_points >= 5:
        cont.y5            = y_points[4]
        cont.y5.frozen     = False
    if n_points >= 6:
        cont.y6            = y_points[5]
        cont.y6.frozen     = False
    if n_points >= 7:
        cont.y7            = y_points[6]
        cont.y7.frozen     = False
    if n_points >= 8:
        cont.y8            = y_points[7]
        cont.y8.frozen     = False

    cont.n_piece           = n_piece
    print(cont)

    # add cont to model
    model = cont


    # =========================================================================

    # find peaks of spectrum

    prominence = (np.max(flux_subset) - np.min(flux_subset)) * peak_cutoff
    peaks, _ = find_peaks(-flux_subset, prominence=prominence)

    plt.plot(wave_subset, flux_subset)
    plt.plot(wave_subset[peaks], flux_subset[peaks], 'x')
    plt.show()


    # =========================================================================

    # Semi-automated voigt profile model generator

    for i in range(len(peaks)):

        # create temporary voigt object
        obj = AstroVoigt1D()
        obj.cent          = wave_subset[peaks[i]]  #  7664.87
# <<<<<<< HEAD
        # obj.cent.frozen = False
        # obj.alpha         = alpha
        # obj.gamma         = gamma
# =======
        # obj.cent.frozen = False
        obj.b_eff         = b_eff
        obj.Gamma         = Gamma
# >>>>>>> 0bd5ecd769e2776ec33c0290a06cd233fa08264e
        obj.scaling       = scaling
        obj.scaling.frozen = False
        print(obj)

        # add object to model
        model += obj

        # reset obj variable
        obj = 0


    d = Data1D('Initial model', wave_subset, flux_subset)

    dplot = DataPlot()
    dplot.prepare(d)
    dplot.plot()

    mplot = ModelPlot()
    mplot.prepare(d, model)
    dplot.plot()
    mplot.overplot()
    # plt.show()

    stat = LeastSq()
    opt = LevMar()
    print(opt)

    vfit = Fit(d, model, stat=stat, method=opt)
    print(vfit)
    vres = vfit.fit()
    print()
    print()
    print('Did the fit succeed? [bool]')
    print(vres.succeeded)
    print()
    print()
    print(vres.format())

        ###########################
    print("__________TEST__________")
    # print()
    # print(vres.parnames)
    # print(vres.parvals)
    print()

    var = n_piece+1
    for i in range(3*len(peaks)):
        print("{}       :      {}".format(vres.parnames[var+i],vres.parvals[var+i]))
    print()   
    print(pd.DataFrame(np.array(vres.parnames[var:])))
    print(pd.DataFrame(np.array(vres.parvals[var:])))
    print() 
    print("___________END__________")

        ############################
    fplot = FitPlot()
    mplot.prepare(d, model)
    fplot.prepare(dplot, mplot)
    fplot.plot()

    title = 'number of spline points: ' + str(n_piece+1)
    plt.title(title)
    plt.plot(wave_subset, flux_subset-model(wave_subset))
    plt.show()

    return None


if __name__ == "__main__":


    # =======================
    # Real data from HD170740
    # =======================


    # load data from file

    # # NaI 5890/5896
# <<<<<<< HEAD
#     file = '/home/ranjan/python/data/DR3_fits/HD170740/RED_564/HD170740_w564_n9_20160612_U.fits'
#     xmin = 5885.
#     xmax = 5898.
# # =======
    # file = '/data/DR3_fits/HD170740/RED_564/HD170740_w564_n9_20160612_U.fits'
    # xmin = 5885.
    # xmax = 5898.
    # scaling = 250.
# >>>>>>> 0bd5ecd769e2776ec33c0290a06cd233fa08264e

    # NaI 3000
    # file = '/home/ranjan/python/data/DR3_fits/HD170740/BLUE_346/HD170740_w346_n6_20160612_B.fits'
    # xmin = 3300.
    # xmax = 3305.
    # scaling = 30.

    # # KI 7665
    # file = '/data/DR3_fits/HD170740/RED_860/HD170740_w860_redl_20140915_O12.fits'
    # xmin = 7662
    # xmax = 7670
    # scaling = 1.

# <<<<<<< HEAD
    # KI 7665
    # file = '/data/DR3_fits/HD170740/RED_860/HD170740_w860_redl_20140915_O12.fits'

    # xmin = 7662
    # xmax = 7670
    # peak_cutoff = 0.5
    # n_points = 5
    # alpha = 0.05
    # gamma = 0.005
    # scaling = 75.




    # ____________________Using line merger___________________

    file1 = '/home/ranjan/python/data/DR3_fits/HD170740/BLUE_346/HD170740_w346_n6_20160612_B.fits'
    xmin1 = 3300.
    xmax1 = 3305.

    file2 = '/home/ranjan/python/data/DR3_fits/HD170740/RED_564/HD170740_w564_n9_20160612_U.fits'
    xmin2 = 5885.
    xmax2 = 5898.

    wave_subsetf, flux_subsetf = line_merger(file1,xmin1,xmax1,file2,xmin2,xmax2)

    plt.plot(wave_subsetf,flux_subsetf)
    plt.show()

    scaling = 30.

    # ________________________________________________________
# =======

    peak_cutoff = 0.5
    n_points = 5

    b_eff=6.47
    Gamma=6.064e9
# >>>>>>> 0bd5ecd769e2776ec33c0290a06cd233fa08264e


    fitter(wave_subsetf, flux_subsetf, peak_cutoff=peak_cutoff, n_points=n_points, 
                    b_eff=b_eff, Gamma=Gamma, scaling=scaling)







print('finished!')