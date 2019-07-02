import numpy as np
import matplotlib.pyplot as plt

from sherpa.data import Data1D
from sherpa.stats import LeastSq
from sherpa.optmethods import LevMar, NelderMead
from sherpa.fit import Fit
from sherpa.plot import DataPlot, ModelPlot, FitPlot

import time 



def fit(star_name, data, model):
    '''A function that will fit a given multi-part model to a given spectrum.


    INPUT:

    spectrum:   [tuple]             In the form (wave, flux)
    model:      [model instance]    Should contain CONTINUUM and ALL LINES - an initial model of the spectrum


    OUTPUT:

    fit_model:  [model instance]    Same model as input but with parameters tuned to match spectrum

    '''



    wave, flux = data




    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # setup to fit / plot

    d = Data1D('2014-09-15', wave, flux)


    # ==========================================
    # Initial guesses

        # Dataset 1
    dplot = DataPlot()
    dplot.prepare(d)
    dplot.plot()

    mplot = ModelPlot()
    mplot.prepare(d, model)
    dplot.plot()
    mplot.overplot()
    plt.show()


    # =========================================
    # Fitting happens here - don't break please
    start = time.time()

    stat = LeastSq()

    opt = LevMar()

    opt.verbose = 3
    opt.ftol = 1e-15
    opt.xtol = 1e-15
    opt.gtol = 1e-15
    opt.epsfcn = 1e-15

    print(opt)

    vfit = Fit(d, model, stat=stat, method=opt)


    print(vfit)
    vres = vfit.fit()

    print()
    print()
    print(vres.format())

    # =========================================
    # Plotting after fit

        # Dataset 1
    fplot = FitPlot()
    mplot.prepare(d, model)
    fplot.prepare(dplot, mplot)
    fplot.plot()

        # residual
    title = '2014-09-15'
    plt.title(title)
    plt.plot(wave, flux-model(wave))

    # plt.xaxis(fontsize = )
    plt.xlabel('Wavelength (AA)', fontsize=12)
    plt.ylabel('Flux', fontsize=12)
    plt.tick_params(axis='both', labelsize=12)

    print()
    print()
    duration = time.time() - start
    print('Time taken: ' + str(duration))


    # print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    # print('RESULTS FOR ' + star_name)
    # print('Line #    cent           b             d           tau_0')
    # print('1         {:.5f}     {:.5f}       {:.5f}     {:.5f}'.format(line1.lam_0.val, line1.b.val, line1.d.val, line1.tau_0.val))
    # print('2         {:.5f}     {:.5f}       {:.5f}     {:.5f}'.format(line2.lam_0.val, line2.b.val, line2.d.val, line2.tau_0.val))
    # print('3         {:.5f}     {:.5f}       {:.5f}     {:.5f}'.format(line3.lam_0.val, line3.b.val, line3.d.val, line3.tau_0.val))
    # # print('4         {:.5f}     {:.5f}       {:.5f}     {:.5f}'.format(obj4.lam_0.val, obj4.b.val, obj4.d.val, obj4.tau_0.val))
    # print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    # print()
    # print(type(line1.lam_0.val))

    plt.show()


    return model
