import matplotlib.pyplot as plt
import time

from sherpa.data import Data1D, DataSimulFit
from sherpa.stats import LeastSq
from sherpa.optmethods import LevMar
from sherpa.fit import Fit, SimulFitModel
from sherpa.plot import DataPlot, ModelPlot, FitPlot, SplitPlot


def fit(star_name, data, model, silent=False, breakdown=False):
    """A function that will fit a given multi-part model to a given spectrum.



    :param star_name: Name of the target star
    :type star_name: string
    :param data: Spectrum data in the form (wave, flux)
    :type data: tuple
    :param model: An unfit spectrum model
    :type model: object
    :param silent:  If true, no plots will generate, defaults to False
    :type silent: bool

    :return: model that is fit to the data
    :rtype: object


    """

    wave, flux = data

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    d = Data1D(star_name, wave, flux)

    # ==========================================
    # Initial guesses

    # Dataset 1
    dplot = DataPlot()
    dplot.prepare(d)
    if silent is False:
        dplot.plot()

    mplot = ModelPlot()
    mplot.prepare(d, model)
    if silent is False:
        dplot.plot()
        mplot.overplot()
        plt.show()

    # =========================================
    # Fitting happens here - don't break please
    start = time.time()

    stat = LeastSq()

    opt = LevMar()

    opt.verbose = 0
    opt.ftol = 1e-15
    opt.xtol = 1e-15
    opt.gtol = 1e-15
    opt.epsfcn = 1e-15

    if silent == False:
        print(opt)

    vfit = Fit(d, model, stat=stat, method=opt)

    if silent == False:
        print(vfit)

    vres = vfit.fit()

    if silent == False:
        print()
        print()
        print(vres.format())

    # =========================================
    # Plotting after fit

    # Dataset 1
    if silent is False:
        fplot = FitPlot()
        mplot.prepare(d, model)
        fplot.prepare(dplot, mplot)
        fplot.plot()

        # residual
        plt.title(star_name)
        plt.plot(wave, flux - model(wave))

        # plt.xaxis(fontsize = )
        plt.xlabel("Wavelength (AA)", fontsize=12)
        plt.ylabel("Flux", fontsize=12)
        plt.tick_params(axis="both", labelsize=12)

    if silent is False:
        duration = time.time() - start
        print()
        print("Time taken: " + str(duration))
        print()

    plt.show()

    if breakdown is True:
        params = []

        cont = model[0]

        if silent is False:
            plt.scatter(wave, flux, marker=".", c="black")
            plt.plot(wave, model(wave), c="C1")

        for line in model:
            if line.name[0] != "(":
                if line.name == "Cont_flux":
                    if silent is False:
                        print(line)
                        plt.plot(wave, line(wave), linestyle="--")
                else:
                    params.append(line)
                    if silent is False:
                        print()
                        print(line)
                        plt.plot(wave, line(wave) * cont(wave), linestyle="--")

        plt.show()

        return model, params

    return model


def multifit(star_name, data_list, model_list, silent=False):
    """A function that will fit 2 models to 2 spectra simultaneously.
        This was created to fit the NaI doublets at ~3300 and ~5890 Angstroms.

    :param star_name: Name of the target star
    :type star_name: string
    :param data_list: List of spectrum data in the form [(wave, flux), (wave, flux),...]
    :type data_list: tuple
    :param model_list:  A list of unfit spectrum models
    :type model_list: model instance
    :param silent:  If true, no plots will generate, defaults to False
    :type silent: bool

    :return: models that are fit to the data
    :rtype: list

    """


    wave1, flux1 = data_list[0]
    wave2, flux2 = data_list[1]

    model1 = model_list[0]
    model2 = model_list[1]

    name_1 = star_name + ' 1'
    name_2 = star_name + ' 2'

    d1 = Data1D(name_1, wave1, flux1)
    d2 = Data1D(name_2, wave2, flux2)

    dall = DataSimulFit("combined", (d1, d2))
    mall = SimulFitModel("combined", (model1, model2))

    # # ==========================================
    # # Initial guesses

    # Dataset 1
    dplot1 = DataPlot()
    dplot1.prepare(d1)
    if silent is False:
        dplot1.plot()

    mplot1 = ModelPlot()
    mplot1.prepare(d1, model1)
    if silent is False:
        dplot1.plot()
        mplot1.overplot()
        plt.show()

        # Dataset 2
    dplot2 = DataPlot()
    dplot2.prepare(d2)
    if silent is False:
        dplot2.plot()

    mplot2 = ModelPlot()
    mplot2.prepare(d2, model2)
    if silent is False:
        dplot2.plot()
        mplot2.overplot()
        plt.show()

    # # =========================================
    # # Fitting happens here - don't break please
    stat = LeastSq()

    opt = LevMar()
    opt.verbose = 0
    opt.ftol = 1e-15
    opt.xtol = 1e-15
    opt.gtol = 1e-15
    opt.epsfcn = 1e-15
    print(opt)

    vfit = Fit(dall, mall, stat=stat, method=opt)
    print(vfit)
    vres = vfit.fit()

    print()
    print()
    print("Did the fit succeed? [bool]")
    print(vres.succeeded)
    print()
    print()
    print(vres.format())

    # # =========================================
    # # Plotting after fit
    if silent is False:
        # Dataset 1
        fplot1 = FitPlot()
        mplot1.prepare(d1, model1)
        fplot1.prepare(dplot1, mplot1)
        fplot1.plot()

        # residual
        title = "Data 1"
        plt.title(title)
        plt.plot(wave1, flux1 - model1(wave1))
        plt.show()

        # Dataset 2
        fplot2 = FitPlot()
        mplot2.prepare(d2, model2)
        fplot2.prepare(dplot2, mplot2)
        fplot2.plot()

        # residual
        title = "Data 2"
        plt.title(title)
        plt.plot(wave2, flux2 - model2(wave2))
        plt.show()

        # both datasets - no residuals
        splot = SplitPlot()
        splot.addplot(fplot1)
        splot.addplot(fplot2)

        plt.tight_layout()
        plt.show()

    return model_list
