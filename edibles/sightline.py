import numpy as np
import matplotlib.pyplot as plt
import bisect
from lmfit import Parameters

from edibles.models import ContinuumModel, VoigtModel
from edibles.utils.edibles_spectrum import EdiblesSpectrum


class Sightline:
    '''A model of the sightline between the telescope and the target star.

    Args:
        Spectrum (EdiblesSpectrum): The input spectrum object
        n_anchors (int): Optional, The number of anchors in the ContinuumSpline
    '''

    def __init__(self, Spectrum, init_cont=True, n_anchors=4):

        self.__dict__.update(Spectrum.__dict__)

        self.wave = Spectrum.wave
        self.flux = Spectrum.flux
        self.Spectrum = Spectrum

        if init_cont:
            cont_model = ContinuumModel(n_anchors=n_anchors)
            cont_pars = cont_model.guess(self.flux, x=self.wave)

            for yname in cont_model.ynames:
                flux_range = np.max(self.flux) - np.min(self.flux)
                ymin = cont_pars[yname].value - (flux_range / 2)
                ymax = cont_pars[yname].value + (flux_range / 2)

                cont_pars[yname].set(min=ymin, max=ymax)

        self.model = cont_model
        self.model_pars = cont_pars

        self.peaks = []

        self.num_sources = 0
        self.n_anchors = n_anchors
        self.source_names = []


    def add_source(self, name, similar={'b': 3}):
        '''Adds a new source of absorption to the sightline.

        The purpose of a source is to house parameters similar across multiple lines- ex: b

        Args:
            name (str): The name of the absorption source
            similar (dict): A dict of parameters that change with the source,
                not the specific line, default: {'b': 3}

        '''

        self.num_sources += 1
        self.source_names.append(name)

        par = Parameters()
        par.add(name + '_b', value=similar['b'], min=0, max=30)

        self.model_pars = self.model_pars + par


    def add_line(self, name, source=None, pars=None, guess_data=None):
        '''Adds a new line to a given absorption source.
        If no source is given, a new one will be created.

        Args:
            name (str): The name of the line
            source (str): the name of the source this line will belong to
            pars (dict): user input parameters
            guess_data (1darray): flux data to guess with

        '''

        if source is None:
            self.add_source(name)
            source = name

        if source is not None:
            if source not in self.source_names:
                print()
                print('Could not find source \'{}\' in source_names.'.format(source))
                print('Creating source \'{}\''.format(source))
                self.add_source(source)

        new_line = VoigtModel(prefix=source + '_' + name + '_')

        if guess_data is not None:
            new_pars = new_line.guess(guess_data, x=self.wave)
        else:
            new_pars = new_line.guess(self.flux, x=self.wave)

        if pars is not None:
            for par in pars:  # lam_0...
                par_name = source + '_' + name + '_' + par  # telluric_line1_lam_0...
                new_pars[par_name].set(value=pars[par])

        b_name = source + '_b'
        new_pars[source + '_' + name + '_b'].set(expr=b_name)

        new_pars[source + '_' + name + '_lam_0'].set(
            min=self.Spectrum.xmin, max=self.Spectrum.xmax
        )

        self.model = self.model * new_line
        self.model_pars = self.model_pars + new_pars

        lambda_name = source + '_' + name + '_lam_0'
        index = bisect.bisect(self.peaks, new_pars[lambda_name])
        self.peaks.insert(index, new_pars[lambda_name])


        if len(self.peaks) > 1:
            for idx in range(len(self.peaks)):
                if idx == 0:
                    self.model_pars[self.peaks[idx].name].set(max=self.peaks[idx + 1].value)
                elif idx == len(self.peaks) - 1:
                    self.model_pars[self.peaks[idx].name].set(min=self.peaks[idx - 1].value)
                else:
                    self.model_pars[self.peaks[idx].name].set(
                        min=self.peaks[idx - 1].value,
                        max=self.peaks[idx + 1].value
                    )


    def fit(self, data=None, params=None,
            x=None, report=False, plot=False, method='leastsq', bary=False):
        '''Fits the sightline models to the sightline data given by the EdiblesSpectrum object.

        Args:
            data (1darray): Flux data to fit
            params (lmfit.parameter.Parameters): Initial parameters to fit
            x (1darray): Wavelength data to fit
            report (bool): default False: If true, prints the report from the fit.
            plot (bool): default False: If true, plots the data and the fit model.
            method (str): The method of fitting. default: leastsq
            bary (bool): If true, creates bary_result instead of result

        '''
        if data is None:
            data = self.flux
        if params is None:
            params = self.model_pars
        if x is None:
            x = self.wave

        if bary:
            self.bary_result = self.model.fit(data=data,
                                              params=params,
                                              x=x,
                                              method=method)
            if report:
                print(self.bary_result.fit_report())
                self.bary_result.params.pretty_print()
            if plot:
                self.bary_result.plot_fit()
                plt.show()

        else:
            self.result = self.model.fit(data=data,
                                         params=params,
                                         x=x,
                                         method=method)
            if report:
                print(self.result.fit_report())
                self.result.params.pretty_print()
            if plot:
                self.result.plot_fit()
                plt.show()


if __name__ == "__main__":


    FILE1 = "/HD170740/RED_860/HD170740_w860_redl_20140915_O12.fits"
    xmin = 7661.5
    xmax = 7669

    sp1 = EdiblesSpectrum(FILE1)
    sp1.getSpectrum(xmin=7661, xmax=7670)

    sightline = Sightline(sp1)

    # Add source
    sightline.add_source('telluric')

    # Add line with auto-guessed params
    sightline.add_line(name='line1', source='telluric')

    # Add line with user defined params
    d = {'d': 0.01, 'tau_0': 0.6, 'lam_0': 7664.8}
    sightline.add_line(name='line2', pars=d, source='telluric')

    # Add line with different source
    d = {'d': 0.01, 'tau_0': 0.1, 'lam_0': 7665.2}
    sightline.add_source('interstellar', similar={'b': 1.9})
    sightline.add_line(name='line3', source='interstellar', pars=d)

    # Add line with no source & user defined pars
    # d = {'d': 0.01, 'tau_0': 0.1, 'lam_0': 7662}
    # sightline.add_line(name='line4', pars=d)

    # ###############################################################
    # Fit and plot
    sightline.fit(report=True, plot=True, method='leastsq')

    out = sightline.model.eval(data=sp1.flux, params=sightline.result.params, x=sp1.wave)
    resid = sp1.flux - out

    plt.plot(sp1.wave, sp1.flux)
    plt.plot(sp1.wave, out)
    plt.plot(sp1.wave, resid)
    plt.show()
    # ###############################################################

    # Add line using guess_pars, and link parameters together
    sightline.add_line(name='line5', source='interstellar', guess_data=resid)
    sightline.model_pars['interstellar_line5_lam_0'].set(expr='interstellar_line3_lam_0 + 0.091')

    # ###############################################################
    # Fit and plot
    sightline.fit(report=True, plot=True, method='leastsq')

    out = sightline.model.eval(data=sp1.flux, params=sightline.result.params, x=sp1.wave)
    resid = sp1.flux - out

    plt.plot(sp1.wave, sp1.flux)
    plt.plot(sp1.wave, out)
    plt.plot(sp1.wave, resid)
    plt.show()
    # ###############################################################
