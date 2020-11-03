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

        self.cont_model = cont_model
        self.cont_model_pars = cont_pars

        self.complete_model = cont_model
        self.all_pars = cont_pars

        self.peaks = []

        self.n_anchors = n_anchors
        self.n_lines = 0
        self.source_names = []

        self.add_source("Telluric", similar={'b': 3})
        self.add_source("Nontelluric", similar=None)


    def add_source(self, name, similar=None):
        '''Adds a new source of absorption to the sightline.

        The purpose of a source is to hold multiple line models
            together, sometiimes with similar parameters

        Args:
            name (str): The name of the absorption source
            similar (dict): A dict of parameters that change with the source,
                not the specific line, default: None, example: similar={'b': 3}

        '''

        self.source_names.append(name)



        if name == "Telluric" and similar is not None:

            par = Parameters()
            for key in similar:
                par.add(name + '_' + key, value=similar[key], min=0, max=30)

            self.telluric_pars = par
            self.all_pars = self.all_pars + par


    def add_line(self, name, source=None, pars=None, guess_data=None):
        '''Adds a new line to a given absorption source.
        If no source is given, a new one will be created.

        Args:
            name (str): The name of the line
            source (str): the name of the source this line will belong to
            pars (dict): user input parameters
            guess_data (1darray): flux data to guess with

        '''

        assert source is not None, "Source must not be None"

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

        if source == "Telluric":
            b_name = source + '_b'
            new_pars[source + '_' + name + '_b'].set(expr=b_name)

        new_pars[source + '_' + name + '_lam_0'].set(
            min=self.Spectrum.xmin, max=self.Spectrum.xmax
        )


        self.complete_model = self.complete_model * new_line
        self.all_pars = self.all_pars + new_pars

        if source == "Telluric":
            try:
                self.telluric_model = self.telluric_model * new_line
            except AttributeError:
                self.telluric_model = new_line

            try:
                self.telluric_pars = self.telluric_pars + new_pars
            except AttributeError:
                print('Something bad is probably happening')
                self.telluric_pars = new_pars

        else:
            try:
                self.nontelluric_model = self.nontelluric_model * new_line
            except AttributeError:
                self.nontelluric_model = new_line
            try:
                self.nontelluric_pars = self.nontelluric_pars + new_pars
            except AttributeError:
                self.nontelluric_pars = new_pars


        lambda_name = source + '_' + name + '_lam_0'
        index = bisect.bisect(self.peaks, new_pars[lambda_name])
        self.peaks.insert(index, new_pars[lambda_name])

        self.most_recent = source + '_' + name
        self.n_lines += 1


    def fit(self, data=None, params=None,
            x=None, report=False, plot=False, weights=None, method='leastsq', bary=False):
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
            params = self.all_pars
        if x is None:
            x = self.wave


        self.result = self.complete_model.fit(data=data,
                                              params=params,
                                              x=x,
                                              weights=weights,
                                              method=method)
        if report:
            print(self.result.fit_report())
            self.result.params.pretty_print()
        if plot:
            self.result.plot_fit()
            plt.show()


        # ##################################################
        # Update parameter values after fit - for use in model separation

        self.all_pars = self.result.params

        # create new parameters object and add to it from the results parameters
        try:
            tell_pars = Parameters()
            for par_name in self.telluric_pars:
                tell_pars.add(self.result.params[par_name])

            # update attribute
            assert len(self.telluric_pars) == len(tell_pars)
            self.telluric_pars = tell_pars
        except AttributeError:
            pass

        try:
            non_tell_pars = Parameters()
            for par_name in self.nontelluric_pars:
                non_tell_pars.add(self.result.params[par_name])
            assert len(self.nontelluric_pars) == len(non_tell_pars)
            self.nontelluric_pars = non_tell_pars
        except AttributeError:
            pass

        try:
            cont_pars = Parameters()
            for par_name in self.cont_model_pars:
                cont_pars.add(self.result.params[par_name])
            assert len(self.cont_model_pars) == len(cont_pars)
            self.cont_model_pars = cont_pars
        except AttributeError:
            pass


    def freeze(self, prefix=None, freeze_cont=True):
        '''Freezes the current params, so you can still add to the
model but the 'old' parameters will not change

        Args:
            prefix (str): Prefix of parameters to freeze, default: None, example: 'Telluric'
            freeze_cont (bool): Freeze the continuum or not, default: True

        '''

        if prefix:
            for par in self.all_pars:
                if prefix in par:
                    self.all_pars[par].set(vary=False)
                if 'y_' in par:
                    self.all_pars[par].set(vary=False)

        else:
            for par in self.all_pars:
                self.all_pars[par].set(vary=False)

        if not freeze_cont:
            for par in self.all_pars:
                if 'y_' in par:
                    self.all_pars[par].set(vary=True)


    def separate(self, data, x):
        '''Separate the sources

        '''

        assert len(self.telluric_pars) > 0
        assert len(self.nontelluric_pars) > 0

        if len(self.source_names) == 2:
            complete_out = self.complete_model.eval(
                data=data,
                params=self.result.params,
                x=x
            )
            telluric_out = self.telluric_model.eval(
                data=data,
                params=self.telluric_pars,
                x=x
            )
            nontelluric_out = self.nontelluric_model.eval(
                data=data,
                params=self.nontelluric_pars,
                x=x
            )
            cont_out = self.cont_model.eval(
                data=data,
                params=self.cont_model_pars,
                x=x
            )

            fig, axs = plt.subplots(1, 2)

            axs[0].plot(x, data, label='Data')
            axs[0].plot(x, complete_out, label='Complete model')
            axs[0].plot(x, cont_out, label='Continuum model')
            axs[0].plot(x, (data - complete_out), label='Residual')
            axs[0].legend()


            norm_data = data / cont_out

            axs[1].plot(x, norm_data, label='Continuum-Normalized data')
            axs[1].plot(x, telluric_out, label='Telluric model')
            axs[1].plot(x, nontelluric_out, label='Non-Telluric model')
            axs[1].plot(x, (norm_data - telluric_out), label='Non-telluric data')
            axs[1].plot(x, (norm_data - nontelluric_out), label='Telluric data')
            axs[1].legend()

            plt.show()


if __name__ == "__main__":


    FILE1 = "/HD170740/RED_860/HD170740_w860_redl_20140915_O12.fits"
    xmin = 7661.5
    xmax = 7669

    sp1 = EdiblesSpectrum(FILE1)
    sp1.getSpectrum(xmin=7661, xmax=7670)

    sightline = Sightline(sp1)


    # Add line with auto-guessed params
    sightline.add_line(name='line1', source='Telluric')

    # Add line with user defined params
    d = {'d': 0.01, 'tau_0': 0.6, 'lam_0': 7664.8}
    sightline.add_line(name='line2', pars=d, source='Telluric')

    # Add line with different source
    d = {'d': 0.01, 'tau_0': 0.1, 'lam_0': 7665.1}
    sightline.add_line(name='line3', source='Nontelluric', pars=d)

    # # ###############################################################
    # # Fit and plot
    sightline.fit(report=True, plot=True, method='leastsq')

    out = sightline.complete_model.eval(data=sp1.flux, params=sightline.result.params, x=sp1.wave)
    resid = sp1.flux - out

    plt.plot(sp1.wave, sp1.flux)
    plt.plot(sp1.wave, out)
    plt.plot(sp1.wave, resid)
    plt.show()

    # sightline.separate(data=sp1.interp_flux, x=sp1.grid)

    # sightline.freeze(prefix='Telluric', freeze_cont=False)

    # Add line using guess_pars, and link parameters together
    sightline.add_line(name='line4', source='Nontelluric', guess_data=resid)

    sightline.fit(report=True, plot=True, method='leastsq')

    sightline.separate(data=sp1.interp_flux, x=sp1.grid)

