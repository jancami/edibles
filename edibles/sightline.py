import numpy as np
import matplotlib.pyplot as plt
import bisect
from lmfit import Parameters
import astropy.constants as cst

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
        self.num_prior_lines = 0
        self.source_names = []

        self.add_source("Telluric", similar={'b': 2})
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

        self.old_complete_model = self.complete_model
        self.complete_model = self.complete_model * new_line

        self.old_all_pars = self.all_pars
        self.all_pars = self.all_pars + new_pars

        self.old_cont_model = self.cont_model
        self.old_cont_pars = self.cont_model_pars

        if source == "Telluric":
            try:
                self.old_telluric_model = self.telluric_model
                self.telluric_model = self.telluric_model * new_line
            except AttributeError:
                self.old_telluric_model = new_line
                self.telluric_model = new_line

            try:
                self.old_telluric_pars = self.telluric_pars
                self.telluric_pars = self.telluric_pars + new_pars
            except AttributeError:
                print('Something bad is probably happening')
                self.old_telluric_pars = new_pars
                self.telluric_pars = new_pars

        else:
            try:
                self.old_nontelluric_model = self.nontelluric_model
                self.nontelluric_model = self.nontelluric_model * new_line
            except AttributeError:
                self.old_nontelluric_model = new_line
                self.nontelluric_model = new_line
            try:
                self.old_nontelluric_pars = self.nontelluric_pars
                self.nontelluric_pars = self.nontelluric_pars + new_pars
            except AttributeError:
                self.old_nontelluric_pars = new_pars
                self.nontelluric_pars = new_pars


        lambda_name = source + '_' + name + '_lam_0'
        index = bisect.bisect(self.peaks, new_pars[lambda_name])
        self.peaks.insert(index, new_pars[lambda_name])

        self.most_recent = source + '_' + name
        self.n_lines += 1


    def fit(self, data=None, old=False, x=None, report=False,
            plot=False, weights=None, method='leastsq', **kwargs):
        '''Fits a model to the sightline data given by the EdiblesSpectrum object.

        Args:
            data (1darray): Flux data to fit
            params (lmfit.parameter.Parameters): Initial parameters to fit
            model (lmfit.model.CompositeModel): The model to fit, default: self.complete_model
            x (1darray): Wavelength data to fit
            report (bool): default False: If true, prints the report from the fit.
            plot (bool): default False: If true, plots the data and the fit model.
            method (str): The method of fitting. default: leastsq

        '''
        if data is None:
            data = self.flux
        if x is None:
            x = self.wave

        if old is True:
            model = self.old_complete_model
            params = self.old_all_pars
        else:
            model = self.complete_model
            params = self.all_pars

        self.result = model.fit(data=data,
                                params=params,
                                x=x,
                                weights=weights,
                                method=method,
                                **kwargs)
        if report:
            print(self.result.fit_report())
            self.result.params.pretty_print()
        if plot:
            self.result.plot_fit()
            plt.show()


        # Update parameter values after fit - for use in model separation
        self.all_pars = self.result.params

        # create new parameters object and add to it from the results parameters
        if old is False:
            try:
                tell_pars = Parameters()
                for par_name in self.telluric_pars:
                    tell_pars.add(self.all_pars[par_name])

                # update attribute
                assert len(self.telluric_pars) == len(tell_pars)
                self.telluric_pars = tell_pars
            except AttributeError:
                pass

            try:
                non_tell_pars = Parameters()
                for par_name in self.nontelluric_pars:
                    non_tell_pars.add(self.all_pars[par_name])
                assert len(self.nontelluric_pars) == len(non_tell_pars)
                self.nontelluric_pars = non_tell_pars
            except AttributeError:
                pass

            try:
                cont_pars = Parameters()
                for par_name in self.cont_model_pars:
                    cont_pars.add(self.all_pars[par_name])
                assert len(self.cont_model_pars) == len(cont_pars)
                self.cont_model_pars = cont_pars
            except AttributeError:
                pass


    def freeze(self, pars=None, prefix=None, freeze_cont=True, unfreeze=False):
        '''Freezes the current params, so you can still add to the
            model but the 'old' parameters will not change

        Args:
            prefix (str): Prefix of parameters to freeze, default: None, example: 'Telluric'
            freeze_cont (bool): Freeze the continuum or not, default: True
            unfreeze (bool): unfreezes all parameters except x values of
                spline anchors, default=False

        '''
        if pars is None:
            pars = self.all_pars

        if unfreeze is False:
            if prefix:
                for par in pars:
                    if prefix in par:
                        pars[par].set(vary=False)

            else:
                for par in pars:
                    pars[par].set(vary=False)

            if not freeze_cont:
                for par in pars:
                    if 'y_' in par:
                        pars[par].set(vary=True)

        if unfreeze is True:
            for par in pars:

                if ('y_' in par):
                    pars[par].set(vary=True)

                if ('Telluric' in par) and (par[-2:] != '_b'):
                    pars[par].set(vary=True)
                pars['Telluric_b'].set(vary=True)

                if ('Nontelluric' in par) and (par[-2:] != '_d'):
                    pars[par].set(vary=True)



    def separate(self, data, x, old=False, plot=True):
        '''Separate the sources that were added to Sightline.

        Args:
            data (1darray): FLux data to use for separation
            x (1darray): Wavelength array to use
            old (bool): If true, uses the older, second-most recent model and parameters
            plot (bool): If true, plots separted spectrum

        '''

        assert len(self.telluric_pars) > 0
        assert len(self.nontelluric_pars) > 0

        if old is True:
            model = self.old_complete_model
            params = self.old_all_pars
            telluric_model = self.old_telluric_model
            telluric_params = self.old_telluric_pars
            nontelluric_model = self.old_nontelluric_model
            nontelluric_params = self.old_nontelluric_pars
            cont_model = self.old_cont_model
            cont_params = self.old_cont_pars

        else:
            model = self.complete_model
            params = self.all_pars
            telluric_model = self.telluric_model
            telluric_params = self.telluric_pars
            nontelluric_model = self.nontelluric_model
            nontelluric_params = self.nontelluric_pars
            cont_model = self.cont_model
            cont_params = self.cont_model_pars

        if len(self.source_names) == 2:
            complete_out = model.eval(
                data=data,
                params=params,
                x=x
            )
            telluric_out = telluric_model.eval(
                data=data,
                params=telluric_params,
                x=x
            )
            nontelluric_out = nontelluric_model.eval(
                data=data,
                params=nontelluric_params,
                x=x
            )
            cont_out = cont_model.eval(
                data=data,
                params=cont_params,
                x=x
            )

            if plot:

                plt.plot(x, data, label='Data', color='k')
                plt.plot(x, complete_out, label='Final model', color='r')
                plt.plot(x, data - complete_out, label='Residual', color='g')
                plt.plot(x, telluric_out * cont_out, label='Telluric model')
                plt.plot(x, nontelluric_out * cont_out, label='Non-telluric model')
                plt.xlabel(r'Wavelength ($\AA$)', fontsize=14)
                plt.ylabel('Flux', fontsize=14)
                plt.legend()

                plt.show()

            return complete_out, telluric_out, nontelluric_out, cont_out


if __name__ == "__main__":


    FILE1 = "/HD170740/RED_860/HD170740_w860_redl_20140915_O12.fits"
    xmin = 7661.75
    xmax = 7669

    sp1 = EdiblesSpectrum(FILE1)
    sp1.getSpectrum(xmin=xmin, xmax=xmax)

    sightline = Sightline(sp1, n_anchors=5)


    # Add line with auto-guessed params
    sightline.add_line(name='line1', source='Telluric')

    # Add line with user defined params
    pars = {'d': 0.01, 'tau_0': 0.6, 'lam_0': 7664.8}
    sightline.add_line(name='line2', pars=pars, source='Telluric')


    # # ###############################################################
    # # Fit and plot
    sightline.fit(report=True, plot=False, method='leastsq')

    out = sightline.complete_model.eval(data=sp1.flux, params=sightline.result.params, x=sp1.wave)
    resid = sp1.flux - out



    # Add line with different source

    lam_0 = 7665.25

    K_Gamma = 3.820e7
    K_d = K_Gamma * lam_0**2 / (4 * np.pi * (cst.c.to("cm/s").value * 1e8))


    pars = {'d': K_d, 'tau_0': 0.07, 'lam_0': lam_0}
    sightline.add_line(name='line3', source='Nontelluric', pars=pars)
    sightline.all_pars['Nontelluric_line3_d'].set(vary=False)

    # sightline.fit(report=True, plot=False, method='leastsq')
    # out = sightline.complete_model.eval(data=sp1.flux, params=sightline.result.params, x=sp1.wave)
    # resid = sp1.flux - out


    lam_0 = 7665.33
    pars = {'d': K_d, 'tau_0': 0.01, 'b': 1, 'lam_0': lam_0}
    sightline.add_line(name='line4', source='Nontelluric', pars=pars)
    sightline.all_pars['Nontelluric_line4_d'].set(vary=False)
    # sightline.fit(report=True, plot=False, method='leastsq')


    lam_0 = 7665.15
    pars = {'d': K_d, 'tau_0': 0.001, 'b': 1, 'lam_0': lam_0}
    sightline.add_line(name='line5', source='Nontelluric', pars=pars)
    sightline.all_pars['Nontelluric_line5_d'].set(vary=False)
    sightline.fit(report=True, plot=False, method='leastsq')





    pars = {'d': 0.01, 'tau_0': 0.01, 'b': 1, 'lam_0': 7662}
    sightline.add_line(name='line6', source='Telluric', pars=pars)
    sightline.fit(report=True, plot=False, method='leastsq')



    pars = {'d': 0.01, 'tau_0': 0.01, 'b': 1, 'lam_0': 7663.7}
    sightline.add_line(name='line7', source='Telluric', pars=pars)
    sightline.fit(report=True, plot=False, method='leastsq')


    pars = {'d': 0.01, 'tau_0': 0.01, 'b': 1, 'lam_0': 7666.5}
    sightline.add_line(name='line8', source='Telluric', pars=pars)
    sightline.fit(report=True, plot=False, method='leastsq')


    pars = {'d': 0.01, 'tau_0': 0.01, 'b': 1, 'lam_0': 7667.5}
    sightline.add_line(name='line9', source='Telluric', pars=pars)
    sightline.fit(report=True, plot=False, method='leastsq')






    out = sightline.complete_model.eval(data=sp1.interp_flux, params=sightline.result.params,
                                        x=sp1.grid)
    resid = sp1.interp_flux - out


    plt.plot(sp1.grid, sp1.interp_flux)
    plt.plot(sp1.grid, out)
    plt.plot(sp1.grid, resid)
    plt.show()

    sightline.separate(data=sp1.interp_flux, x=sp1.grid)
