import matplotlib.pyplot as plt
from lmfit import Parameters

from edibles.models import ContinuumModel, VoigtModel
from edibles.utils.edibles_spectrum import EdiblesSpectrum


class Sightline():
    '''A model of the sightline between the telescope and the target star.


    Args:
        spectrum (EdiblesSpectrum): The input spectrum object


    '''

    def __init__(self, spectrum, init_cont=True, n_anchors=4):

        self.__dict__.update(spectrum.__dict__)

        self.wave = spectrum.wave
        self.flux = spectrum.flux


        if init_cont:
            cont_model = ContinuumModel(n_anchors=n_anchors)
            cont_pars = cont_model.guess(self.flux, x=self.wave)


        self.model = cont_model
        self.model_pars = cont_pars

        self.num_sources = 0
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
        par.add(name + '_b', value=similar['b'], min=0)

        self.model_pars = self.model_pars + par

        # print(list(par.valuesdict().keys())[0])


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

        par_name = source + '_b'
        new_pars[source + '_' + name + '_b'].set(expr=par_name)

        self.model = self.model * new_line
        self.model_pars = self.model_pars + new_pars



    def fit(self, report=False, plot=False):
        '''Fits the sightline models to the sightline data given by the EdiblesSpectrum object.

        Args:
            report (bool): default False: If true, prints the report from the fit.
            plot (bool): default False: If true, plots the data and the fit model.

        '''

        self.result = self.model.fit(data=self.flux, params=self.model_pars, x=self.wave)

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

    # Add line with undefined source
    d = {'d': 0.01, 'tau_0': 0.1, 'lam_0': 7665.2}

    sightline.add_line(name='line3', source='interstellar')

    # Add line with no source & user defined pars
    d = {'d': 0.01, 'tau_0': 0.1, 'lam_0': 7662}
    sightline.add_line(name='line4', pars=d)

    # ###############################################################
    # Fit and plot
    sightline.fit(report=True, plot=True)

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
    sightline.fit(report=True, plot=True)

    out = sightline.model.eval(data=sp1.flux, params=sightline.result.params, x=sp1.wave)
    resid = sp1.flux - out

    plt.plot(sp1.wave, sp1.flux)
    plt.plot(sp1.wave, out)
    plt.plot(sp1.wave, resid)
    plt.show()
    # ###############################################################
