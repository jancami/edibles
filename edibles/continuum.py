import matplotlib.pyplot as plt

from edibles.utils.edibles_spectrum import EdiblesSpectrum
from edibles.models import ContinuumModel

class Continuum():
    '''A class that has multiple methods for fitting different types of continua.

    Args:
        Spectrum (EdiblesSpectrum): The input EiblesSpectrum data
        method (str): The method of continuum fitting. spline (default), alphashape, polynomial


    '''


    def __init__(self, Spectrum, method='spline', *args, **kwargs):

        self.method = method
        self.Spectrum = Spectrum

        if method == 'spline':
            self.spline(*args, **kwargs)
        elif method == 'alphashape':
            self.alphashape(*args, **kwargs)
        elif method == 'polynomial':
            self.polynomial(*args, **kwargs)



    def spline(self, *args, **kwargs):
        print("method: ", self.method)

        n_anchors = kwargs['n_anchors']

        model = ContinuumModel(n_anchors=n_anchors)
        cont_pars = model.guess(self.Spectrum.flux, x=self.Spectrum.wave)


        result = model.fit(data=self.Spectrum.flux, params=cont_pars, x=self.Spectrum.wave)

        result.plot_fit()
        plt.show()
        return model


    def alphashape(self):
        print("method: ", self.method)

    def polynomial(self):
        print("method: ", self.method)







if __name__ == '__main__':

    # sp = EdiblesSpectrum("/HD170740/RED_860/HD170740_w860_redl_20140915_O12.fits")
    sp = EdiblesSpectrum("/HD23466/BLUE_346/HD23466_w346_blue_20180731_O11.fits")

    subset = sp.getSpectrum(xmin=3270, xmax=3305)

    # Continuum(x, y, method='spline', anchors=4)
    model = Continuum(sp, method='spline', n_anchors=4)


