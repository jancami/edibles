import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np

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

        print(result.params)

        result.plot_fit()
        plt.show()
        return model


    def alphashape(self):
        print("method: ", self.method)

    def polynomial(self):
        print("method: ", self.method)



def copy_to_FITS(filename, x_points, y_points, name, initial=False):


    hdulist = fits.open(filename)
    spectrum = hdulist['PRIMARY']

    # INITIAL
    if len(hdulist) == 1:

        print('Adding INITIAL continuum points.')



        columns = fits.ColDefs(
            [fits.Column(name='initial_x', format='D10.4', array=x_points),
             fits.Column(name='initial_y', format='D10.4', array=y_points)]
        )
        continuum_hdu = fits.BinTableHDU.from_columns(columns)


        hdulist.append(continuum_hdu)




        hdulist[1].name = 'CONTINUUM'
        hdulist['CONTINUUM'].header['METHOD'] = 'spline'



        # TODO: update filename to save to

        hdulist.writeto('test.fits', overwrite=True)

        hdulist.info()


    # ADD LOCAL SOLINE POINTS
    elif len(hdulist) > 1:

        print('Adding NEW continuum points.')

        continuum_hdu = hdulist['CONTINUUM']
        data = continuum_hdu.data


        new_column = fits.ColDefs(
            [fits.Column(name=name + '_x', format='D', array=x_points),
             fits.Column(name=name + '_y', format='D', array=y_points)]
        )

        data = fits.BinTableHDU.from_columns(data.columns + new_column)


        # hdulist['CONTINUUM'] = fits.BinTableHDU.from_columns(cols + new_column)

        data.name = 'CONTINUUM'


        hdulist = fits.HDUList([spectrum, data])
        hdulist.info()
        print('\n\n')


        # TODO: update filename to save to

        hdulist.writeto('test.fits', overwrite=True)



    print('Saved to file: ', filename)
    print('\n\n\n')



if __name__ == '__main__':

    # sp = EdiblesSpectrum("/HD170740/RED_860/HD170740_w860_redl_20140915_O12.fits")
    # sp = EdiblesSpectrum("/HD23466/BLUE_346/HD23466_w346_blue_20180731_O11.fits")

    # subset = sp.getSpectrum(xmin=3270, xmax=3305)

    # # Continuum(x, y, method='spline', anchors=4)
    # model = Continuum(sp, method='spline', n_anchors=4)



    init_x = [1, 2, 3]
    init_y = [4004.56767336737, 5, 6]

    filename = 'HD23466_w346_blue_20180731_O11.fits'


    copy_to_FITS(filename=filename, name='init', x_points=init_x, y_points=init_y)






    next_x = [7, 8, 9]
    next_y = [34, 678, 23.646345]


    filename = 'test.fits'

    copy_to_FITS(filename=filename, name='test', x_points=next_x, y_points=next_y)
