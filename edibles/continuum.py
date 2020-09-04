import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
import os
from datetime import datetime

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

        self.model = ContinuumModel(n_anchors=n_anchors)
        cont_pars = self.model.guess(self.Spectrum.flux, x=self.Spectrum.wave)

        self.result = self.model.fit(
            data=self.Spectrum.flux,
            params=cont_pars,
            x=self.Spectrum.wave
        )

        print(self.result.params)

        self.result.plot_fit()
        plt.show()

    def alphashape(self):
        print("method: ", self.method)

    def polynomial(self):
        print("method: ", self.method)


def copy_to_FITS(filename, x_points, y_points, comment=False):
    # filename = "/Users/haoyufan/test.fits"

    # TODO: automatically read in filename, from edibles-spectrum class?
    hdulist = fits.open(filename)
    spectrum = hdulist['PRIMARY']

    # create continuum_hdu to add, by initializing or appending to existing data
    try:
        continuum_hdu = hdulist["SPLINE-CONTINUUM"]
        continuum_header = continuum_hdu.header
        print("Adding NEW spline anchor points")
        # data = continuum_hdu.data
        batch_idx = str(len(continuum_hdu.data.columns) // 2)
        new_column = fits.ColDefs(
            [fits.Column(name='x_' + batch_idx, format='D', array=x_points),
             fits.Column(name='y_' + batch_idx, format='D', array=y_points)]
        )
        continuum_hdu = fits.BinTableHDU.from_columns(continuum_hdu.data.columns + new_column)
        continuum_hdu.header = continuum_header
    except:
        print('Adding INITIAL spline anchor points.')
        batch_idx = "0"
        columns = fits.ColDefs(
            [fits.Column(name='x_' + batch_idx, format='D10.4', array=x_points),
             fits.Column(name='y_' + batch_idx, format='D10.4', array=y_points)]
        )
        continuum_hdu = fits.BinTableHDU.from_columns(columns)
        continuum_hdu.header["DICT"] = " "

    # edit header of continuum_hdu
    continuum_hdu.name = "SPLINE-CONTINUUM"
    user_name = input("Please type your user name:")
    continuum_hdu.header["USER_" + batch_idx] = user_name
    continuum_hdu.header["DATE_" + batch_idx] = datetime.today().strftime('%Y-%m-%d')
    if comment:
        comment_str = input("Please type your comment:")
    else:
        comment_str = ""
    continuum_hdu.header['CMT_' + batch_idx] = comment_str
    continuum_hdu.header["DICT"] = continuum_hdu.header["DICT"] + \
                                   batch_idx + ":" + user_name + "-" + datetime.today().strftime('%Y')

    # create new hdulist and write to file
    # TODO: update filename to save to, read from edibles_spectrum class?
    if len(hdulist) == 1:
        hdulist.append(continuum_hdu)
    else:
        hdulist[1] = continuum_hdu
    hdulist.writeto(filename, overwrite=True)
    hdulist.info()
    hdulist.close()


if __name__ == '__main__':
    # sp = EdiblesSpectrum("/HD170740/RED_860/HD170740_w860_redl_20140915_O12.fits")
    sp = EdiblesSpectrum("/HD23466/BLUE_346/HD23466_w346_blue_20180731_O11.fits")

    subset = sp.getSpectrum(xmin=3270, xmax=3305)

    # Continuum(x, y, method='spline', anchors=4)
    cont = Continuum(sp, method='spline', n_anchors=4)

    print(cont.result.params)



    # init_x = [1, 2, 3]
    # init_y = [4004.56767336737, 5, 6]

    # filename = 'HD23466_w346_blue_20180731_O11.fits'

    # copy_to_FITS(filename=filename, name='init', x_points=init_x, y_points=init_y)

    # next_x = [7, 8, 9]
    # next_y = [34, 678, 23.646345]

    # filename = 'test.fits'

    # copy_to_FITS(filename=filename, name='test', x_points=next_x, y_points=next_y)
