from astropy.io import fits

from edibles.edibles import DATADIR


def printHeader(input_fits):
    """
    A function to print out the header of a FITS file

    Parameters
    ----------

    input_fits : str
        path to FITS file starting from DATADIR


    Returns
    -------
    Prints the header.
    
    """

    path = DATADIR + input_fits

    hdu = fits.open(path)
    print(hdu[0].header)


if __name__ == "__main__":

    input_fits = "/HD170740/RED_860/HD170740_w860_redl_20140915_O12.fits"

    printHeader(input_fits)
