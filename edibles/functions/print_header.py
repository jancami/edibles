from astropy.io import fits


def print_header(input_fits):
    '''
    usage:
    header("/path/to/fits_file.fits")
    '''

    hdu = fits.open(input_fits)
    print(hdu[0].header)

if __name__ == '__main__':
	print_header('/data/DR3_fits/HD170740/RED_860/HD170740_w860_n12_20160613_L.fits')
