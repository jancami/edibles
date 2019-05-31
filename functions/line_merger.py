import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits


def line_merger(file1,xmin1,xmax1,file2,xmin2,xmax2):




    hdu1 = fits.open(file1)

    spec_flux1 = hdu1[0].data
    crval1 = hdu1[0].header["CRVAL1"]
    cdelt1 = hdu1[0].header["CDELT1"]
    nwave = len(spec_flux1)
    wave = np.arange(0, nwave, 1)
    spec_wave1 = (wave) * cdelt1 + crval1

    fig1 = plt.plot(spec_wave1, spec_flux1)
    print(plt.show(fig1))

    min_idx1 = (np.abs(spec_wave1 - xmin1)).argmin()
    max_idx1 = (np.abs(spec_wave1 - xmax1)).argmin()
    wave_subset1 = spec_wave1[min_idx1:max_idx1]
    flux_subset1 = spec_flux1[min_idx1:max_idx1]

    


    hdu2 = fits.open(file2)

    spec_flux2 = hdu2[0].data
    crval1 = hdu2[0].header["CRVAL1"]
    cdelt1 = hdu2[0].header["CDELT1"]
    nwave = len(spec_flux2)
    wave = np.arange(0, nwave, 1)
    spec_wave2 = (wave) * cdelt1 + crval1

    fig2 = plt.plot(spec_wave2, spec_flux2)
    print(plt.show(fig2))

    min_idx2 = (np.abs(spec_wave2 - xmin2)).argmin()
    max_idx2 = (np.abs(spec_wave2 - xmax2)).argmin()
    wave_subset2 = spec_wave2[min_idx2:max_idx2]
    flux_subset2 = spec_flux2[min_idx2:max_idx2]

    # fig1 = plt.plot(wave_subset1,flux_subset1)
    # fig2 = plt.plot(wave_subset2,flux_subset2)
    # print(plt.show(fig1),plt.show(fig2))


    x_diff = xmin2 - xmax1
    wave_subset1 += x_diff

    y_diff = spec_flux1[max_idx1] - spec_flux2[min_idx2]
    flux_subset1 -= y_diff

    fig1 = plt.plot(wave_subset1,flux_subset1)
    fig2 = plt.plot(wave_subset2,flux_subset2)
    print(plt.show(fig1),plt.show(fig2))

    wave_subset = np.append( wave_subset1 , wave_subset2)
    # wave_subset.append(wave_subset2)

    flux_subset = np.append(flux_subset1 , flux_subset2)
    # flux_subset.append(flux_subset2)

    # plt.plot(wave_subset,flux_subset)
    # plt.show()

    return wave_subset, flux_subset



if __name__ == "__main__":

    file1 = '/data/DR3_fits/HD170740/BLUE_346/HD170740_w346_n6_20160612_B.fits'
    xmin1 = 3300.
    xmax1 = 3305.

    file2 = '/data/DR3_fits/HD170740/RED_564/HD170740_w564_n9_20160612_U.fits'
    xmin2 = 5885.
    xmax2 = 5898.

    wave_subset, flux_subset = line_merger(file1,xmin1,xmax1,file2,xmin2,xmax2)

    plt.plot(wave_subset,flux_subset)
    plt.show()
