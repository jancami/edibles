import numpy as np
import matplotlib.pyplot as plt
from edibles.utils import voigt_profile as vp
from edibles.utils.benchmark_voigt.parameter_modelling import file_reader
import astropy.constants as cst
import scipy.signal as ss
from edibles import PYTHONDIR
from pathlib import Path
import pandas as pd

def find_nearest(array, value):
    """
    finds the value in an array which is closest to the desired value
    paramaters :
        array : list/array/numpy array
            an array of vaules which we plan to extrcat a single values from
        value : int
            the value which we wish to be close to
    return:
        idx : int
            the position in the array of the nearest value
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


#define path to file which contains the names of each star and the resolution of intraments used in survey
folder = Path(PYTHONDIR + "/data")
filename = folder / "voigt_benchmarkdata" / 'parameter_modelling_data' / "files_for_auto_voigt.txt"

# state what headers the desired data is under
Headers = ["star_name","resolution"]

# read in the data
file = pd.read_csv(
    filename,
     delim_whitespace=True,
     header=None,
     names=Headers,
     engine="python",
    )

# collect the data from each header in easy to iterate form
resolution = file['resolution']
star_name = file["star_name"]
central_wavelength = 7698.974

# loop over all data
for i in range(len(star_name)):

    # extract the raw data from the star file
    wavelengths, normflux = np.asarray(file_reader(star_name[i]+'.7698.txt'))
    print(star_name[i])

    # asign guesses to b and N values, value does not matter as the fitting routine should fix this
    # so we are more worried about getting an accurate v_rad value
    b= np.array([1.00])
    N = np.array([1e12])

    # need good initial guess for v_rad
    # calc flux weighted wavelength and convert to v_rad
    wave_weight = sum((1 - normflux) * wavelengths / sum(1-normflux))
    #print(wave_weight)

    # convert the guess wave position into velocity using eq below
    v_rad_guess = (wave_weight - central_wavelength)/central_wavelength * cst.c.to("km/s").value
    #print(v_rad_guess)
    v_rad = np.array([v_rad_guess])

    # assuming that the first 50 values in the normflux array are part of the continuum
    # and have used the standard dev of these points as an overall uncertainty on all data collected
    errors = np.std(normflux[:50])

    # use function to fit 1 component to data
    fit = vp.fit_multi_voigt_absorptionlines(wavegrid=wavelengths,
                                        ydata=normflux,
                                        restwave=central_wavelength,
                                        f=3.393e-1,
                                        gamma=3.8e7,
                                        b=b,
                                        N=N,
                                        v_rad=v_rad,
                                        v_resolution= resolution[i],
                                        n_step=25,
                                        std_dev= errors)
    # plotting results

    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    #ax.plot(wavelengths, normflux, c = 'k', label = 'data')
    #ax.plot(wavelengths, fit.best_fit, 'b', label = '1 component')
    #ax.plot(wavelengths[:50], normflux[:50], c= 'orange', marker= '*', markersize = 5)
    #ax.plot(wavelengths, normflux/fit.best_fit, c='r', label = 'flux/best fit')
    #plt.show()

    # assigning 'previous_bic' (bic = bayesian info criterium) this may be used in
    # component fitting in order to tell our program when we have a good fit
    previous_bic = fit.summary()['bic']
    previous_aic = fit.summary()['aic']

    # define reduced chi and delta_bic to allow us to enter the loop
    previous_red_chi = fit.summary()['redchi']
    delta_bic = previous_bic
    delta_aic = previous_aic
    delta_red_chi = -1

    # minimum accepted previous_red_chi for working data
    # omiper.m95 previous_red_chi > 1
    # zetoph.lf previous_red_chi > 1.18
    # zetoph.k94 previous_red_chi > 1.47
    while (delta_bic < 0) | (delta_aic < 0) | (previous_red_chi > 1.47) | (delta_red_chi < 0):

        # retrieve the previous reduced chi squared
        previous_red_chi = fit.summary()['redchi']

        previous_fit = fit

        # array of residuals obtained by dividing the normflux by the flux from the last fit
        residual_array = normflux/fit.best_fit

        # find the peak position of the residuals plot
        peak_position = wavelengths[np.argwhere(residual_array == np.max(residual_array))[0,0]]

        # create an array of all minima in the residuals plot
        mins = wavelengths[ss.argrelextrema(residual_array, np.less)]

        #print('peak ', peak_position)
        #print('mins', find_nearest(mins[np.argwhere(mins > peak_position)], peak_position))

        # find where the 2 nearest local minimum to the peak is
        lower_near_lambda = mins[find_nearest(mins[np.argwhere(mins < peak_position)], peak_position)]
        upper_near_lambda = mins[find_nearest(mins[np.argwhere(mins > peak_position)], peak_position) + np.argwhere(mins > peak_position)[0,0]]

        #print('lower lambda', lower_near_lambda)
        #print('upper lambda', upper_near_lambda)

        # convert the previously fitted v_rads into wavelength for easy comparison with the upper and lower minimas
        min_lambdas = (v_rad * central_wavelength/ cst.c.to("km/s").value) + central_wavelength

        # replace the v_rad nearest to the peak of the residuals data with the lower minima from above and then add
        # the upper minima to the  end of the array and sort the array into size order using .sort
        min_lambdas[find_nearest(min_lambdas, peak_position)] = lower_near_lambda
        min_lambdas = np.append(min_lambdas, upper_near_lambda)
        min_lambdas.sort()

        # convert all lambdas back to velocities for fitting
        v_rad = (min_lambdas - central_wavelength)/central_wavelength * cst.c.to("km/s").value

        # make sure the b and N arrays are the right length
        b = np.array([1]* len(v_rad))
        N = np.array([1e12]* len(v_rad))

        # check that we are adding the right amount of components
        print('components to be fitted', len(v_rad))

        # fit all of the components
        fit = vp.fit_multi_voigt_absorptionlines(wavegrid=wavelengths,
                                                   ydata=normflux,
                                                   restwave=central_wavelength,
                                                   f=3.393e-1,
                                                   gamma=3.8e7,
                                                   b=b,
                                                   N=N,
                                                   v_rad=v_rad,
                                                   v_resolution=resolution[i],
                                                   n_step=25,
                                                   std_dev = errors)

        #fit.params.pretty_print()

        # calculate the difference in the bic from this fit and the previous fit, when this value is positive it means
        # that we have added too many components and should use the fit before this as th ideal fit
        delta_bic= fit.summary()['bic'] - previous_bic
        previous_bic = fit.summary()['bic']
        delta_red_chi = fit.summary()['redchi'] - previous_red_chi
        delta_aic = fit.summary()['aic'] - previous_aic
        previous_aic = fit.summary()['aic']

        #print('reduced chi-square: ',fit.summary()['redchi'])
        #print('bic', fit.summary()['bic'])
        #print('delta bic', fit.summary()['bic'] - previous_bic)

        # This is used for debugging just to check what is being fitted and when
        residual_array = normflux / fit.best_fit

        # find the peak position of the residuals plot
        print('peak',wavelengths[np.argwhere(residual_array == np.max(residual_array))[0, 0]])
        print('lambdas',(v_rad * central_wavelength / cst.c.to("km/s").value) + central_wavelength)
        print('reduced chi-square: ', fit.summary()['redchi'])
        print('delta bic', delta_bic)
        print('delta aic', delta_aic)

        #print(fit.fit_report())

        #fig = plt.figure()
        #ax = fig.add_subplot(111)
        #ax.plot(wavelengths, normflux, mfc = 'none',  c= 'k', marker = 'D', label = 'Data')
        #ax.plot(wavelengths, fit.best_fit, c = 'b', label = 'Fit')
        #ax.plot(wavelengths, normflux/fit.best_fit, label = 'residuals')
        #ax.plot(np.array([wavelengths[np.argwhere(residual_array == np.max(residual_array))[0, 0]]] * 50), np.linspace(0, 2, 50), linestyle= 'dashed')
        #plt.legend()
        #plt.show()

    # using the previous fit instead of the most recent one due use of bic
    fit = previous_fit

    # stating the number of components, this seemed like the easiest way to calculate it
    components = (len(b) - 1)
    b = np.zeros(components)
    b_err = np.zeros(components)
    N = np.zeros(components)
    N_err = np.zeros(components)
    v_rad = np.zeros(components)
    v_rad_err = np.zeros(components)

    for j in range(components):
        b[j] = fit.params[f'b{j}'].value
        b_err[j] = fit.params[f'b{j}'].stderr
        N[j] = fit.params[f'N{j}'].value
        N_err[j] = fit.params[f'N{j}'].stderr
        v_rad[j] = fit.params[f'v_rad{j}'].value
        v_rad_err[j] = fit.params[f'v_rad{j}'].stderr

    # print results to see if we are satisfied
    print('final no components:', components)
    print('reduced chi-square: ',fit.summary()['redchi'])
    #print('bic', fit.summary()['bic'])
    print('delta bic', delta_bic)
    #fit.params.pretty_print()

    # plot final fit alongside the raw data
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(wavelengths, normflux, mfc = 'none',  c= 'k', marker = 'D', label = 'Data')
    ax.plot(wavelengths, fit.best_fit, c = 'b', label = 'Fit')
    ax.plot(wavelengths, normflux/fit.best_fit, label = 'residuals')
    fig.suptitle(star_name[i])
    plt.legend()

    # get parameters from Welty and Hobbs paper for comparison
    folder = Path(PYTHONDIR + "/data")
    filename = folder / "voigt_benchmarkdata" / 'parameter_modelling_data' / (star_name[i] + '.txt')

    Headers = ["b", "N", "v_rad"]
    star_parameters = pd.read_csv(filename,
                                  delim_whitespace=True,
                                  header=None,
                                  names=Headers,
                                  engine="python",
                                  )

    # collect the data from each header in easy to iterate form
    b_WH = np.asarray(star_parameters["b"])
    N_WH = np.asarray(star_parameters["N"])
    v_rad_WH = np.asarray(star_parameters["v_rad"])

    # this makes sure that the arrays are the same size as the fitted array
    # just to allow for comparisons whilst bug fixing
    while (components - len(b_WH)) > 0:
        b_WH = np.append(b_WH,0)
        N_WH = np.append(N_WH, 0)
        v_rad_WH = np.append(v_rad_WH,0)

    while (len(b_WH) - components) > 0:
        new = len(b_WH)
        fit.params.add(f'b{new}', value = 0, vary = False)
        fit.params.add(f'N{new}', value = 0, vary = False)
        fit.params.add(f'v_rad{new}', value = 0, vary = False)

    # set the v_rad into increasing order and make sure the associated b and N values
    # move with them to the correct positions in their respective arrays
    v_rad_order = v_rad
    v_rad_err_order = v_rad_err
    b_order = b
    b_err_order = b_err
    N_order = N
    N_err_order = N_err
    v_rad.sort()
    # reorder the arrays
    for j in range(components):
        b[np.argwhere(v_rad == v_rad_order[j])[0,0]] = b_order[j]
        b_err[np.argwhere(v_rad == v_rad_order[j])[0, 0]] = b_err_order[j]
        N[np.argwhere(v_rad == v_rad_order[j])[0,0]] = N_order[j]
        N_err[np.argwhere(v_rad == v_rad_order[j])[0, 0]] = N_err_order[j]
        v_rad_err[np.argwhere(v_rad == v_rad_order[j])[0, 0]] = v_rad_err_order[j]


    # create and empty array of the right shape for our table that allows strings as an input
    table_data = np.asarray([['                                           ' for x in range(2)] for y in range(components* 3)])
    row_labels = np.array(['           ' for x in range( components * 3)])
    col_labels = ['Welty & Hobbs', 'Fitting components']

    # fill table data with all necessary data that we want on the table
    for x in range(2):
        if x == 0:
            for y in range(components):
                table_data[y][x] = f'{b_WH[y]:#.3f}'
                table_data[y + components][x] = f'{N_WH[y]:#.3g}'
                table_data[y + components * 2][x] = f'{v_rad_WH[y]:#.3f}'
                row_labels[y] = f'b{y}'
                row_labels[y + components] = f'N{y}'
                row_labels[y + components * 2] = f'v_rad{y}'
        else:
            for y in range(components):
                table_data[y][x] = '{0:#.3f} \u00B1 {1:#.3f}'.format(b[y], b_err[y])
                table_data[y + components][x] = '{0:#.3g} \u00B1 {1:#.3g}'.format(N[y],N_err[y])
                table_data[y + components * 2][x] = '{0:#.3f} \u00B1 {1:#.3f}'.format(v_rad[y], v_rad_err[y])

    # show table
    fig_3 = plt.figure()
    ax = fig_3.add_subplot(111)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    table = ax.table(table_data, rowLabels = row_labels, colLabels = col_labels, loc = 'center')
    plt.box(on=False)
    table.scale(1,1)
    fig_3.suptitle(star_name[i])
    plt.show()




