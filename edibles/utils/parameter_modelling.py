import numpy as np
from scipy.special import wofz
from scipy.interpolate import interp1d
import astropy.constants as cst
import matplotlib.pyplot as plt
from edibles import EDIBLES_PYTHONDIR
from edibles.utils.edibles_oracle import EdiblesOracle
from edibles.utils.edibles_spectrum import EdiblesSpectrum
from pathlib import Path
import pandas as pd
from scipy.ndimage import gaussian_filter
from lmfit import Parameters, minimize,Model
import voigt_profile as vp

def file_reader (star_name):
    """
    Will read in any file from the voigt_benchmark folder once given the name of the file.
    Is used in order to read in data for o per, sigma sco, zeta per, and zeta oph stars in example 5
    :paramaters :
        star_name: string
            the file name containing the data that needs to be analysed
    :return:
        data["Wavelength"] : 1d ndarray
            Contains all the measured wavelength values for a particular file
        data["Normfluxval"] : 1d ndarray
            Contains normailised flux values corresponding to each wavelength reading
    """
    # define where the desired file is
    folder = EDIBLES_PYTHONDIR / "data"
    filename = folder / "voigt_benchmarkdata" / star_name

    # state what headers the desired data is under
    Headers = ["Wavelength", "Normfluxval"]

    # read in the data
    data = pd.read_csv(
        filename,
        delim_whitespace=True,
        skiprows=[0],
        header=None,
        names=Headers,
        engine="python",
    )
    return data["Wavelength"], data["Normfluxval"]

def wavegrid(b_array, N_array,v_rad_array, v_resolution_array, centre_lambda, offset, range, no_points):
    """
    Creates an array of wavelengths with predicted normalised flux vals for the parameters entered
    parameters
    ----------
        b_array: 1-d ndarray
            array containing calculated b values for a given sightline
        N_array: 1-d ndarray
            array containing calculated coluum density (N) values for a given sightline
        v_rad_array: 1-d ndarray
            array containing calculated radial velocity of the clouds in the given sightline
        v_resolution_array: 1-d ndarray
            array containg the resolutin velocity of the equipment used to obtain the data
        centre_lambda: float
            expected ideal absorption wavelength of molecule being studied
        offset: float
            used to adjust where the wavegrid array starts,
            allows us to see the model for data where the DIB is not exactly at centre lambda
        range: float
            desired range for the points to cover (as multiples of the number of data points from the data being fitted to)
        no_points: float
            number of points in the data that is being fitted to, this is done to allow easy data analyse later on.
             (e.g: reduced chisquared)

    returns
    ----------
        wavegrid_array : 1-d numpy ndarray
            array of wavelengths which an absoption line has been fitted to
        AbsorptionLine : 1-d numpy array
            contains the corresponding normalised flux values for the range of wavelengths contained in wavegrid_array
    """
    # create array of wavelength to fit an absoption line to
    wavegrid_array = np.arange(no_points) * range + (centre_lambda - offset)
    # calculate absoption line for wavelength using given parameters
    AbsorptionLine = vp.voigt_absorption_line(
        wavegrid_array,
        lambda0=centre_lambda,
        b=b_array,
        N=N_array,
        f=3.393e-1,
        gamma=3.8e7,
        v_rad=v_rad_array,
        v_resolution=v_resolution_array,
    )
    return wavegrid_array, AbsorptionLine

def reduced_chi_squared(observed_value, observed_error, expected_value):
    """
    calculates reduced chi squared of data entered

    parameters
    ----------
    observed_value : 1-D numpy array
        value from the real life datd
    observed_error : 1-D numpy array
        error on real data
    expected_value : 1-D numpy array
        value at ths point from fitted function

    returns
    ----------
   (np.sum(((observed_value - expected_value)/observed_error)**2))/(len(observed_value)-2
       reduced chi squared
    """
    return (np.sum(((observed_value - expected_value) / observed_error) ** 2)) / (len(observed_value) - 2)

def errors(low_lim, high_lim, data):
    """
    calculates the error (standard dev) of data entered within a certain slice of data defined using upper and lower limits.
    For this to work upper and lower limits must slice a section of the spectrum which is flat, this error will then be used for
    all data points. Whilst this is likely not the most accurate to calculate the error on the data it will suffice for our
    analyses of the data.

    parameters
    ----------
        low_lim: float
            lower limit of flat section of spectrum
        high_lim: float
            upper limit of flat section of spectrum
        data: 2-d numpy array
            contains the wavelength and flux data for the spectrum
    returns
    ----------
        error_array : 1-d numpy ndarray
            array of errors which will be used in the calculation of the reduced chi squared
    """
    # create an empty array for the errors to be stored in
    error_vals = []
    # compile all data points within the desired flat spectrum range into the error_vals array
    for i in range(len(data[0])):
        if (data[0,i] > low_lim and data[0,i] < high_lim):
            error_vals = np.append(error_vals, data[0, i])
    # get the standar dev of the points to use as the error
    error = np.std(error_vals)
    # assign every point this error
    error_array = np.array([error] * len(data[0]))
    return error_array

def o_per():
    """
    This function will use the parameters given in welty & hobbs, 2001, apjs 133, 345 to fit a voigt profile to data for
    the o per star as given in file in voigt_benchmarkdata. This model and the actual data are ploted on the same graph.
    reduced chi squared is then calculated for both of the models.
    currently attempting to create code which will gte accurate parameters from the data in the voigt_benchmarkdata
    """
    # arrays containing the Broadening parameters, coluum densities and radial velocitys, respectively for the different
    # components responsible for the spectrum
    b_o_per = [0.60, 0.44, 0.72, 0.62, 0.60]
    N_o_per = [12.5e10, 10e10, 44.3e10, 22.5e10, 3.9e10]
    v_rad_o_per = [10.5, 11.52, 13.45, 14.74, 15.72]

    # read in data and convert it from a tuple into a numpy array
    o_per_m = np.asarray(file_reader(files[0]))
    # fit a voigt profile to data using parameters defined above
    fit_o_per_m = wavegrid(b_o_per, N_o_per, v_rad_o_per, 0.56, 7698.974, 0, 0.002, len(o_per_m[0]))
    # calculate error between these wavelength which have been identified as the flat part of the spectrum by eye
    error_o_per_m = errors(7699.5, 7699.8, o_per_m)
    # calculate reduced chi squared of fit
    print('o_per_m reduced \u03C7\u00B2', reduced_chi_squared(o_per_m[1], error_o_per_m, fit_o_per_m[1]))

    # repeat process for same star but with data from another survey
    o_per_k = np.asarray(file_reader(files[1]))
    fit_o_per_k = wavegrid(b_o_per, N_o_per, v_rad_o_per, 1.40, 7698.974, 0, 0.002, len(o_per_k[0]))
    error_o_per_k = errors(7699.5, 7699.8, o_per_k)
    print('o_per_k reduced \u03C7\u00B2', reduced_chi_squared(o_per_k[1], error_o_per_k, fit_o_per_k[1]))

    fig_1 = plt.figure(figsize=(7, 8))
    axes = fig_1.add_subplot(211), fig_1.add_subplot(212)
    # plot data and corresponding fits on the same graph but on different subplots
    axes[0].plot(fit_o_per_m[0], fit_o_per_m[1], c='b', marker="1", label='fit')
    axes[1].plot(fit_o_per_k[0], fit_o_per_k[1], c='b', marker='1', label='fit')
    # set x limit to allow for better analyse of data by eye
    for i in range(2):
        axes[i].set_xlim(7699, 7699.75)
        axes[i].set_xlabel('Wavelength \u00C5')
        axes[i].set_ylabel('Normalised Flux')

    axes[0].plot(o_per_m[0], o_per_m[1], color="k", marker="D", fillstyle='none', markersize=0.5, label='o per (m95)')
    axes[1].plot(o_per_k[0], o_per_k[1], color="k", marker="D", fillstyle='none', markersize=0.5, label='o per (k94)')
    fig_1.legend()

    # attempting to obtain parameter values from the data
    parameters = Parameters()
    #parms = {'b': 0.6, 'N': 12.5e10, 'v_rad':10.5}
    parameters.add('b', value = 0.6)
    parameters.add('N', value = 12.5e10)
    parameters.add('v_rad', value = 10.5)
    #parameters.add('v_resolution', value = 0.56, vary=False)
    #parameters.add('f', value=3.393e-1, vary=False)
    #parameters.add('gamma', value=3.8e7, vary=False)
    #parameters.add('lambda0', value=7698.974, vary=False)
    #parameters.add('n_step', value=25, vary=False)
    #parameters.add('v_res', value=0.56)
    #parameters.add('f', value=3.393e-1)
    #parameters.add('gamma', value=3.8e7)
    #parameters.add('lambda0', value=7698.974)

    #parameters.add('v_res', min=0.3, max=1.5)

    #model = Model(vp.voigt_absorption_line)
    #fit = model.fit(o_per_m[1], parameters)
    #wavegrid, lambda0=0.0, f=0.0, gamma=0.0, b=0.0, N=0.0, v_rad=0.0, v_resolution=0.0
    #print(fit)
    #print(type(parameters))

    model = Model(vp.voigt_absorption_line, indepedent_vars = ['wavegrid','f', 'gamma', 'v_resolution', 'lambda0'], param_names=['b','N','v_rad','f', 'gamma', 'v_resolution', 'lambda0'])
    params = model.make_params(#wavegrid={'value': o_per_m[0], 'vary':False},
                               #f={'value': 3.393e-1, 'vary':False},
                               #gamma= {'value':3.8e7 , 'vary':False},
                               #v_resolution ={'value': 0.56, 'vary':False},
                               #lambda0= {'value': 7698.974, 'vary':False},
                               #n_step = {'value': 25, 'vary':False},
                               b=0.6, N=12.5e10, v_rad=10.5)
    #print(o_per_m[0])
    fit = model.fit(o_per_m[1], params, wavegrid=o_per_m[0])
    print(fit)


    #reduced_chi = 0
    #runs = 1
    #number_of_components = 1
    #while (reduced_chi > 5 or reduced_chi < 0.5):

    #    b_array = np.array([1])
    #    N_array = np.array([0] * number_of_components)
    #    v_rad_array = np.array([0] * number_of_components)
    #    v_res_array = np.array([0] * number_of_components)

    #    parameters.add('b', min = 0.3, max = 1.5 )
        #parameters.add('b', min = 0.2, max = 1.5 )
    #    print(parameters['b'])
    #    parameters.add('N', min=0.1e10, max=50e10)
        #parameters.add('N', min=0.1e10, max=50e10)
    #    parameters.add('v_rad', min=-20, max=20)
        #parameters.add('v_rad', min=-20, max=20)
    #    parameters.add('v_res', min=0.3, max=1.5)
        #parameters.add('v_res', min=0.3, max=1.5)

   #     fit = minimize(wavegrid_lmfit, parameters)
    #    print(fit.params['b'])
   #     b = fit.params['b'].value
   #     N = fit.params['N'].value
   #     v_rad = fit.params['v_rad'].value
    #    v_res = fit.params['v_res'].value
    #    fitted_vals = wavegrid(b, N, v_rad, v_res, 7698.974, 0, 0.002, 399)
    #    reduced_chi = reduced_chi_squared(o_per_m[1], error_o_per_m, fitted_vals[1])
    #    #print('reduced_chi', reduced_chi)
    #    runs += 1
    #print(b)
    #print(N)
    #print(v_rad)
    #print(v_res)
    #print(runs)
    #axes[0].plot(fitted_vals[0], fitted_vals[1], label = 'minimised fit ', c= 'r')
    #print('reduced_chi', reduced_chi)

def sigsco():
    """
    This function will use the parameters given in welty & hobbs, 2001, apjs 133, 345 to fit a voigt profile to data for
    the sigma sco star as given in file in voigt_benchmarkdata. This model and the actual data are ploted on the same graph.
    reduced chi squared is then calculated for both of the models.
    """
    # arrays containing the Broadening parameters, coluum densities and radial velocitys, respectively for the different
    # components responsible for the spectrum
    b_sigsco = [0.70, 1.23, 0.75, 0.54, 0.60]
    N_sigsco = np.array([0.4e10, 0.8e10, 8.1e10, 3.8e10, 0.6e10])
    v_rad_sigsco = [-14.23, -8.75, -6.26, -4.62, -3.13]
    # read in data and convert it from a tuple into a numpy array
    sigsco = np.asarray(file_reader(files[2]))
    # fit a voigt profile to data using parameters defined above
    fit_sigsco = wavegrid(b_sigsco, N_sigsco, v_rad_sigsco, 1.20, 7698.974, 0.6, 0.002, len(sigsco[0]))
    # calculate error between these wavelength which have been identified as the flat part of the spectrum by eye
    error_sigsco = errors(7699.0, 7699.4, sigsco)
    # calculate reduced chi squared of fit
    print('sigsco reduced \u03C7\u00B2', reduced_chi_squared(sigsco[1], error_sigsco, fit_sigsco[1]))

    # plot data and corresponding fits on the same graph but on different subplots
    fig_2 = plt.figure()
    ax = fig_2.add_subplot(111)
    ax.plot(sigsco[0], sigsco[1], color="k", marker="D", fillstyle='none', markersize=0.5, label='sigsco')
    ax.plot(fit_sigsco[0], fit_sigsco[1], c='b', label='fit')
    # set x limit to allow for better analyse of data by eye
    ax.set_xlim(7698, 7699.5)
    ax.set_xlabel('Wavelength \u00C5')
    ax.set_ylabel('Normalised Flux')
    fig_2.legend()

def zetoph():
    """
    This function will use the parameters given in welty & hobbs, 2001, apjs 133, 345 to fit a voigt profile to data for
    the zeta oph star as given in file in voigt_benchmarkdata. This model and the actual data are ploted on the same graph.
    reduced chi squared is then calculated for both of the models.
    """
    # arrays containing the Broadening parameters, coluum densities and radial velocitys, respectively for the different
    # components responsible for the spectrum
    b_zetoph = [0.96, 0.80, 0.57, 0.43, 0.58]
    N_zetoph = np.array([1e10, 1.2e10, 40.9e10, 27.2e10, 1.1e10])
    v_rad_zetoph = [-19.09, -16.50, -14.98, -13.96, -12.73]
    # read in data and convert it from a tuple into a numpy array
    zetoph_k = np.asarray(file_reader(files[3]))
    # fit a voigt profile to data using parameters defined above
    fit_zetoph_k = wavegrid(b_zetoph, N_zetoph, v_rad_zetoph, 1.40, 7698.974, 0.55, 0.0015, len(zetoph_k[0]))
    # calculate error between these wavelength which have been identified as the flat part of the spectrum by eye
    error_zetoph_k = errors(7698.8, 7699, zetoph_k)
    # calculate reduced chi squared of fit
    print('Zetoph_k reduced \u03C7\u00B2', reduced_chi_squared(zetoph_k[1], error_zetoph_k, fit_zetoph_k[1]))

    # repeat process for same star but with data from another survey
    zetoph_lf = np.asarray(file_reader(files[4]))
    fit_zetoph_lf = wavegrid(b_zetoph, N_zetoph, v_rad_zetoph, 0.40, 7698.974, 0.6, 0.0014, len(zetoph_lf[0]))
    error_zetoph_lf = errors(7698.8, 7699, zetoph_lf)
    print('Zetoph_lf reduced \u03C7\u00B2', reduced_chi_squared(zetoph_lf[1], error_zetoph_lf, fit_zetoph_lf[1]))

    # repeat process for same star but with data from another survey
    zetoph_m = np.asarray(file_reader(files[5]))
    fit_zetoph_m = wavegrid(b_zetoph, N_zetoph, v_rad_zetoph, 0.56, 7698.974, 0.55, 0.0015, len(zetoph_m[0]))
    error_zetoph_m = errors(7698.8, 7699, zetoph_m)
    print('Zetoph_m reduced \u03C7\u00B2', reduced_chi_squared(zetoph_m[1], error_zetoph_m, fit_zetoph_m[1]))

    # plot data and corresponding fits on the same graph but on different subplots
    fig_3 = plt.figure()
    axes = fig_3.add_subplot(311), fig_3.add_subplot(312), fig_3.add_subplot(313)
    axes[0].plot(zetoph_m[0], zetoph_m[1], color="k", marker="D", fillstyle='none', markersize=0.5,
                 label='zetoph')
    axes[0].plot(fit_zetoph_m[0], fit_zetoph_m[1], c='b', label='fit')

    axes[1].plot(zetoph_lf[0], zetoph_lf[1], color="k", marker="D", fillstyle='none', markersize=0.5,
                 label='zetoph_2')
    axes[1].plot(fit_zetoph_lf[0], fit_zetoph_lf[1], c='b', label='fit')

    axes[2].plot(zetoph_k[0], zetoph_k[1], color="k", marker="D", fillstyle='none', markersize=0.5,
                 label='zetoph_3')
    axes[2].plot(fit_zetoph_k[0], fit_zetoph_k[1], c='b', label='fit')
    # set x limit to allow for better analyse of data by eye
    for i in range(3):
        axes[i].set_xlim(7698.4, 7699.4)
        axes[i].set_xlabel('Wavelength \u00C5')
        axes[i].set_ylabel('Normalised Flux')
    fig_3.legend()

def zetper():
    """
    This function will use the parameters given in welty & hobbs, 2001, apjs 133, 345 to fit a voigt profile to data for
    the zeta per star as given in file in voigt_benchmarkdata. This model and the actual data are ploted on the same graph.
    reduced chi squared is then calculated for both of the models.
    """
    # arrays containing the Broadening parameters, coluum densities and radial velocitys, respectively for the different
    # components responsible for the spectrum
    b_zetper = [1, 0.68, 0.8, 0.5]
    N_zetper = [2.2e10, 26.0e10, 42.7e10, 0.1e10]
    v_rad_zetper = [11.48, 13.25, 14.54, 16.40]
    # read in data and convert it from a tuple into a numpy array
    zetper = np.asarray(file_reader(files[6]))
    # fit a voigt profile to data using parameters defined above
    fit_zetper = wavegrid(b_zetper, N_zetper, v_rad_zetper, 0.56, 7698.974, 0.1, 0.002, len(zetper[0]))
    # calculate error between these wavelength which have been identified as the flat part of the spectrum by eye
    error_zetper = errors(7699.5, 7700, zetper)
    # calculate reduced chi squared of fit
    print('Zetper reduced \u03C7\u00B2', reduced_chi_squared(zetper[1], error_zetper, fit_zetper[1]))
    # plot data and corresponding fits on the same graph but on different subplots
    fig_4 = plt.figure()
    ax = fig_4.add_subplot(111)
    ax.plot(zetper[0], zetper[1], color="k", marker="D", fillstyle='none', markersize=0.5, label='zetper')
    ax.plot(fit_zetper[0], fit_zetper[1], c='b')
    # set x limit to allow for better analyse of data by eye
    ax.set_xlim(7698.5, 7700)
    ax.set_xlabel('Wavelength \u00C5')
    ax.set_ylabel('Normalised Flux')
    fig_4.legend()


#############################################################
#
# EXAMPLE 5: analysis of  interstellar Potassium line in 4 different stars using data from multiple surveys
# attempts to also find parameters from the raw data through the use of minimisation functions
#
#############################################################

# array of files with the data that will be studied in the example
files = ["omiper.m95.7698.txt", "omiper.k94.7698.txt", "sigsco.k00a.7698.txt", "zetoph.k94.7698.txt",
         "zetoph.lf.7698.txt", "zetoph.m95.7698.txt", "zetper.m95.7698.txt"]

# run analyse of the 4 stars
o_per()
#sigsco()
#zetoph()
#zetper()


plt.show()



