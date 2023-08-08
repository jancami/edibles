import numpy as np
import matplotlib.pyplot as plt
from edibles import PYTHONDIR
from pathlib import Path
import pandas as pd
import voigt_profile as vp
import os.path

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
    folder = Path(PYTHONDIR + "/data")
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

def errors(low_lim, high_lim, wavelength, normflux):
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
    for i in range(len(wavelength)):
        if (wavelength[i] > low_lim and wavelength[i]< high_lim):
            error_vals = np.append(error_vals, normflux)
    # get the standar dev of the points to use as the error
    error = np.std(error_vals)
    # assign every point this error
    error_array = np.array([error] * len(wavelength))
    return error_array

def compare_data(file_name,data, b, N, v_rad):


    save_path = 'c:/Users/user/edibles/edibles/data/voigt_benchmarkdata/parameter_modelling_data/'
    file_path = os.path.join(save_path, file_name)
    file = open(file_path, 'w+')


    txt = 'welty & hobbs parameteres: \n b: {0}, \n N: {1}, \n v_rad: {2} \n' \
          'Fit parameters: {3}'.format(b,N,v_rad,data.fit_report())
    file.write(txt)
    file.close()

def o_per():
    """
    This function will use the parameters given in welty & hobbs, 2001, apjs 133, 345 to fit a voigt profile to data for
    the o per star as given in file in voigt_benchmarkdata. This model and the actual data are ploted on the same graph.
    reduced chi squared is then calculated for both of the models.
    currently attempting to create code which will gte accurate parameters from the data in the voigt_benchmarkdata
    """
    # Create figure early on to allow for plotting throughout the function, allows for a nicer set out
    fig_1 = plt.figure(figsize=(7, 8))
    axes = fig_1.add_subplot(211), fig_1.add_subplot(212)
    # arrays containing the Broadening parameters, coluum densities and radial velocitys, respectively for the different
    # components responsible for the spectrum
    b_o_per = [0.60, 0.44, 0.72, 0.62, 0.60]
    N_o_per = [12.5e10, 10e10, 44.3e10, 22.5e10, 3.9e10]
    v_rad_o_per = [10.5, 11.52, 13.45, 14.74, 15.72]

    # read in data and convert it from a tuple into a numpy array
    o_per_m_wavelengths, o_per_m_normflux = np.asarray(file_reader(files[0]))
    # fit a voigt profile to data using parameters defined above
    fit_flux_o_per_m = vp.voigt_absorption_line(
        o_per_m_wavelengths,
        lambda0=7698.974,
        b=b_o_per,
        N=N_o_per,
        f=3.393e-1,
        gamma=3.8e7,
        v_rad=v_rad_o_per,
        v_resolution=0.56,
    )

    multi_fit_flux_m = vp.fit_multi_voigt_absorptionlines(wavegrid= o_per_m_wavelengths,
                                       ydata= o_per_m_normflux,
                                       restwave= 7698.974,
                                       f= 3.393e-1,
                                       gamma= 3.8e7,
                                       b= b_o_per,
                                       N= N_o_per,
                                       v_rad= v_rad_o_per,
                                       v_resolution= 0.56,
                                       n_step= 25)
    # calculate error between these wavelength which have been identified as the flat part of the spectrum by eye
    error_o_per_m = errors(7699.5, 7699.8, o_per_m_wavelengths, o_per_m_normflux)
    # calculate reduced chi squared of fit
    print('o_per_m reduced \u03C7\u00B2', reduced_chi_squared(o_per_m_normflux, error_o_per_m, fit_flux_o_per_m))

    # repeat process for same star but with data from another survey
    o_per_k_wavelengths, o_per_k_normflux = np.asarray(file_reader(files[1]))
    fit_flux_o_per_k = vp.voigt_absorption_line(
        o_per_k_wavelengths,
        lambda0=7698.974,
        b=b_o_per,
        N=N_o_per,
        f=3.393e-1,
        gamma=3.8e7,
        v_rad=v_rad_o_per,
        v_resolution=1.40,
    )


    multi_fit_flux_k = vp.fit_multi_voigt_absorptionlines(wavegrid= o_per_k_wavelengths,
                                       ydata= o_per_k_normflux,
                                       restwave= 7698.974,
                                       f= 3.393e-1,
                                       gamma= 3.8e7,
                                       b= b_o_per,
                                       N= N_o_per,
                                       v_rad= v_rad_o_per,
                                       v_resolution= 1.40,
                                       n_step= 25)

    error_o_per_k = errors(7699.5, 7699.8, o_per_k_wavelengths, o_per_k_normflux)
    print('o_per_k reduced \u03C7\u00B2', reduced_chi_squared(o_per_k_normflux, error_o_per_k, fit_flux_o_per_k))

    # plot data and corresponding fits on the same graph but on different subplots
    axes[0].plot(o_per_m_wavelengths, fit_flux_o_per_m, c='b', marker="1", label='fit')
    axes[1].plot(o_per_k_wavelengths, fit_flux_o_per_k, c='b', marker='1', label='fit')
    # set x limit to allow for better analyse of data by eye
    for i in range(2):
        axes[i].set_xlim(7699, 7699.75)
        axes[i].set_xlabel('Wavelength \u00C5')
        axes[i].set_ylabel('Normalised Flux')

    axes[0].plot(o_per_m_wavelengths, o_per_m_normflux, color="k", marker="D", fillstyle='none', markersize=0.5, label='o per (m95)')
    axes[1].plot(o_per_k_wavelengths, o_per_k_normflux, color="k", marker="D", fillstyle='none', markersize=0.5, label='o per (k94)')
    axes[0].plot(o_per_m_wavelengths, multi_fit_flux_m.best_fit, c = 'green')
    axes[1].plot(o_per_k_wavelengths,multi_fit_flux_k.best_fit, c= 'green')
    fig_1.legend()

    compare_data('omiper.m95.fit_data.txt', multi_fit_flux_m, b_o_per, N_o_per, v_rad_o_per)
    compare_data('ompier.k94.fit_data.txt', multi_fit_flux_k, b_o_per, N_o_per, v_rad_o_per)

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
    sigsco_wavelengths, sigsco_normflux = np.asarray(file_reader(files[2]))
    # fit a voigt profile to data using parameters defined above
    fit_flux_sigsco = vp.voigt_absorption_line(
        sigsco_wavelengths,
        lambda0=7698.974,
        b=b_sigsco,
        N=N_sigsco,
        f=3.393e-1,
        gamma=3.8e7,
        v_rad=v_rad_sigsco,
        v_resolution=1.20,
    )


    multi_fit_flux_sig= vp.fit_multi_voigt_absorptionlines(wavegrid= sigsco_wavelengths,
                                       ydata= sigsco_normflux,
                                       restwave= 7698.974,
                                       f= 3.393e-1,
                                       gamma= 3.8e7,
                                       b= b_sigsco,
                                       N= N_sigsco,
                                       v_rad= v_rad_sigsco,
                                       v_resolution= 1.2,
                                       n_step= 25)

    # calculate error between these wavelength which have been identified as the flat part of the spectrum by eye
    error_sigsco = errors(7699.0, 7699.4, sigsco_wavelengths, sigsco_normflux)
    # calculate reduced chi squared of fit
    print('sigsco reduced \u03C7\u00B2', reduced_chi_squared(sigsco_normflux, error_sigsco, fit_flux_sigsco))

    # plot data and corresponding fits on the same graph but on different subplots
    fig_2 = plt.figure()
    ax = fig_2.add_subplot(111)
    ax.plot(sigsco_wavelengths, sigsco_normflux, color="k", marker="D", fillstyle='none', markersize=0.5, label='sigsco')
    ax.plot(sigsco_wavelengths, fit_flux_sigsco, c='b', label='fit')
    ax.plot(sigsco_wavelengths, multi_fit_flux_sig.best_fit, c= 'green')
    # set x limit to allow for better analyse of data by eye
    ax.set_xlim(7698, 7699.5)
    ax.set_xlabel('Wavelength \u00C5')
    ax.set_ylabel('Normalised Flux')
    fig_2.legend()

    compare_data('sigsco.fit_data.txt', multi_fit_flux_sig, b_sigsco, N_sigsco, v_rad_sigsco)

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
    zetoph_k_wavelengths, zetoph_k_normflux = np.asarray(file_reader(files[3]))
    # fit a voigt profile to data using parameters defined above

    fit_flux_zetoph_k = vp.voigt_absorption_line(
        zetoph_k_wavelengths,
        lambda0=7698.974,
        b=b_zetoph,
        N=N_zetoph,
        f=3.393e-1,
        gamma=3.8e7,
        v_rad=v_rad_zetoph,
        v_resolution=1.40,
    )


    multi_fit_flux_k = vp.fit_multi_voigt_absorptionlines(wavegrid= zetoph_k_wavelengths,
                                       ydata= zetoph_k_normflux,
                                       restwave= 7698.974,
                                       f= 3.393e-1,
                                       gamma= 3.8e7,
                                       b= b_zetoph,
                                       N= N_zetoph,
                                       v_rad= v_rad_zetoph,
                                       v_resolution= 1.40,
                                       n_step= 25)

    # calculate error between these wavelength which have been identified as the flat part of the spectrum by eye
    error_zetoph_k = errors(7698.8, 7699, zetoph_k_wavelengths, zetoph_k_normflux)
    # calculate reduced chi squared of fit
    print('Zetoph_k reduced \u03C7\u00B2', reduced_chi_squared(zetoph_k_normflux, error_zetoph_k, fit_flux_zetoph_k))

    # repeat process for same star but with data from another survey
    zetoph_lf_wavelengths, zetoph_lf_normflux = np.asarray(file_reader(files[4]))
    fit_flux_zetoph_lf = vp.voigt_absorption_line(
        zetoph_lf_wavelengths,
        lambda0=7698.974,
        b=b_zetoph,
        N=N_zetoph,
        f=3.393e-1,
        gamma=3.8e7,
        v_rad=v_rad_zetoph,
        v_resolution=0.40,
    )

    multi_fit_flux_lf = vp.fit_multi_voigt_absorptionlines(wavegrid=zetoph_lf_wavelengths,
                                                          ydata=zetoph_lf_normflux,
                                                          restwave=7698.974,
                                                          f=3.393e-1,
                                                          gamma=3.8e7,
                                                          b=b_zetoph,
                                                          N=N_zetoph,
                                                          v_rad=v_rad_zetoph,
                                                          v_resolution=0.40,
                                                          n_step=25)

    error_zetoph_lf = errors(7698.8, 7699, zetoph_lf_wavelengths, zetoph_lf_normflux)
    print('Zetoph_lf reduced \u03C7\u00B2', reduced_chi_squared(zetoph_lf_normflux, error_zetoph_lf, fit_flux_zetoph_lf))

    # repeat process for same star but with data from another survey
    zetoph_m_wavelengths, zetoph_m_normflux = np.asarray(file_reader(files[5]))
    fit_flux_zetoph_m = vp.voigt_absorption_line(
        zetoph_m_wavelengths,
        lambda0=7698.974,
        b=b_zetoph,
        N=N_zetoph,
        f=3.393e-1,
        gamma=3.8e7,
        v_rad=v_rad_zetoph,
        v_resolution=0.56,
    )

    multi_fit_flux_m = vp.fit_multi_voigt_absorptionlines(wavegrid=zetoph_m_wavelengths,
                                                          ydata=zetoph_m_normflux,
                                                          restwave=7698.974,
                                                          f=3.393e-1,
                                                          gamma=3.8e7,
                                                          b=b_zetoph,
                                                          N=N_zetoph,
                                                          v_rad=v_rad_zetoph,
                                                          v_resolution=0.56,
                                                          n_step=25)

    error_zetoph_m = errors(7698.8, 7699, zetoph_m_wavelengths, zetoph_m_normflux)
    print('Zetoph_m reduced \u03C7\u00B2', reduced_chi_squared(zetoph_m_normflux, error_zetoph_m, fit_flux_zetoph_m))

    # plot data and corresponding fits on the same graph but on different subplots
    fig_3 = plt.figure()
    axes = fig_3.add_subplot(311), fig_3.add_subplot(312), fig_3.add_subplot(313)
    axes[0].plot(zetoph_m_wavelengths, zetoph_m_normflux, color="k", marker="D", fillstyle='none', markersize=0.5,
                 label='zetoph')
    axes[0].plot(zetoph_m_wavelengths, fit_flux_zetoph_m, c='b', label='fit')
    axes[0].plot(zetoph_m_wavelengths, multi_fit_flux_m.best_fit, c= 'green')

    axes[1].plot(zetoph_lf_wavelengths, zetoph_lf_normflux, color="k", marker="D", fillstyle='none', markersize=0.5,
                 label='zetoph_2')
    axes[1].plot(zetoph_lf_wavelengths, fit_flux_zetoph_lf, c='b', label='fit')
    axes[1].plot(zetoph_lf_wavelengths, multi_fit_flux_lf.best_fit, c= 'green')

    axes[2].plot(zetoph_k_wavelengths, zetoph_k_normflux, color="k", marker="D", fillstyle='none', markersize=0.5,
                 label='zetoph_3')
    axes[2].plot(zetoph_k_wavelengths, fit_flux_zetoph_k, c='b', label='fit')
    axes[2].plot(zetoph_k_wavelengths, multi_fit_flux_k.best_fit, c= 'green')
    # set x limit to allow for better analyse of data by eye
    for i in range(3):
        axes[i].set_xlim(7698.4, 7699.4)
        axes[i].set_xlabel('Wavelength \u00C5')
        axes[i].set_ylabel('Normalised Flux')
    fig_3.legend()

    compare_data('zetoph.m95.fit_data.txt', multi_fit_flux_m, b_zetoph, N_zetoph, v_rad_zetoph)
    compare_data('zetoph.k94.fit_data.txt', multi_fit_flux_k, b_zetoph, N_zetoph, v_rad_zetoph)
    compare_data('zetoph.lf.fit_data.txt', multi_fit_flux_lf, b_zetoph, N_zetoph, v_rad_zetoph)

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
    zetper_wavelengths, zetper_normflux = np.asarray(file_reader(files[6]))
    # fit a voigt profile to data using parameters defined above
    fit_flux_zetper = vp.voigt_absorption_line(
        zetper_wavelengths,
        lambda0=7698.974,
        b=b_zetper,
        N=N_zetper,
        f=3.393e-1,
        gamma=3.8e7,
        v_rad=v_rad_zetper,
        v_resolution=0.56,
    )

    multi_fit_flux_zet = vp.fit_multi_voigt_absorptionlines(wavegrid=zetper_wavelengths,
                                                          ydata=zetper_normflux,
                                                          restwave=7698.974,
                                                          f=3.393e-1,
                                                          gamma=3.8e7,
                                                          b=b_zetper,
                                                          N=N_zetper,
                                                          v_rad=v_rad_zetper,
                                                          v_resolution=0.56,
                                                          n_step=25)

    # calculate error between these wavelength which have been identified as the flat part of the spectrum by eye
    error_zetper = errors(7699.5, 7700, zetper_wavelengths, zetper_normflux)
    # calculate reduced chi squared of fit
    print('Zetper reduced \u03C7\u00B2', reduced_chi_squared(zetper_normflux, error_zetper, fit_flux_zetper))

    # plot data and corresponding fits on the same graph but on different subplots
    fig_4 = plt.figure()
    ax = fig_4.add_subplot(111)
    ax.plot(zetper_wavelengths, zetper_normflux, color="k", marker="D", fillstyle='none', markersize=0.5, label='zetper')
    ax.plot(zetper_wavelengths, fit_flux_zetper, c='b')
    #ax.plot(zetper_wavelengths, zetper_lmfit_flux, c='r')
    ax.plot(zetper_wavelengths, multi_fit_flux_zet.best_fit, c= 'green')
    # set x limit to allow for better analyse of data by eye
    ax.set_xlim(7698.5, 7700)
    ax.set_xlabel('Wavelength \u00C5')
    ax.set_ylabel('Normalised Flux')
    fig_4.legend()


    compare_data('zetper.fit_data.txt', multi_fit_flux_zet, b_zetper, N_zetper, v_rad_zetper)

folder = Path(PYTHONDIR + "/data")
filename = folder / "voigt_benchmarkdata" / 'parameter_modelling_data' / "files_for_parameter_modelling.txt"

    # state what headers the desired data is under
Headers = ["star_name", "file_name", "star_file", "resolution"]

        # read in the data
file = pd.read_csv(
    filename,
     delim_whitespace=True,
     header=None,
     names=Headers,
     engine="python",
    )

files = file['file_name']
resolution = file['resolution']
star_data = file["star_file"]
star_name = file["star_name"]

loop = True
if not(loop):
    # run analyse of the 4 stars
    o_per()
    sigsco()
    zetoph()
    zetper()

else:

    #for i in range(len(files)):
    for i in 0:

        wavelengths, normflux = np.asarray(file_reader(files[i]))

        folder = Path(PYTHONDIR + "/data")
        filename = folder / "voigt_benchmarkdata" / 'parameter_modelling_data' / star_data[i]


        Headers = ["b", "N", "v_rad"]
        star_parameters =pd.read_csv(filename,
         delim_whitespace = True,
         header=None,
         names=Headers,
         engine="python",
        )

        b = np.asarray(star_parameters["b"])
        N = np.asarray(star_parameters["N"])
        v_rad = np.asarray(star_parameters["v_rad"])

        WH_flux= vp.voigt_absorption_line(
            wavelengths,
            lambda0=7698.974,
            b=star_parameters["b"],
            N=star_parameters["N"],
            f=3.393e-1,
            gamma=3.8e7,
            v_rad=star_parameters["v_rad"],
            v_resolution=0.56,
        )

        continuum = np.array([])
        for j in range(len(normflux)):
            if normflux[j] > 0.98:
                continuum = np.append(continuum,normflux[j])
        error = np.array([np.std(continuum)] * len(normflux))
        print('reduced chi squared for {0} from the WH model is {1}'.format(star_name[i],(reduced_chi_squared(normflux, error, WH_flux))))
        fit_flux = vp.fit_multi_voigt_absorptionlines(wavegrid= wavelengths,
                                           ydata= normflux,
                                           restwave= 7698.974,
                                           f= 3.393e-1,
                                           gamma= 3.8e7,
                                           b= star_parameters['b'],
                                           N= star_parameters['N'],
                                           v_rad= star_parameters['v_rad'],
                                           v_resolution= resolution[i],
                                           n_step= 25)


        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(wavelengths, normflux, color="k", marker="D", fillstyle='none', label= 'Data')
        ax.plot(wavelengths, fit_flux.best_fit, c='r', marker = '*', label= 'best-fit')
        ax.plot(wavelengths, WH_flux, c='g', label='Wetly & Hobbs')
        fig.suptitle(star_name[i])
        ax.set_xlabel('Wavelength \u00C5')
        ax.set_ylabel('Normalised Flux')
        plt.legend()

        save_path = 'c:/Users/user/edibles/edibles/data/voigt_benchmarkdata/parameter_modelling_data/'
        file_name = star_name[i] + '.fit_data.txt'
        file_path = os.path.join(save_path, file_name)
        file = open(file_path, 'w+')

        txt = 'Welty & Hobbs parameteres: \n b: {0}, \n N: {1}, \n v_rad: {2} \n' \
              'Fit parameters: {3}'.format(b, N, v_rad, fit_flux.fit_report())
        file.write(txt)
        file.close()

print('done')
plt.show()
