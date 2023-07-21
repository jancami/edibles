import numpy as np
from scipy.special import wofz
from scipy.interpolate import interp1d
import astropy.constants as cst
import matplotlib.pyplot as plt
from edibles import PYTHONDIR
from edibles.utils.edibles_oracle import EdiblesOracle
from edibles.utils.edibles_spectrum import EdiblesSpectrum
from pathlib import Path
import pandas as pd
from scipy.ndimage import gaussian_filter
from lmfit import Parameters, minimize,Model

from scipy.optimize import fmin


def voigt_profile(x, sigma, gamma):
    """
    Function to return the value of a (normalized) Voigt profile centered at x=0
    and with (Gaussian) width sigma and Lorentz damping (=HWHM) gamma.

    The Voigt profile is computed using the scipy.special.wofz, which returns
    the value of the Faddeeva function.


    WARNING
    scipy.special.wofz is not compaible with np.float128 type parameters.

    Args:
        x (float64): Scalar or array of x-values
        sigma (float64): Gaussian sigma component
        gamma (float64): Lorentzian gamma (=HWHM) component

    Returns:
        ndarray: Flux array for given input

    """

    z = (x + 1j * gamma) / sigma / np.sqrt(2)

    return np.real(wofz(z)) / sigma / np.sqrt(2 * np.pi)

def voigt_optical_depth(wave, lambda0=0.0, b=0.0, N=0.0, f=0.0, gamma=0.0, v_rad=0.0):
    """
    Function to return the value of a Voigt optical depth profile at a given wavelength, for a line
    centered at lambda0.

    Args:
        wave (float64): Wavelength at which optical depth is to be calculatd (in Angstrom)
        lambda0 (float64): Central (rest) wavelength for the absorption line, in Angstrom.
        b (float64): The b parameter (Gaussian width), in km/s.
        N (float64): The column density (in cm^{-2})
        f (float64): The oscillator strength (dimensionless)
        gamma (float64): Lorentzian gamma (=HWHM) component
        v_rad (float64): Radial velocity of absorption line (in km/s)

    Returns:
        float64: Optical Depth at wave.

    """

    # All we have to do is proper conversions so that we feed the right numbers into the call
    # to the VoigtProfile -- see documentation for details.
    nu = cst.c.to("angstrom/s").value / wave
    nu0 = cst.c.to("angstrom/s").value / lambda0
    sigma = (b * 1e13) / lambda0 / np.sqrt(2)
    gamma_voigt = gamma / 4 / np.pi
    tau_factor = (N * np.pi * cst.e.esu ** 2 / cst.m_e.cgs / cst.c.cgs * f).value

    # print("Nu0 is:        " + "{:e}".format(nu0))
    # print("Sigma is:      " + "{:e}".format(sigma))
    # print("Gamma is:      " + "{:e}".format(gamma_voigt))
    # print("Tau_factor is: " + "{:e}".format(tau_factor))

    # Transform this into a frequency grid centered around nu0

    ThisVoigtProfile = voigt_profile(nu - nu0, sigma, gamma_voigt)
    tau = tau_factor * ThisVoigtProfile

    return tau

# Use this method to reproduce data from overleaf documents


def voigt_absorption_line(
    wavegrid, lambda0=0.0, f=0.0, gamma=0.0, b=0.0, N=0.0, v_rad=0.0, v_resolution=0.0, n_step=25, debug=False
):
    """
    Function to return a complete Voigt Absorption Line Model, smoothed to the specified
    resolution and resampled to the desired wavelength grid.
    This can in fact be a set of different absorption lines -- same line, different
    cloud components or different line for single cloud component.

    Args:
        wavegrid (float64): Wavelength grid (in Angstrom) on which the final result is desired.
        lambda0 (float64): Central (rest) wavelength for the absorption line, in Angstrom.
        b (float64): The b parameter (Gaussian width), in km/s.
        N (float64): The column density (in cm^{-2})
        f (float64): The oscillator strength (dimensionless)
        gamma (float64): Lorentzian gamma (=HWFM) component
        v_rad (float64): Radial velocity of absorption line (in km/s)
        v_resolution (float64): Instrument resolution in velocity space (in km/s)
        n_step (int): no. of point per FWHM length, governing sampling rate and efficiency
        debug (bool): If True, info on the calculation will be displayed

    Returns:
        ndarray: Normalized flux for specified grid & parameters.

    """
    # Let's convert the parameters to numpy arrays first.
    lambda0_array = np.array(lambda0, ndmin=1)
    f_array = np.array(f, ndmin=1)
    gamma_array = np.array(gamma, ndmin=1)
    b_array = np.array(b, ndmin=1)
    N_array = np.array(N, ndmin=1)
    v_rad_array = np.array(v_rad, ndmin=1)

    # How many lines are passed on?
    n_lines = lambda0_array.size
    if debug:
        print("Number of lines: " + "{:d}".format(n_lines))

    # How many cloud components do we have?
    n_components = N_array.size
    if debug:
        print("Number of lines: " + "{:d}".format(n_lines))
        print("Number of components: " + "{:d}".format(n_components))

    # We will consider 3 different cases here:
    # 1. A single line, but multiple components.
    #    --> Each component represents a cloud; use the same spectral line parameters for each
    #        component. We will set that up here, and then do a recursive call to
    #        voigt_absorption_line.
    # 2. Multiple lines, but a single component.
    #    --> E.g. the Na doublet lines. For each line component, we will use the same
    #        cloud components. We will set that up here, and then do a recursive call
    #        to voigt_absorption_line.
    # 3. Multiple lines, and multiple components.
    #    Here, we consider 2 different cases:
    #    A) If n_lines == n_components, we will treat each combination as a unique, single line.
    #       We could have gotten here from the recursive calls in 1. and 2. Exception: if the
    #       keyword /MultiCloud is set. In that case we will proceed as in B).
    #    B) If n_lines <> n_components, we will interpret this as meaning that each
    #       component will produce each of the lines. This means we now have to create
    #       n_lines * n_components entries for the recursive call. We will set that up here.

    if (n_lines == 1) & (n_components != 1):
        # Case 1 from above. Replicate all the line parameters for each of the components,
        # and issue the call to self.
        lambda0_use = np.repeat(lambda0, n_components)
        f_use = np.repeat(f, n_components)
        gamma_use = np.repeat(gamma, n_components)
        # we should really check whether we have the correct number of other parameters....
        interpolated_model = voigt_absorption_line(
            wavegrid,
            lambda0=lambda0_use,
            f=f_use,
            gamma=gamma_use,
            b=b,
            N=N,
            v_rad=v_rad,
            v_resolution=v_resolution,
            n_step=n_step
        )
    elif (n_components == 1) & (n_lines != 1):
        # Case 2 from above. Replicate all the sightline parameters for each of the lines,
        # and issue the call to self.
        b_use = np.repeat(b, n_lines)
        N_use = np.repeat(N, n_lines)
        v_rad_use = np.repeat(v_rad, n_lines)
        interpolated_model = voigt_absorption_line(
            wavegrid,
            lambda0=lambda0,
            f=f,
            gamma=gamma,
            b=b_use,
            N=N_use,
            v_rad=v_rad_use,
            v_resolution=v_resolution,
            n_step=n_step
        )
    elif n_lines == n_components:
        # Case 3A from above.
        if debug:
            print("Number of components: ", n_lines)
        # We can process each line/component now with its own set of parameters.
        # We will loop over each line, create the proper wavelength grid, then get the
        # corresponding optical depth profile, and then decide how to combine everything.
        #
        # One thing to ensure though is that the step size in velocity space is the same for
        # each component -- otherwise the Gaussian smoothing at the end will go wrong.
        #
        # And we want the v_Grid is sufficiently finely sampled, so that:
        # 1. There are at least 7 data points within each FWHM, i.e. a oversample ratio of 3??
        # 2. dv is no larger than the step of input x-grid.

        Voigt_FWHM = VoigtFWHM(lambda0_array, gamma_array, b_array)
        FWHM2use = np.min(np.append(Voigt_FWHM, v_resolution))
        xgrid_test = np.asarray(wavegrid)
        dv_xgrid = np.median(xgrid_test[1:] - xgrid_test[0:-1]) / \
            np.mean(xgrid_test) * cst.c.to("km/s").value
        n_step_dv = np.ceil(FWHM2use / dv_xgrid)

        if n_step < np.max([7, n_step_dv]):
            n_step = np.max([7, n_step_dv])
            print("n_step too small. To avoid under-sampling, n_step reset to %d" % (n_step))
        v_stepsize = FWHM2use / n_step

        # We will also need to be able to add each optical depth profile, so we need a common
        # wavelength grid to interpolate on.
        # We use pm 8.5 * FWHM for each line, corresponding to pm 20*b assuming pure Gaussian,
        # and see what wavelength limits to consider.
        bluewaves = lambda0_array * (
            1.0 + (v_rad_array - 8.5 * Voigt_FWHM) / cst.c.to("km/s").value
        )
        redwaves = lambda0_array * (
            1.0 + (v_rad_array + 8.5 * Voigt_FWHM) / cst.c.to("km/s").value
        )
        #print("Bluewaves:", bluewaves)
        #print("Waves    :", lambda0_array)
        #print("Redwaves :", redwaves)
        #print("b_array  :", b_array)
        minwave = bluewaves.min()
        maxwave = redwaves.max()
        minwave = min(minwave, wavegrid.min())
        maxwave = max(maxwave, wavegrid.max())

        # print(v_rad_array)
        #print("Wave range: ", minwave, maxwave)
        n_v = int(
            np.ceil((maxwave - minwave) / minwave * cst.c.to("km/s").value / v_stepsize)
        )
        refgrid = minwave * (1.0 + np.arange(n_v) * v_stepsize / cst.c.to("km/s").value)
        allcomponents = np.zeros(shape=(n_v, n_lines))
        for lineloop in range(n_lines):
            dv = getVGrid(
                lambda0_array[lineloop],
                gamma_array[lineloop],
                b_array[lineloop],
                v_resolution,
                n_step)
            thiswavegrid = lambda0_array[lineloop] * (1.0 + dv / cst.c.to("km/s").value)
            tau = voigt_optical_depth(
                thiswavegrid,
                lambda0=lambda0_array[lineloop],
                b=b_array[lineloop],
                N=N_array[lineloop],
                f=f_array[lineloop],
                gamma=gamma_array[lineloop],
                v_rad=v_rad_array[lineloop],
            )
            if debug:
                print("Max tau:", tau.max())
            # Shift to the proper wavelength given the radial velocity
            vel = dv + v_rad_array[lineloop]
            thiswavegrid = lambda0_array[lineloop] * (
                1.0 + vel / cst.c.to("km/s").value
            )
            # Interpolate to reference grid
            interpolationfunction = interp1d(
                thiswavegrid, tau, kind="cubic", fill_value="extrapolate"
            )
            tau_grid = interpolationfunction(refgrid)
            tau_grid[np.where(refgrid > np.max(thiswavegrid))] = 0
            tau_grid[np.where(refgrid < np.min(thiswavegrid))] = 0
            # plt.plot(thiswavegrid,tau,marker="+")
            #plt.plot(refgrid,tau_grid, color='red')
            # plt.show()

            allcomponents[:, lineloop] = tau_grid

        # Now add up all the optical depth components.
        tau = np.sum(allcomponents, axis=1)

        # plt.plot(refgrid,tau)
        # plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
        # plt.show()

        # Do the radiative transfer
        AbsorptionLine = np.exp(-tau)
        # plt.plot(refgrid,AbsorptionLine)
        # plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
        # plt.show()

        # Apply a Gaussian instrumental smoothing function!
        # Calculate sigma -- in units of step size!
        smooth_sigma = fwhm2sigma(v_resolution) / v_stepsize
        if debug:
            print("Smoothing sigma is: " + "{:e}".format(smooth_sigma))

        # One thing to watch out for is that the smoothing width is large compared to

        gauss_smooth = gaussian_filter(AbsorptionLine, sigma=smooth_sigma)
        interpolationfunction = interp1d(
            refgrid, gauss_smooth, kind="cubic", bounds_error=False, fill_value=(1, 1)
        )
        interpolated_model = interpolationfunction(wavegrid)
        #plt.plot(refgrid, AbsorptionLine, marker='o')
        #plt.plot(refgrid, gauss_smooth, color='green', marker='D')
        #plt.plot(wavegrid,interpolated_model, color='red', marker='1')
        # plt.show()
    else:
        # Create arrays that hold all the lines for all the components.
        # We need in total n_components * n_lines array elements.
        lambda0_use = np.repeat(lambda0, n_components)
        f_use = np.repeat(f, n_components)
        gamma_use = np.repeat(gamma, n_components)
        # For the eomponents, do some dimensional juggling....
        b_use = np.concatenate(np.repeat([b], n_lines, axis=0), axis=0)
        N_use = np.concatenate(np.repeat([N], n_lines, axis=0), axis=0)
        v_rad_use = np.concatenate(np.repeat([v_rad], n_lines, axis=0), axis=0)
        # print(lambda0_use)
        # print(N_use)
        interpolated_model = voigt_absorption_line(
            wavegrid,
            lambda0=lambda0_use,
            f=f_use,
            gamma=gamma_use,
            b=b_use,
            N=N_use,
            v_rad=v_rad_use,
            v_resolution=v_resolution,
            n_step=n_step
        )
        #"voigt_absorption_line Panic: This option has not been implemented yet.... "
        # )

    return interpolated_model


def fwhm2sigma(fwhm):
    """
    Simple function to convert a Gaussian FWHM to Gaussian sigma.
    """
    sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    return sigma


def VoigtFWHM(lambda0, gamma, b):
    """
    Calculate FWHM of Voigt using the formula by Olivero.
    See the bottom of wiki Voigt Profile page

    :param lambda0: center wavelength in AA
    :param gamma: Lorentzian gamma (=HWHM), in frequency
    :param b: Gaussian b in km/s.
    :return: Gau_FWHM, Lor_FWHM, V_FWHM, all in km/s
    """

    #gamma_AA = gamma * lambda0**2 / cst.c.to("angstrom/s").value
    gamma_kms = gamma * lambda0 * 1e-13
    Lor_FWHM = gamma_kms * 2

    # FWHM = 2*sigma*sqrt(2*ln2)
    Gau_FWHM = 2.35482 * b

    #fV = fL/2 + sqrt(fL^2/4 + fG^2)
    Voigt_FWHM = 0.5 * Lor_FWHM + np.sqrt(Lor_FWHM**2 / 4 + Gau_FWHM ** 2)

    return Voigt_FWHM


def getVGrid(lambda0, gamma, b, v_resolution, n_step):
    """
    Calculate v grid to be used for Voigt profile
    Size of grid is determined by the greater of FWHM_Voigt and v_resolution

    :return: dv, velocity grid
    """

    Voigt_FWHM = VoigtFWHM(lambda0, gamma, b)
    FWHM2use = np.max([Voigt_FWHM, v_resolution])
    v_stepsize = FWHM2use / n_step
    dv = np.arange(
        start=-8.5*FWHM2use, stop=8.5*FWHM2use, step=v_stepsize
    )
    return dv
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
    AbsorptionLine = voigt_absorption_line(
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

def wavegrid_lmfit(params):
    """
    currently being used to test fitting process for attempting to obtain parameers from absoption curve

    :param params:
    :param data:
    :return:

    """
    wavegrid_array = np.arange(1000) * 0.001 + 7698.974
    b_array = params['b']
    N_array = params['N']
    v_rad_array = params['v_rad']
    v_resolution_array = params['v_res']
    AbsorptionLine = voigt_absorption_line(
        wavegrid_array,
        lambda0=7698.974,
        b=b_array,
        N=N_array,
        f=3.393e-1,
        gamma=3.8e7,
        v_rad=v_rad_array,
        v_resolution=v_resolution_array,
    )
    return AbsorptionLine

def wavegrid_lmfit(params):
    """
    currently being used to test fitting process for attempting to obtain parameers from absoption curve

    :param params:
    :param data:
    :return:

    """
    wavegrid_array = np.arange(1000) * 0.001 + 7698.974
    b_array = params['b']
    N_array = params['N']
    v_rad_array = params['v_rad']
    AbsorptionLine = voigt_absorption_line(
        wavegrid_array,
        lambda0=7698.974,
        b=b_array,
        N=N_array,
        f=3.393e-1,
        gamma=3.8e7,
        v_rad=v_rad_array,
        v_resolution=0.56,
    )
    return AbsorptionLine


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
    parameters.add('b', value = 0.6)
    parameters.add('N', value = 12.5e10)
    parameters.add('v_rad', value = 10.5)
    #parameters.add('v_res', value=0.56)
    #parameters.add('f', value=3.393e-1)
    #parameters.add('gamma', value=3.8e7)
    #parameters.add('lambda0', value=7698.974)

    #parameters.add('v_res', min=0.3, max=1.5)

    model = Model(voigt_absorption_line)
    fit = model.fit(o_per_m[1], Parameters, v_res= 0.56, f=3.393e-1, gamma=3.8e7, lambda0= 7698.974)
    #wavegrid, lambda0=0.0, f=0.0, gamma=0.0, b=0.0, N=0.0, v_rad=0.0, v_resolution=0.0
    print(fit)
    #print(type(parameters))

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


#  --------------------------------------------------------------------------------------


if __name__ == "__main__":

    from math import *

    """
    This is a series of examples / tests to see if the implementation of the various Voigt
    functions has been done correctly. It includes tests for the basic Voigt function, as
    well as for the various forms of the normalized Voigt profiles. One test aims to
    reproduce the high-resolution profile for the K line of omi Per (form Welty et al.)
    """
    show_example = 5

    if show_example == 1:
        #############################################################
        #
        # EXAMPLE 1: a single line.
        #
        #############################################################
        wavegrid = np.arange(1000) * 0.01 + 5000
        AbsorptionLine = voigt_absorption_line(
            wavegrid,
            lambda0=5003.0,
            b=0.75,
            N=1e11,
            f=1,
            gamma=1e7,
            v_rad=-10.0,
            v_resolution=0.26,
        )
        plt.plot(wavegrid, AbsorptionLine, marker="1")
        plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
        plt.show()

    elif show_example == 2:
        #############################################################
        #
        # EXAMPLE 2: Single line, multiple cloud components.
        #
        #############################################################
        # This, in fact, is reproducing Dan Welty's high-resolution
        # observations of the K line of o Per.
        folder = Path(PYTHONDIR + "/data")
        filename = folder / "voigt_benchmarkdata" / "omiper.m95.7698.txt"
        Headers = ["Wavelength", "Normfluxval"]
        omiper_data = pd.read_csv(
            filename,
            delim_whitespace=True,
            skiprows=[0],
            header=None,
            names=Headers,
            engine="python",
        )

        lambda0 = 7698.974  # K wavelength in angstrom.
        f = 3.393e-1  # oscillator strength
        gamma = 3.8e7  # Gamma for K line
        b = [0.60, 0.44, 0.72, 0.62, 0.60]  # b parameter in km/s
        N = np.array([12.5, 10.0, 44.3, 22.5, 3.9]) * 1e10  # Column density
        v_rad = (
            np.array([10.50, 11.52, 13.45, 14.74, 15.72]) + 0.1
        )  # Radial velocity of cloud in km/s
        v_resolution = 0.56  # Instrumental resolution in km / s

        AbsorptionLine = voigt_absorption_line(
            omiper_data["Wavelength"],
            lambda0=lambda0,
            b=b,
            N=N,
            f=f,
            gamma=gamma,
            v_rad=v_rad,
            v_resolution=v_resolution,
        )

        plt.plot(
            omiper_data["Wavelength"],
            omiper_data["Normfluxval"],
            marker="D",
            fillstyle="none",
            color="black",
        )
        plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
        plt.xlim(7699.15, 7699.45)
        plt.plot(omiper_data["Wavelength"], AbsorptionLine, color="orange", marker="+")
        plt.show()

    elif show_example == 3:
        #############################################################
        #
        # EXAMPLE 3: Na doublet: multiple lines, single cloud.
        #
        #############################################################

        # Let's see if we can reproduce the Na lines for HD 145502
        lambda0 = [3302.369, 3302.978]
        f = [8.26e-03, 4.06e-03]
        gamma = [6.280e7, 6.280e7]
        b = [0.8]
        N = [1.8e13]
        v_rad = 18.9
        v_resolution = 5.75

        pythia = EdiblesOracle()
        List = pythia.getFilteredObsList(
            object=["HD 145502"], MergedOnly=True, Wave=3302.0
        )
        test = List.values.tolist()
        filename = test[0]
        print(filename)
        wrange = [3301.5, 3304.0]
        sp = EdiblesSpectrum(filename)
        wave = sp.wave
        flux = sp.flux
        idx = np.where((wave > wrange[0]) & (wave < wrange[1]))
        wave = wave[idx]
        flux = flux[idx]
        flux = flux / np.median(flux)
        AbsorptionLine = voigt_absorption_line(
            wave,
            lambda0=lambda0,
            b=b,
            N=N,
            f=f,
            gamma=gamma,
            v_rad=v_rad,
            v_resolution=v_resolution,
        )
        plt.plot(wave, flux)
        plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
        plt.plot(wave, AbsorptionLine, color="orange", marker="*")
        plt.show()

    elif show_example == 4:
        #############################################################
        #
        # EXAMPLE 4: Na doublet: multiple lines, multiple cloud.
        #
        #############################################################

        # Let's see if we can reproduce the Na lines for HD 183143
        lambda0 = [3302.369, 3302.978]
        f = [8.26e-03, 4.06e-03]
        gamma = [6.280e7, 6.280e7]
        b = [1.0, 1.4, 1.4]
        N = [1e13, 1.5e14, 5.0e14]
        v_rad = [1.0, 8.0, 22.0]
        v_resolution = 5.75

        pythia = EdiblesOracle()
        List = pythia.getFilteredObsList(object=["HD 183143"], MergedOnly=True, Wave=3302.0)
        test = List.values.tolist()
        filename = test[0]
        print(filename)
        wrange = [3301.5, 3304.0]
        sp = EdiblesSpectrum(filename)
        wave = sp.wave
        flux = sp.flux
        idx = np.where((wave > wrange[0]) & (wave < wrange[1]))
        wave = wave[idx]
        flux = flux[idx]
        flux = flux / np.median(flux)
        AbsorptionLine = voigt_absorption_line(
            wave,
            lambda0=lambda0,
            b=b,
            N=N,
            f=f,
            gamma=gamma,
            v_rad=v_rad,
            v_resolution=v_resolution,
        )
        plt.plot(wave, flux)
        plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
        plt.plot(wave, AbsorptionLine, color="orange", marker="*")
        plt.show()
    elif show_example == 5:
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

#calculate chi squares and bests parameters to reduce chi squared
#read up on lmfit attempt fitting routine