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


def multi_voigt_absorption_line( **params_list):
    """
    This function is essentially a wrapper around voigt_absorption_line in voigt_profile.py, that exists
    to allow it to create an instance of the Model class and create models with multiple Voigt commponents
    (both transitions and clouds). 
    This function will parse the list of parameters, and reformat them to call the voigt_absorption_line function. 
    It will then call that function, and return the resulting model. 

    Args:
        wavegrid: the wavelength grid on which to calculate the models. 
        n_trans: the number of unique transitions to calculate. 
                 If n_trans == 1, just one rest wavelength, f and gamma value is passed on. 
                 If n_trans == 2, we have a doublet and so on. 
        lambda0, lambda1, ...: rest wavelengths for each transition. 
        f0, f1, ...:           oscillator strengths for each transition. 
        gamma0, gamma1, ...:   Lorentz broadening parameter for each transition. 
        n_components: the number of different cloud components for which to calculate the profiles. 
        b0, b1, ....:  Doppler b-values for each cloud component. 
        N0, N1, ....:  column densities for each cloud component. 
        v_rad0, v_rad1, ...:  radial velocities for each cloud component/transition. 
        v_resolution: the velocity resolution (in km/s) of the desired final result. 
        n_step: the number of steps to sample the Voigt profile (default: 25). 
        debug:   Boolean: print debug info while running or not. 
    Returns:
    np.array
        Model array with normalized Voigt profiles corresponding to the specified parameters. 
    """

    # We should probably do parameter checking and set some defaults when parameters are missing, or 
    # issue an error or warning message. 
    #print('Entering multi_voigt_absorption_line....')
    n_trans = params_list['n_trans']
    wavegrid = params_list['wavegrid']
    n_components = params_list['n_components']
    v_resolution = params_list['v_resolution']
    n_step = params_list['n_step']
    
    # Initialize the lambda,f,gamma arrays for voigt_absorption_line
    all_lambda = np.empty(n_trans)
    all_f = np.empty(n_trans)
    all_gamma = np.empty(n_trans)
    #print(all_lambdas)
    #print(all_f)
    #print(all_gamma)
    #print('----')
    for i in range(n_trans):
        all_lambda[i] = params_list[f'lambda{i}']
        all_f[i] = params_list[f'f{i}']
        all_gamma[i] = params_list[f'gamma{i}']
    #print("Lambdas: ",all_lambda)
    #print("f:       ",all_f)
    #print("Gamma :  ",all_gamma)    
    # Initialize the b,N,v_rad arrays for voigt_absorption_line
    all_b = np.empty(n_components)
    all_N = np.empty(n_components)
    all_v_rad = np.empty(n_components)
    #print(all_b)
    #print(all_N)
    #print(all_v_rad)
    #print('----')
    for i in range(n_components):
        all_b[i] = params_list[f'b{i}']
        all_N[i] = params_list[f'N{i}']
        all_v_rad[i] = params_list[f'v_rad{i}']
    #print('After: ')
    #print("b:     ",all_b)
    #print("N:     ",all_N)
    #print("v_rad: ",all_v_rad)
    #print('----')

    # Now call voigt_absorption_line with these parameters.... 

    model = voigt_absorption_line(wavegrid, lambda0=all_lambda, f=all_f, gamma=all_gamma, b=all_b, N=all_N, v_rad=all_v_rad, 
                                  v_resolution=v_resolution, n_step=n_step, debug=False)
    
    return model


def fit_multi_voigt_absorptionlines(wavegrid=np.array, ydata=np.array, restwave=np.array, f=np.array, gamma=np.array, 
             b=np.array, N=np.array, v_rad=np.array, v_resolution=0., n_step=0):
    """
    This function will take an observed spectrum contained in (wavegrid, ydata) and fit a set of Voigt profiles to
    it. The transitions to consider are specified by restwave, f, and gamma, and can be single floats or numpy arrays 
    containing the information for all transitions that we want to consider. 
    There can then be 1 or more clound components to consider in the sightline. The components are specified by b, N and v_rad
    that can be again floats or numpy arrays. 

    This function essentially creates a Model instance from the multi_voigt_absorption_line function, and most of the work done
    here is to "translate" the parameters from arrays to unique parameters (using the Parameters class) that will then be parsed 
    by multi_voigt_absorption_line.
    """
    
    # We should probably do lots of parameter checking first!!! To be done later.... 


    # How many transitions do we have? 
    restwave = np.asarray(restwave)
    n_trans = restwave.size
    #print("Fit_test: ", n_trans, " transitions....")

    # And how many cloud components? 
    b = np.asarray(b)
    n_components = b.size

    # Now create the Parameters class. 
    params = Parameters()
    # Create the parameters for the transitions. Those should *not* be free parameters for lmfit. 
    if n_trans == 1:
        params.add('lambda0', value=restwave, vary=False)
        params.add('f0', value=f, vary=False)
        params.add('gamma0', value=gamma, vary=False)
    else:
        #print('Multiple transitions found....')
        for i in range(n_trans):
            #print(i)
            #print(restwave[i])
            # The statement below creates variables lambda0, lambda1, lambda2, .... and similarly f0, f1, .. and gamma0, gamma1, ... 
            params.add(f'lambda{i}', value=restwave[i], vary=False)
            params.add(f'f{i}', value=f[i], vary=False) 
            params.add(f'gamma{i}', value=gamma[i], vary=False) 

    # Also create parameters for the other keywords -- these too should *not* be free parameters. 
    params.add('v_resolution', value=v_resolution, vary=False)
    params.add('n_step', value=n_step, vary=False)
    params.add('n_trans', value=n_trans, vary=False)
    params.add('n_components', value=n_components, vary=False)
   
    # Create the parameters for the clouds. 
    if n_components == 1:
        #print('Single Cloud Component:')
        params.add('b0', value=b,min=0)
        params.add('N0', value=N,min=0)
        params.add('v_rad0', value=v_rad)
    else:
        #print('Multiple Clouds:')
        for i in range(n_components):
            #print(i)
            #print(b[i])
            # Same idea as above -- will create b0, b1, ... and N0, N1, ... and v_rad0, v_rad1, .... 
            params.add(f'b{i}', value=b[i],min=0)
            params.add(f'N{i}', value=N[i],min=0)
            params.add(f'v_rad{i}', value=v_rad[i])
    #print(params)

    # Now create the Model instance. 
    voigtmod = Model(multi_voigt_absorption_line, independent_vars=['wavegrid'])
    #print("Resolution: ", v_resolution)

    # and do the fitting with the parameters we have created. 
    result=voigtmod.fit(ydata, params, wavegrid=wavegrid)

    return result



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

    elif show_example == 1:
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
