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
from edibles.utils import voigt_profile as vp

def read_benchmark_data (filename):
    """
    Will read in the specified filename from the voigt_benchmark folder.
    Is used in order to read in data for o per, sigma sco, zeta per, and zeta oph stars in example 5
    :paramaters :
        filename: string
            the file name containing the data that needs to be analysed
    :return:
        data["Wavelength"] : 1d ndarray
            Contains all the measured wavelength values for a particular file
        data["Normfluxval"] : 1d ndarray
            Contains normailised flux values corresponding to each wavelength reading
    """
    # define where the desired file is
    folder = Path(PYTHONDIR + "/data")
    full_filename = folder / "voigt_benchmarkdata" / filename

    # state what headers the desired data is under
    Headers = ["Wavelength", "Normfluxval"]

    # read in the data
    data = pd.read_csv(
        full_filename,
        delim_whitespace=True,
        skiprows=[0],
        header=None,
        names=Headers,
        engine="python",
    )
    return data["Wavelength"], data["Normfluxval"]

def multi_voigt_absorption_line( **params_list):
    """
    This function will calculate the resulting Voigt absorption line profiles for multiple lines. 
    The parameter list will be parsed for all the parameters -- see below. 
    Essentially, we will reformat the parameters to then issue the proper call to
    voigt_absorption_line. 

    Args:
        wavegrid: the wavelength grid on which to calculate the models. 
        n_trans: the number of unique transitions to calculate. The default is 1, meaning that just
                 one set of lambda0 and f values are passed on. 
                 n_trans = 2 for a doublet and so on. 
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

    print('Entering multi_voigt_absorption_line....')
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
    print("Lambdas: ",all_lambda)
    print("f:       ",all_f)
    print("Gamma :  ",all_gamma)    
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
    print("b:     ",all_b)
    print("N:     ",all_N)
    print("v_rad: ",all_v_rad)
    print('----')

    model = vp.voigt_absorption_line(wavegrid, lambda0=all_lambda, f=all_f, gamma=all_gamma, b=all_b, N=all_N, v_rad=all_v_rad, 
                                  v_resolution=v_resolution, n_step=n_step, debug=False)
    
    return model

def fit_test(wavegrid=np.array, ydata=np.array, restwave=np.array, f=np.array, gamma=np.array, 
             b=np.array, N=np.array, v_rad=np.array, v_resolution=0.56, n_step=25):
    """
    Fit a multi-voigt profile for the K line at 7698AA to the data contained in ydata
    Lots of things hardcoded right now -- just to get it going.... 
    testfit = fit_test(wavegrid=wave, ydata=normflux, restwave=7698.974, f=3.393e-1, gamma=3.8e7, 
                       b=b_o_per, N=N_o_per, v_rad=v_rad_o_per, v_resolution=0.56, n_step=25)
    """
    
    # How many transitions do we have? 
    restwave = np.asarray(restwave)
    n_trans = restwave.size

    # And how many cloud components? 
    b = np.asarray(b)
    n_components = b.size

    voigtmod = Model(multi_voigt_absorption_line, independent_vars=['wavegrid'])
    params = Parameters()

    # Create the parameters for the transitions. Those should not be free parameters. 
    if n_trans == 1:
        params.add('lambda0', value=restwave, vary=False)
        params.add('f0', value=f, vary=False)
        params.add('gamma0', value=gamma, vary=False)
    else:
        for i in range(n_trans):
            print(i)
            print(restwave[i])
        params.add(f'lambda{i}', value=restwave[i])


    # Create the parameters for the clouds. 
    if n_components == 1:
        print('Single Cloud Component:')
        params.add('b0', value=b)
        params.add('N0', value=N)
        params.add('_vrad0', value=v_rad)
    else:
        print('Multiple Clouds:')
        for i in range(n_components):
            print(i)
            print(b[i])
            params.add(f'b{i}', value=b[i])
            params.add(f'N{i}', value=N[i])
            params.add(f'v_rad{i}', value=v_rad[i])
    #print(params)
    #print(v_resolution,n_step)
    result=voigtmod.fit(ydata, params, wavegrid=wavegrid, v_resolution=v_resolution, n_step=n_step, n_trans=n_trans, 
                        n_components=n_components)
    return result


if __name__ == "__main__":
    # Let's try to get things to work for o Per first. This a single line / multiple component case.  
    # These should be the parameters (from Welty et al.)
    b_o_per = [0.60, 0.44, 0.72, 0.62, 0.60]
    N_o_per = [12.5e10, 10e10, 44.3e10, 22.5e10, 3.9e10]
    v_rad_o_per = [10.5+0.15, 11.52+0.15, 13.45+0.15, 14.74+0.15, 15.72+0.15]

    # array of files with the data that will be studied in the example
    files = ["omiper.m95.7698.txt", "omiper.k94.7698.txt", "sigsco.k00a.7698.txt", "zetoph.k94.7698.txt",
             "zetoph.lf.7698.txt", "zetoph.m95.7698.txt", "zetper.m95.7698.txt"]
    wave,normflux = read_benchmark_data(files[0])
    # Let's first get the combined Voigt profiles using all absorption components from Welty. 
    welty_fit=vp.voigt_absorption_line(wave,lambda0=7698.974,f=3.393e-1,gamma=3.8e7,b=b_o_per,N=N_o_per,v_rad=v_rad_o_per,v_resolution=0.56,n_step=25)

    # Then let's see if we can fit a single Voigt profile using the LMFIT model class. 
    model=Model(vp.voigt_absorption_line, independent_vars=['wavegrid'])
    model.set_param_hint('lambda0', value=7698.974, vary=False)
    model.set_param_hint('f', value=3.393e-1, vary=False)
    model.set_param_hint('gamma', value=3.8e7, vary=False)
    model.set_param_hint('v_resolution', value=0.56, vary=False)
    model.set_param_hint('n_step', value=25, vary=False)
    model.set_param_hint('N', value=1.25e11)
    model.set_param_hint('b', value=0.60)
    model.set_param_hint('v_rad', value=12.5)
    params=model.make_params()
    params.pretty_print()
    print(' ')
    result = model.fit(normflux,params,wavegrid=wave)
    result.params.pretty_print()

    
    testmodel = multi_voigt_absorption_line(wavegrid=wave, n_trans=1, lambda0=7698.974, f0=3.393e-1, gamma0=3.8e7, 
                                            n_components=2, b0=0.60, b1=0.44, N0=1.25e11, N1=1e11, v_rad0=10.5, v_rad1=11.52, v_resolution=0.56, n_step=25)
    print(wave)
    fitresult = fit_test(wavegrid=wave, ydata=normflux, restwave=7698.974, f=3.393e-1, gamma=3.8e7, 
                       b=b_o_per, N=N_o_per, v_rad=v_rad_o_per, v_resolution=0.56, n_step=25)

    plt.plot(wave, normflux, color='black', marker="s", fillstyle='none')
    plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
    plt.xlim(7698.96,7699.65)
    #plt.plot(wave, welty_fit, color='r')
    #plt.plot(wave, result.init_fit, color='orange')
    #plt.plot(wave, result.best_fit, color='g')
    #plt.plot(wave, testmodel, color='limegreen')
    plt.plot(wave, welty_fit, color='red')
    plt.plot(wave, fitresult.best_fit, color='b')
    plt.show()

    # Now let's try a doublet -- in this case, the Na I lines, again for o Per. 

