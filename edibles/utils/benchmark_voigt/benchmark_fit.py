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
from edibles.utils.voigt_profile import *

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
    welty_fit=voigt_absorption_line(wave,lambda0=7698.974,f=3.393e-1,gamma=3.8e7,b=b_o_per,N=N_o_per,v_rad=v_rad_o_per,v_resolution=0.56,n_step=25)

    # Then let's see if we can fit a single Voigt profile using the LMFIT model class. 
    model=Model(voigt_absorption_line, independent_vars=['wavegrid'])
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
    #print(wave)
    fitresult = fit_multi_voigt_absorptionlines(wavegrid=wave, ydata=normflux, restwave=7698.974, f=3.393e-1, gamma=3.8e7, 
                       b=b_o_per, N=N_o_per, v_rad=v_rad_o_per, v_resolution=0.56, n_step=25)

    #plt.plot(wave, normflux, color='black', marker="s", fillstyle='none')
    #plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
    #plt.xlim(7698.96,7699.65)
    ##plt.plot(wave, welty_fit, color='r')
    ##plt.plot(wave, result.init_fit, color='orange')
    ##plt.plot(wave, result.best_fit, color='g')
    ##plt.plot(wave, testmodel, color='limegreen')
    #plt.plot(wave, welty_fit, color='red')
    #plt.plot(wave, fitresult.best_fit, color='b')
    #plt.show()

    # Now let's try to fake a doublet, and use 4 cloud components. 
    b_fake = [0.60, 0.44, 0.62, 0.60]
    N_fake = [12.5e10, 10e10, 22.5e10, 3.9e10]
    v_rad_fake = [10.5+0.15, 11.52+0.15, 14.74+0.15, 15.72+0.15]
    fitresult = fit_multi_voigt_absorptionlines(wavegrid=wave, ydata=normflux, restwave=[7698.974, 7699.074], f=[3.393e-1, 5e-1], gamma=[3.8e7, 3.8e7], 
                       b=b_fake, N=N_fake, v_rad=v_rad_fake, v_resolution=0.56, n_step=25)
    fitresult.params.pretty_print()
    plt.plot(wave, normflux, color='black', marker="s", fillstyle='none')
    plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
    plt.xlim(7698.96,7699.65)
    plt.plot(wave, welty_fit, color='red')
    plt.plot(wave, fitresult.best_fit, color='b')
    plt.show()

    # The real doublet test: let's try lithium. 
    #7LiI,6707.761,0.4982,3.69E+07,10
    #7LiI,6707.912,0.2491,3.69E+07,10
    #6LiI,6707.921,0.4982,3.69E+07,10
    #6LiI,6708.072,0.2491,3.69E+07,10
    # We really need to consider the 2 doublets. The problem that we will have is that
    # the 7Li and 6Li should have different column densities, but the same b and v_rad values. 
    # That's hard to code in the Voigt profile functions. An easy workaround is to 
    # "fake" this by changing the f values for 6Li and 7Li to represent the ratio of the 
    # column densities. In the example below, I divide the f values of 6Li by 5. This then
    # means that for the final fit, the 7Li column density will be given by the best-fit N value
    # (it comes out to 1.7e10), and this then means that the 6Li column 
    # density is 1.7e10/5 = 3.4e9. 
    restwave = [6707.761, 6707.912, 6707.921, 6708.072]
    f = [0.4982, 0.2491, 0.4982/5, 0.2491/5]
    gamma = [3.69e7, 3.69e7, 3.69e7, 3.69e7]
    pythia = EdiblesOracle()
    List = pythia.getFilteredObsList(object=["HD 147889"], MergedOnly=True, Wave=6708.)
    test = List.values.tolist()
    filename = test[0]
    print(filename)
    sp = EdiblesSpectrum(filename)
    wave = sp.wave
    flux = sp.flux
    wrange = [6707.0, 6709.5]
    idx = np.where((wave > wrange[0]) & (wave < wrange[1]))
    wave = wave[idx]
    flux = flux[idx]
    norm = flux / np.median(flux)

    fitresult = fit_multi_voigt_absorptionlines(wavegrid=wave, ydata=norm, restwave=restwave, f=f, gamma=gamma, 
                       b=0.4, N=1e10, v_rad=20., v_resolution=3.75, n_step=25)
    fitresult.params.pretty_print()

    plt.plot(wave,norm)
    plt.xlim(6707.0,6709.5)
    #plt.ylim(2700,2900)
    plt.plot(wave, fitresult.best_fit, color='b')
    plt.show()