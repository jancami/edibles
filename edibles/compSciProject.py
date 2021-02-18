import numpy as np
import statistics
import matplotlib.pyplot as plt
import os.path
from edibles import DATADIR

from scipy.signal import find_peaks

from edibles.models import ContinuumModel, VoigtModel
from edibles.utils.edibles_spectrum import EdiblesSpectrum
from edibles.utils.edibles_oracle import EdiblesOracle
from  edibles.utils.voigt_profile import voigt_absorption_line
from edibles.optimizeNA import NAModelOptimize

def calcPeaks(data):
    return statistics.median(data)

def PlotFile(filename):

    sp = EdiblesSpectrum(filename)
    median = calcPeaks(sp.flux)
    plt.plot(sp.wave, sp.flux)
    plt.plot(sp.wave, [median] * len(sp.wave))
    plt.show()

def oneList(list1):
    wrange1 = [3301.5, 3304.0]

    for x in range(0,len(list1)):
        if os.path.isfile(DATADIR + list1[x]):
            sp1 = EdiblesSpectrum(list1[x])

            wave1 = sp1.wave
            flux1 = sp1.flux
            idx = np.where((wave1 > wrange1[0]) & (wave1 < wrange1[1]))

            wave = list(wave1[idx])
            flux = list(flux1[idx])/ np.median(list(flux1[idx]))

            return wave, flux
    print("failed")
    return [],[]

def normailzeLists(list1,list2):
    wrange1 = [3301.5, 3304.0]
    wrange2 = [5888.5, 5899.0]

    for x in range(0,min(len(list1),len(list2))):
        if os.path.isfile(DATADIR + list1[x]) and os.path.isfile(DATADIR + list2[x]):
            sp1 = EdiblesSpectrum(list1[x])
            sp2 = EdiblesSpectrum(list2[x])

            wave1 = sp1.bary_wave
            wave2 = sp2.bary_wave

            idx = np.where((wave1 > wrange1[0]) & (wave1 < wrange1[1]))
            idx2 = np.where((wave2 > wrange2[0]) & (wave2 < wrange2[1]))

            flux1 = sp1.flux[idx]/ np.median(sp1.flux[idx])
            flux2 = sp2.flux[idx2]/ np.median(sp2.flux[idx2])

            wave = list(wave1[idx]) +  list(wave2[idx2])
            flux = list(flux1) + list(flux2)

            return wave, flux
    print("failed")
    return [],[]

def NaModel2(filename):
    method = 'least_squares'

    lambda0 = [3302.369, 3302.978, 5889.950, 5895.924]
    f = [8.26e-03, 4.06e-03, 6.49E-01,3.24E-01]
    gamma = [6.280e7, 6.280e7,6.280e7,6.280e7]

    v_resolution = 5.75

    v_rad = [18.9]
    b = [1.1]
    N = [18.9e13]

    print(filename)
    pythia = EdiblesOracle()

    List1 = pythia.getFilteredObsList(
            object=[filename], MergedOnly=False, Wave=3302.0
    ).values.tolist()

    List2 = pythia.getFilteredObsList(
            object=[filename], MergedOnly=False, Wave=5889.0
    ).values.tolist()

    wave, flux = normailzeLists(List1,List2)
    #wave, flux = oneList(List1)

    model = NAModelOptimize()
    print("guess")
    pars = model.guess(flux, np.array(v_rad).flatten(), b, N, wavegrid = wave, lambda0=lambda0, f=f, gamma=gamma, v_resolution=v_resolution)
    print("fit")
    result = model.fit(data=flux, params=pars, wavegrid=wave, lambda0=lambda0, f=f, gamma=gamma, v_resolution=v_resolution, method=method)
    #fit can have weight argument
    #cross-correlation to find dips
    print("Report start")
    print(result.fit_report())
    print("Report End")

    AbsorptionLine = voigt_absorption_line(
        wave,
        lambda0=lambda0,
        b=result.params['b'].value,
        N=result.params['N'].value,
        f=f,
        gamma=gamma,
        v_rad=result.params['v_rad'].value,
        v_resolution=v_resolution,
    )
    test = flux - AbsorptionLine

    lambda0 = [3302.369, 3302.978, 5889.950, 5895.924].extend(lambda0)
    f = [8.26e-03, 4.06e-03, 6.49E-01,3.24E-01].extend(f)
    gamma = [6.280e7, 6.280e7,6.280e7,6.280e7].extend(gamma)

    v_resolution = 5.75

    v_rad = np.array([result.params['v_rad'].value,result.params['v_rad'].value]).flatten()
    b = np.array([result.params['b'].value,result.params['b'].value]).flatten()
    N = np.array([result.params['N'].value,result.params['N'].value]).flatten()

    model = NAModelOptimize()
    print("guess")
    pars = model.guess(test, v_rad, b, N, wavegrid = wave, lambda0=lambda0, f=f, gamma=gamma, v_resolution=v_resolution)
    print("fit")
    result = model.fit(data=flux, params=pars, wavegrid=wave, lambda0=lambda0, f=f, gamma=gamma, v_resolution=v_resolution, method=method)
    #fit can have weight argument

    AbsorptionLine = voigt_absorption_line(
        wave,
        lambda0=lambda0,
        b=result.params['b'].value,
        N=result.params['N'].value,
        f=f,
        gamma=gamma,
        v_rad=result.params['v_rad'].value,
        v_resolution=v_resolution,
    )

    plt.plot(wave, flux)
    plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
    plt.plot(wave, AbsorptionLine, color="orange", marker="*")
    plt.show()

#PlotFile("/HD170740/RED_860/HD170740_w860_redl_20140915_O12.fits")
NaModel2("HD 183143")
