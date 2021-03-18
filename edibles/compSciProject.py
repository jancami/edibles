import numpy as np
import statistics
import matplotlib.pyplot as plt
import os.path
import time
from edibles import DATADIR

from edibles.utils.edibles_spectrum import EdiblesSpectrum
from edibles.utils.edibles_oracle import EdiblesOracle
from  edibles.utils.voigt_profile import voigt_absorption_line
from edibles.optimizeNA import NAModelOptimize

def normailzeLists(list1,list2):
    wrange1 = [3301.5, 3304.0]
    wrange2 = [5888.5, 5899.0]
    file1 = None
    file2 = None

    for y in range(0, len(list1)):
        if os.path.isfile(DATADIR + list1[y]):
            file1 = list1[y]

    for x in range(0, len(list2)):
        if os.path.isfile(DATADIR + list2[x]):
            file2 = list2[x]

    if file1 and file2:
        print(list1[x],list2[x])
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
        #return wave, flux
        return list(wave1[idx]),list(flux1),list(wave2[idx2]),list(flux2)
    print("failed")
    return [],[]

def NaModel(filename):
    method = 'least_squares'
    startTime = time.time()
    lambda0 = [3302.369, 3302.978, 5889.950, 5895.924]
    f = [8.26e-03, 4.06e-03, 6.49E-01,3.24E-01]
    gamma = [6.280e7, 6.280e7,6.280e7,6.280e7]

    v_resolution = 7.5

    v_rad = np.array([18.9])
    b = np.array([1.0])
    N = np.array([1.8e13])

    print(filename)
    pythia = EdiblesOracle()

    List1 = pythia.getFilteredObsList(
            object=[filename], MergedOnly=True, Wave=3302.0
    ).values.tolist()

    List2 = pythia.getFilteredObsList(
            object=[filename], MergedOnly=True, Wave=5889.0
    ).values.tolist()

    wave, flux, wave2, flux2  = normailzeLists(List1,List2)

    temp = 1
    test = flux
    oldModel = None
    oldPars = None
    v_rad = [0.0]
    while temp < 3:
        model = NAModelOptimize(n_components=temp)
        pars = model.guess(flux,wave,v_rad)

        result = model.fit(data=test, params=pars,  x=wave)

        print("Report start")
        print(result.fit_report())
        print("Report End")

        for i in range (0, temp):
            v_rad[i] = result.params['V_off_Cloud' + str(i)].value
            v_rad.append(0.0)

        #test = flux - AbsorptionLine

        lambda0.extend([3302.369, 3302.978, 5889.950, 5895.924])
        f.extend([8.26e-03, 4.06e-03, 6.49E-01,3.24E-01])
        gamma.extend([6.280e7, 6.280e7,6.280e7,6.280e7])

        temp += 1
        if not oldModel:
            oldModel = model
            oldPars = pars
        else:
            #compareModels(testModel, oldModel)
            oldModel = model
            oldPars = pars

    print(result.covar)  #look into f test
    print('Total Time:', time.time() - startTime)

    AbsorptionLine3000 = model.eval(params=result.params, x=wave)
    AbsorptionLine5000 = model.eval(params=result.params, x=wave2)
    AbsorptionLineComb = model.eval(params=result.params, x=(wave+wave2))

    plot1 = plt.figure(1)
    plt.plot(wave, flux)
    plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
    plt.plot(wave, AbsorptionLine3000, color="orange", marker="*")

    plot2 = plt.figure(2)
    plt.plot(wave2, flux2)
    plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
    plt.plot(wave2, AbsorptionLine5000, color="orange", marker="*")

    plot3 = plt.figure(3)
    plt.plot(wave+wave2, flux+flux2)
    plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
    plt.plot(wave+wave2, AbsorptionLineComb, color="orange", marker="*")
    plt.show()


def sample(filename):
    method = 'least_squares'

    lambda0 = [3302.369, 3302.978, 5889.950, 5895.924,3302.369, 3302.978, 5889.950, 5895.924]
    f = [8.26e-03, 4.06e-03, 6.49E-01,3.24E-01,8.26e-03, 4.06e-03, 6.49E-01,3.24E-01]
    gamma = [6.280e7, 6.280e7,6.280e7,6.280e7,6.280e7, 6.280e7,6.280e7,6.280e7]

    v_resolution = 7.5

    v_rad = np.array([-13.04865242,-13.04865242,-13.04865242,-13.04865242,1.4420571,1.4420571,1.4420571,1.4420571])
    b = np.array([1.085406,1.085406,1.085406,1.085406,1.8801332,1.8801332,1.8801332,1.8801332])
    N = np.array([3.28705860e+14,3.28705860e+14,3.28705860e+14,3.28705860e+14,9.33830424e+15,9.33830424e+15,9.33830424e+15,9.33830424e+15])

    v_rad = np.array([4.50011803,4.50011803,4.50011803,4.50011803,-10.1235709,-10.1235709,-10.1235709,-10.1235709])
    b = np.array([1.13228303,1.13228303,1.13228303,1.13228303,1.38759193,1.38759193,1.38759193,1.38759193])
    N = np.array([5.0017e+15,5.0017e+15,5.0017e+15,5.0017e+15,4.6447e+15,4.6447e+15,4.6447e+15,4.6447e+15])

    print(filename)
    pythia = EdiblesOracle()

    List1 = pythia.getFilteredObsList(
            object=[filename], MergedOnly=True, Wave=3302.0
    ).values.tolist()

    List2 = pythia.getFilteredObsList(
            object=[filename], MergedOnly=True, Wave=5889.0
    ).values.tolist()

    wave, flux = normailzeLists(List1,List2)

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

#sample("HD 183143")
#NaModel("HD 183143")
#NaModel("HD 22951")
#NaModel("HD 39680")
#NaModel("HD 49787")
#NaModel("HD 45314")
NaModel("HD 54662")
#NaModel("HD 80558")
