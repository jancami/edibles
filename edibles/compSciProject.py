import numpy as np
import matplotlib.pyplot as plt
import os.path
import time
from edibles import DATADIR

from edibles.utils.edibles_spectrum import EdiblesSpectrum
from edibles.utils.edibles_oracle import EdiblesOracle
from edibles.optimizeNA import NAModelOptimize
from edibles.compareModels import CompareModelBeyes

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
        sp1 = EdiblesSpectrum(list1[x])
        sp2 = EdiblesSpectrum(list2[x])

        wave1 = sp1.bary_wave
        wave2 = sp2.bary_wave

        idx = np.where((wave1 > wrange1[0]) & (wave1 < wrange1[1]))
        idx2 = np.where((wave2 > wrange2[0]) & (wave2 < wrange2[1]))

        flux1 = sp1.flux[idx]/ np.median(sp1.flux[idx])
        flux2 = sp2.flux[idx2]/ np.median(sp2.flux[idx2])
        return list(wave1[idx]),list(flux1),list(wave2[idx2]),list(flux2)
    raise Exception("Files Not Found")

def NaModel(filename):
    startTime = time.time()

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

    wave, flux, wave2, flux2 = normailzeLists(List1,List2)

    test = flux
    print("fit 3300:")
    oldModel, oldPars, oldResult, olditer= newComponent('L_', test, wave)
    AbsorptionLine5000 = oldModel.eval(params=oldResult.params, x=(wave2))
    test = flux2 / AbsorptionLine5000
    test = [1 if x > 1 else x for x in test]

    print("fit 5800:")
    oldModel2, oldPars2, oldResult2, olditer2= newComponent('R_', test, wave2)

    print("final correction:")
    if oldModel2 and oldPars2:
        model = oldModel * oldModel2
        pars = oldPars + oldPars2
        result = model.fit(data=flux+flux2, params=pars,  x=wave+wave2)
        print(olditer,olditer2)
        if CompareModelBeyes(oldResult, result, olditer, olditer+olditer2):
            model = oldModel
            pars = oldPars
            result = oldResult

    else:
        model = oldModel
        pars = oldPars
        result = oldResult


    print('Total Time:', time.time() - startTime)

    print("Final Report start")
    print(result.fit_report())
    print("Final Report End")

    AbsorptionLineComb = model.eval(params=pars, x=(wave+wave2))

    f4 = plt.figure(2)
    plt.plot(wave+wave2, flux+flux2)
    plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
    plt.plot(wave+wave2, AbsorptionLineComb, color="orange", marker="*")
    plt.show()

def newComponent(prefix, test, wave):
    temp = 1
    oldModel = None
    oldPars = None
    v_rad = [0.0]

    while temp < 10:
        print(prefix, 'iteration:', temp)
        model = NAModelOptimize(n_components=temp, prefix=prefix)
        pars = model.guess(test, wave, v_rad)

        result = model.fit(data=test, params=pars,  x=wave)
        print(result.covar)  #look into f test

        print("Report start")
        print(result.fit_report())
        print("Report End")

        for i in range (0, temp):
            v_rad[i] = result.params[prefix+'V_off_Cloud' + str(i)].value
        v_rad.append(0.0)

        if not oldModel:
            oldModel = model
            oldPars = result.params
            oldResult = result

        else:
            if CompareModelBeyes(oldResult, result, temp-1,temp):
                break
            else:
                oldModel = model
                oldPars = result.params
                oldResult = result

        temp += 1
    return oldModel, oldPars, oldResult, temp-1

NaModel("HD 183143")
#NaModel("HD 22951")
#NaModel("HD 39680")
#NaModel("HD 49787")
#NaModel("HD 45314")
#NaModel("HD 54662")
#NaModel("HD 80558")
