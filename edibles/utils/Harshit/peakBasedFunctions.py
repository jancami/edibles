# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.8
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# +
#importing necessary libraries

import numpy as np
import matplotlib.pyplot as plt
from lmfit.models import VoigtModel
import copy


# +
#fits voigt profile to individual peaks
#
#parameters ->
#peakData1 - data to fit voigt profile to (Nx2 array)
#sd1 - standard deviation value to be given if needed (if not needed, simply don't give that argument)
#absOrEm - 0 (default) if you want to fit an absorption line, 1 if you want to fit an emission line
#plot - 0 (default) if you don't want to plot anything, 1 if you want to plot only fit, 2 if you want to plot both data and fit
#retMod - False (default) if you want only parameters returned, True if you want to return the fitted model itself, 

def voigtUniPeak(peakData1, sd1 = 1, absOrEm = 0, plot = 2, retMod = False):
    assert absOrEm == 0 or absOrEm == 1, 'Please enter valid value for parameter absOrEm (0 or 1)'
    assert plot == 0 or plot == 1 or plot == 2, 'Please enter valid value for parameter plot (0, 1 or 2)'
    
    yForFit1 = (1 - 2*absOrEm)*(1 - peakData1[:, 1])
    
    xForFit1 = copy.deepcopy(peakData1[:, 0])
    mod1 = VoigtModel()
    params1 = mod1.guess(yForFit1, x = xForFit1)
    res1 = mod1.fit(yForFit1, params1, x = xForFit1, weights = 1/sd1)
    
    if plot == 1:
        plt.plot(xForFit1, (1 - res1.best_fit), label = 'Voigt profile fit')
        plt.legend()
    elif plot == 2:
        plt.plot(xForFit1, (1 - yForFit1), label = 'Data')
        plt.plot(xForFit1, (1 - res1.best_fit), label = 'Voigt profile fit')
        plt.legend()
    
    if retMod:
        return res1
    else:
        return {'Centre': res1.params['center'].value, 
                'FWHM': res1.params['fwhm'].value, 
                'ChiSq': res1.chisqr, 
                'RedChiSq': res1.redchi}


# +
#importing data

#data5000 = np.loadtxt('Pentacene_air_snr5000.txt')

# +
#peak parameters found by hand (required only for testing purposes)

"""
peakmins = np.zeros(5)
peakmaxs = np.zeros(5)
peakcens = np.zeros(5)

peakmins[0] = 5264.222303
peakmaxs[0] = 5270.360895
peakcens[0] = 5267.600206

peakmins[1] = 5287.032308
peakmaxs[1] = 5290.810678
peakcens[1] = 5288.886049

peakmins[2] = 5303.716283
peakmaxs[2] = 5308.335695
peakcens[2] = 5305.721707

peakmins[3] = 5335.136806
peakmaxs[3] = 5342.695359
peakcens[3] = 5337.921532

peakmins[4] = 5358.495471
peakmaxs[4] = 5365.142247
peakcens[4] = 5361.090278
"""

# +
#extracting wavelengths and intensities within ranges from lab data

"""
peak1data = data5000[np.logical_and(data5000[:, 0]>=peakmins[0], data5000[:, 0]<=peakmaxs[0])]
peak2data = data5000[np.logical_and(data5000[:, 0]>=peakmins[1], data5000[:, 0]<=peakmaxs[1])]
peak3data = data5000[np.logical_and(data5000[:, 0]>=peakmins[2], data5000[:, 0]<=peakmaxs[2])]
peak4data = data5000[np.logical_and(data5000[:, 0]>=peakmins[3], data5000[:, 0]<=peakmaxs[3])]
peak5data = data5000[np.logical_and(data5000[:, 0]>=peakmins[4], data5000[:, 0]<=peakmaxs[4])]
"""


# +
#peak1Params = voigtUniPeak(peak1data, 0.0002402523653397399)
#print(peak1Params)

# +
#print(np.linspace(0, data5000[:, 0].size, 6, dtype = int))
#for i in range(6): print(i)
#print(data5000[0:7759, 0].size)

# +
#fits multiple peaks given the number of peaks and the standard deviation (noise) in the data

def voigtMultiPeakNG(peakData2, nosPeak1, sd2):
    if (not isinstance(nosPeak1, int)) or nosPeak1 < 1:
        print('Please enter valid number of peaks (>=1)')
    else:
        xForFit2 = peakData2[:, 0]
        yForFit2 = 1 - peakData2[:, 1]
        
        idxs1 = np.linspace(0, peakData2[:, 0].size, (nosPeak1+1), dtype = int)
        #print(idxs1)
        
        vmarr1 = np.empty(shape = nosPeak1, dtype = object)
        
        vmarr1[0] = VoigtModel(prefix = 'VM1_')
        params2 = vmarr1[0].guess(yForFit2[idxs1[0]:idxs1[1]], x = xForFit2[idxs1[0]:idxs1[1]])
        
        if nosPeak1 > 1:
            for it1 in range(nosPeak1-1):
                pref1 = 'VM' + str(it1+2) + '_'
                vmarr1[it1+1] = VoigtModel(prefix = pref1)
                params2.update(vmarr1[it1+1].guess(yForFit2[idxs1[it1+1]:idxs1[it1+2]], x = xForFit2[idxs1[it1+1]:idxs1[it1+2]]))
            
            mod2 = np.sum(vmarr1)
        else:
            mod2 = vmarr1[0]
        
        res2 = mod2.fit(yForFit2, params2, x = xForFit2, weights = 1/sd2)
        plt.plot(xForFit2, (1 - yForFit2), label = 'Data')
        plt.plot(xForFit2, (1 - res2.best_fit), label = 'Multipeak Voigt profile fit')
        plt.legend()
        
        tbr1 = {}
        
        for it2 in range(nosPeak1):
            cenkey1 = 'Centre' + str(it2+1)
            cenval1 = 'VM' + str(it2+1) + '_center'
            fwhmkey1 = 'FWHM' + str(it2+1)
            fwhmval1 = 'VM' + str(it2+1) + '_fwhm'
            tbr1.update({cenkey1: res2.params[cenval1].value})
            tbr1.update({fwhmkey1: res2.params[fwhmval1].value})
        
        tbr1.update({'ChiSq': res2.chisqr})
        tbr1.update({'RedChiSq': res2.redchi})
        
        return tbr1


# +
#allPeakParams = voigtMultiPeakNG(data5000, 5, 0.0002402523653397399)
#print(allPeakParams)
# +
#fits multiple peaks given the number of peaks, the cuts for the peaks (1D array) 
#and the standard deviation (noise) in the data

def voigtMultiPeakCuts(peakData3, nosPeak2, cuts1, sd3):
    if (not isinstance(nosPeak2, int)) or nosPeak2 < 1:
        print('Please enter valid number of peaks (>=1)')
    elif not (cuts1.size == nosPeak2 + 1) :
        print('Please provide right number of cuts (number of peaks + 1)')
    else:
        xForFit3 = peakData3[:, 0]
        yForFit3 = 1 - peakData3[:, 1]
        
        #sorter1 = np.argsort(peakData3[:, 0])
        #idxs2 = sorter[np.searchsorted(peakData3[:, 0], cuts1, sorter = sorter1)]
        #print(idxs2)
        xRanges1 = np.empty(shape = nosPeak2, dtype = object)
        yRanges1 = np.empty(shape = nosPeak2, dtype = object) 
        for it5 in range(nosPeak2):
            xRanges1[it5] = xForFit3[np.logical_and(xForFit3 >= cuts1[it5], xForFit3 <= cuts1[it5+1])]
            yRanges1[it5] = yForFit3[np.logical_and(xForFit3 >= cuts1[it5], xForFit3 <= cuts1[it5+1])]
        #print(xRanges1)
        
        vmarr2 = np.empty(shape = nosPeak2, dtype = object)
        
        vmarr2[0] = VoigtModel(prefix = 'VM1_')
        params3 = vmarr2[0].guess(yRanges1[0], x = xRanges1[0])
        
        if nosPeak2 > 1:
            for it3 in range(nosPeak2-1):
                pref2 = 'VM' + str(it3+2) + '_'
                vmarr2[it3+1] = VoigtModel(prefix = pref2)
                params3.update(vmarr2[it3+1].guess(yRanges1[it3], x = xRanges1[it3]))
            
            mod3 = np.sum(vmarr2)
        else:
            mod3 = vmarr2[0]
        
        res3 = mod3.fit(yForFit3, params3, x = xForFit3, weights = 1/sd3)
        plt.plot(xForFit3, (1 - yForFit3), label = 'Data')
        plt.plot(xForFit3, (1 - res3.best_fit), label = 'Multipeak Voigt profile fit')
        plt.legend()
        
        tbr2 = {}
        
        for it4 in range(nosPeak2):
            cenkey2 = 'Centre' + str(it4+1)
            cenval2 = 'VM' + str(it4+1) + '_center'
            fwhmkey2 = 'FWHM' + str(it4+1)
            fwhmval2 = 'VM' + str(it4+1) + '_fwhm'
            tbr2.update({cenkey2: res3.params[cenval2].value})
            tbr2.update({fwhmkey2: res3.params[fwhmval2].value})
        
        tbr2.update({'ChiSq': res3.chisqr})
        tbr2.update({'RedChiSq': res3.redchi})
        
        return tbr2
# +
#fits multiple peaks given the number of peaks, the ranges for the peaks (2D array) 
#and the standard deviation (noise) in the data

def voigtMultiPeakRanges(peakData3, nosPeak2, ranges1, sd3):
    if (not isinstance(nosPeak2, int)) or nosPeak2 < 1:
        print('Please enter valid number of peaks (>=1)')
    elif not (ranges1.shape[0] == nosPeak2) :
        print('Please provide right number of range (= number of peaks)')
    else:
        xForFit3 = peakData3[:, 0]
        yForFit3 = 1 - peakData3[:, 1]
        
        #sorter1 = np.argsort(peakData3[:, 0])
        #idxs2 = sorter[np.searchsorted(peakData3[:, 0], cuts1, sorter = sorter1)]
        #print(idxs2)
        xRanges1 = np.empty(shape = nosPeak2, dtype = object)
        yRanges1 = np.empty(shape = nosPeak2, dtype = object) 
        for it5 in range(nosPeak2):
            xRanges1[it5] = xForFit3[np.logical_and(xForFit3 >= ranges1[it5, 0], xForFit3 <= ranges1[it5, 1])]
            yRanges1[it5] = yForFit3[np.logical_and(xForFit3 >= ranges1[it5, 0], xForFit3 <= ranges1[it5, 1])]
        #print(xRanges1)
        
        vmarr2 = np.empty(shape = nosPeak2, dtype = object)
        
        vmarr2[0] = VoigtModel(prefix = 'VM1_')
        params3 = vmarr2[0].guess(yRanges1[0], x = xRanges1[0])
        params3['VM1_center'].set(min = ranges1[0, 0], max = ranges1[0, 1])
        
        if nosPeak2 > 1:
            for it3 in range(nosPeak2-1):
                pref2 = 'VM' + str(it3+2) + '_'
                vmarr2[it3+1] = VoigtModel(prefix = pref2)
                params3.update(vmarr2[it3+1].guess(yRanges1[it3], x = xRanges1[it3]))
                params3[pref2+'center'].set(min = ranges1[it3, 0], max = ranges1[it3, 1])
            mod3 = np.sum(vmarr2)
        else:
            mod3 = vmarr2[0]
        
        res3 = mod3.fit(yForFit3, params3, x = xForFit3, weights = 1/sd3)
        plt.plot(xForFit3, (1 - yForFit3), label = 'Data')
        plt.plot(xForFit3, (1 - res3.best_fit), label = 'Multipeak Voigt profile fit')
        plt.legend()
        
        tbr2 = {}
        
        for it4 in range(nosPeak2):
            cenkey2 = 'Centre' + str(it4+1)
            cenval2 = 'VM' + str(it4+1) + '_center'
            fwhmkey2 = 'FWHM' + str(it4+1)
            fwhmval2 = 'VM' + str(it4+1) + '_fwhm'
            tbr2.update({cenkey2: res3.params[cenval2].value})
            tbr2.update({fwhmkey2: res3.params[fwhmval2].value})
        
        tbr2.update({'ChiSq': res3.chisqr})
        tbr2.update({'RedChiSq': res3.redchi})
        
        return tbr2
# +
#fits many single peaks at once in the given ranges of wavelengths

def voigtNUniPeak(peakData4, ranges2, sd4):
    tbr3 = {}
    
    N1 = ranges2.shape[0]
    
    for it6 in range(N1):
        #hol1 = voigtUniPeak(peakData4[np.logical_and(peakData4[:, 0] >= ranges2[it6, 0], 
        #                                             peakData4[:, 0] <= ranges2[it6, 1])], sd4)
        yForFit1 = 1 - peakData4[np.logical_and(peakData4[:, 0] >= ranges2[it6, 0], peakData4[:, 0] <= ranges2[it6, 1]), 1]
        xForFit1 = peakData4[np.logical_and(peakData4[:, 0] >= ranges2[it6, 0], peakData4[:, 0] <= ranges2[it6, 1]), 0]
        mod1 = VoigtModel()
        params1 = mod1.guess(yForFit1, x = xForFit1)
        res1 = mod1.fit(yForFit1, params1, x = xForFit1, weights = 1/sd4)
        #plt.plot(xForFit1, (1 - yForFit1), label = 'Data')
        plt.plot(xForFit1, (1 - res1.best_fit), label = 'Fitted peak ' + str(it6+1))
        cenkey3 = 'Centre' + str(it6+1)
        fwhmkey3 = 'FWHM' + str(it6+1)
        chikey1 = 'ChiSq' + str(it6+1)
        rchikey1 = 'RedChiSq' + str(it6+1)
        tbr3.update({cenkey3: res1.params['center'].value})
        tbr3.update({fwhmkey3: res1.params['fwhm'].value})
        tbr3.update({chikey1: res1.chisqr})
        tbr3.update({rchikey1: res1.redchi})
    
    return tbr3
# -

