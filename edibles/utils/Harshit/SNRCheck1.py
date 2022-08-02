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

import numpy as np
import matplotlib.pyplot as plt
from edibles import PYTHONDIR
from astropy.modeling.functional_models import Voigt1D
import copy
from stackingFunctions import widthNormLinStacker
from stackingFunctions import stackCheck
from peakBasedFunctions import voigtUniPeak

# +
SNR = 20.0
signal = 1.0
minWave = 3000.0
maxWave = 5500.0

reso = 70000.0

currX = minWave
x = np.array([minWave])

while currX < maxWave:
    currX = currX + currX/reso
    x = np.append(x, currX)

#x = np.linspace(3000.0, 5500.0, 100000)
y = np.random.normal(signal, 1/SNR, x.shape)

fig1, ax1 = plt.subplots(figsize = (12,5))

ax1.plot(x, y, label = 'Random noise with SNR ' + str(SNR))
ax1.set(xlabel = r'Wavelength ($\AA$)', ylabel = 'Relative flux')
ax1.legend()
#plt.close(fig1)
fig1.show()
#plt.close(fig1)
# -

print(np.std(y))

# +
#select molecule as follows -
# 2-Methylnaphthalene -> 1
# Acenaphthene -> 2
# Benzo[ghi]perylene -> 3
# Pentacene -> 4
# Perylene -> 5
# Phenanthrene -> 6
# Pyrene -> 7
# Phenalenyl -> 8
# Any other molecule -> 0 and input the molecule name (according to parameters file)
moleculeNo = 5

#getting molecule name from given molecule number above
#if on jupyter, just run this part

if moleculeNo == 1:
    molName = '2MethylNaphthalene'
elif moleculeNo == 2:
    molName = 'Acenaphthene'
elif moleculeNo == 3:
    molName = 'Benzoghiperylene'
elif moleculeNo == 4:
    molName = 'Pentacene'
elif moleculeNo == 5:
    molName = 'Perylene'
elif moleculeNo == 6:
    molName = 'Phenanthrene'
elif moleculeNo == 7:
    molName = 'Pyrene'
elif moleculeNo == 8:
    molName = 'Phenalenyl'
elif moleculeNo == 0:
    molName = input('Enter molecule file name (as in parameters file):\n')

# +
#loading the parameters of the given molecule
#if on jupyter, just run this part

paramFile = PYTHONDIR + '\\utils\\Harshit\\Lab Spectra Parameters\\' + molName + 'Params.txt'
molParam = np.loadtxt(paramFile)
P = molParam.shape[0]
# -

print(molParam)

"""
x1 = np.linspace(-5, 5, 100)
fwhm = 1
fL = 2*fwhm/3.6013
fG = 2.355*fwhm/3.6013
v1 = Voigt1D(x_0 = 3, fwhm_L = fL, fwhm_G = fG)
y1 = v1(x1)/0.6556557489878598
print(v1(3))

plt.plot(x1, y1, label = 'Voigt test')
plt.legend()
"""

# +
iters = 10

fig2, axs1 = plt.subplots(iters, 2, figsize = (12, 5*iters))
plt.tight_layout(rect=[0, 0.03, 1, 0.95])

for it2 in range(iters):
    y = np.random.normal(signal, 1/SNR, x.shape)
    
    #strPeak = np.random.randint(P)
    ynew = copy.deepcopy(y)
    
    if(iters < 2):
        desax1 = axs1[0]
        desax2 = axs1[1]
    else:
        desax1 = axs1[it2, 0]
        desax2 = axs1[it2, 1]
    
    for it1 in range(P):
        fwhm1 = molParam[it1, 1]
        fL1 = 2*fwhm1/3.6013
        fG1 = 2.355*fwhm1/3.6013
        v2 = Voigt1D(x_0 = molParam[it1, 0], fwhm_L = fL1, fwhm_G = fG1)
        #print(v2(molParam[it1, 0]))
        """
        if it1 < strPeak:
            curr = np.random.randint(P)
        elif it1 == strPeak:
            curr = P
        else:
            curr = np.random.randint(P+1)
        """
        #curr = np.random.random_sample() #sets all peaks to have random relative depths (max is 1 open)
        curr = 1.0 #sets all peak depths to be equal
        ynew = ynew - 0.01*curr*v2(x)/v2(molParam[it1, 0])

    desax1.plot(x, ynew, label = 'Full simulated spectrum')
    desax1.legend()
    tit1 = 'Simulated spectrum ' + str(it2+1)
    desax1.set_title(tit1)
    desax1.set(xlabel = r'Wavelength ($\AA$)', ylabel = 'Relative flux')
    desax1.legend()
    
    data = np.array([x, ynew]).transpose()
    #print(data.shape)
    
    stack = widthNormLinStacker(data, molParam, hide = True, extent = 5)
    
    dict1 = stackCheck(stack, flatReg = [-2.5, 2.5])
    
    res1 = dict1['Voigt model']
    res2 = dict1['Null model']
    rc1 = dict1['Red chi of voigt']
    print('Reduced chi square for voigt fit (iteration ' + str(it2+1) + ') is ' + str(rc1))
    rc2 = dict1['Red chi of null']
    print('Reduced chi square for null hypothesis (iteration ' + str(it2+1) + ') is ' + str(rc2))
    BCI1 = dict1['BIC of voigt']
    print('BIC for voigt fit (iteration ' + str(it2+1) + ') is ' + str(BCI1))
    BCI2 = dict1['BIC of null']
    print('BIC for null hypothesis (iteration ' + str(it2+1) + ') is ' + str(BCI2))
    fval = dict1['f']
    print('f value is ' + str(fval))
    print('****')
    #lk = dict1['Likelihood']
    #print('Likelihood of it being a detection is ' + str(lk*100) + r'%')
    
    #suppSD = np.std(y)/np.sqrt(P)
    #print('Supposed uncertainity is ' + str(suppSD))
    #calcSD = np.std(stack[np.logical_or(stack[:, 0] > 2.5, stack[:, 0] < -2.5), 1])
    #print('Calculated uncertainity is ' + str(calcSD))
    #res = voigtUniPeak(stack, sd = calcSD, plot = 0, retMod = True, centre = 0.0, sigma = 2/3.6013)
    #print('Reduced chi square for voigt fit (iteration ' + str(it2+1) + ') is ' + str(res.redchi))
    #res2 = voigtUniPeak(stack, sd = calcSD, plot = 0, retMod = True, amp = 0.0)
    #print('Reduced chi square for null hypothesis (iteration ' + str(it2+1) + ') is ' + str(res2.redchi))
    
    #likelihood = np.exp(-res.redchi)/(np.exp(-res.redchi) + np.exp(-res2.redchi))
    #print('Likelihood of it being a detection is ' + str(likelihood*100) + r'%')
    
    desax2.plot(stack[:, 0], stack[:, 1], label = 'Stack')
    desax2.plot(stack[:, 0], (1 - res1.best_fit), label = 'Voigt fit')
    desax2.plot(stack[:, 0], (1 - res2.best_fit), label = 'Null hypothesis')
    lowerY = np.ones(stack[:, 0].shape) - np.std(stack[np.logical_or(stack[:, 0] > 2.5, stack[:, 0] < -2.5), 1])
    upperY = np.ones(stack[:, 0].shape) + np.std(stack[np.logical_or(stack[:, 0] > 2.5, stack[:, 0] < -2.5), 1])
    desax2.plot(stack[:, 0], lowerY, 'r', label = 'One sigma limit')
    desax2.plot(stack[:, 0], upperY, 'r')
    desax2.legend()
    tit2 = 'Stack for simulated spectrum ' + str(it2+1)
    desax2.set_title(tit2)
    desax2.set(xlabel = 'Relative wavelength', ylabel = 'Relative flux')
    desax2.legend()

plt.subplots_adjust(hspace=0.3, wspace=0.2)


# -



# +
#check for all SNRs for with spectra with single molecules only

SNRs = [1, 2, 5, 10, 20, 50, 100, 200, 500]
mols = ['2MethylNaphthalene', 'Acenaphthene', 'Benzoghiperylene', 'Pentacene', 'Perylene', 'Phenanthrene', 'Pyrene', 'Phenalenyl']
indIter = 10

avgLk = np.zeros((len(SNRs), len(mols)))

for it5, mol in enumerate(mols):
    parFile = PYTHONDIR + '\\utils\\Harshit\\Lab Spectra Parameters\\' + molName + 'Params.txt'
    molPar = np.loadtxt(parFile)
    peaks = molPar.shape[0]
    for it6, snr in enumerate(SNRs):
        for it3 in range(indIter):
            y = np.random.normal(signal, 1/snr, x.shape)
            ynew = copy.deepcopy(y)
            for it4 in range(peaks):
                fwhm = molPar[it4, 1]
                fL = 2*fwhm/3.6013
                fG = 2.355*fwhm/3.6013
                voi = Voigt1D(x_0 = molPar[it4, 0], fwhm_L = fL, fwhm_G = fG)
                ynew = ynew - 0.01*voi(x)/voi(molPar[it4, 0])
            
            data = np.array([x, ynew]).transpose()
            stack = widthNormLinStacker(data, molParam, hide = True, extent = 5)
            dict2 = stackCheck(stack, flatReg = [-2.5, 2.5], retMods = False)
            like = dict2['Likelihood']
            avgLk[it6, it5] = avgLk[it6, it5] + like/indIter
            print('Doing mol ' + str(it5) + ', SNR ' + str(it6))
# -

plt.plot(SNRs, avgLk[:, 4], 'o-')

print(avgLk[:, 4])

# +
#check for all SNRs for with spectra with single molecules only

SNRs = [10, 12, 14, 16, 20]
mols = ['2MethylNaphthalene', 'Acenaphthene', 'Benzoghiperylene', 'Pentacene', 'Perylene', 'Phenanthrene', 'Pyrene', 'Phenalenyl']
indIter = 30

fails = np.zeros((len(SNRs), len(mols)))

for it7, mol in enumerate(mols):
    parFile = PYTHONDIR + '\\utils\\Harshit\\Lab Spectra Parameters\\' + molName + 'Params.txt'
    molPar = np.loadtxt(parFile)
    peaks = molPar.shape[0]
    for it8, snr in enumerate(SNRs):
        for it9 in range(indIter):
            y = np.random.normal(signal, 1/snr, x.shape)
            ynew = copy.deepcopy(y)
            for it10 in range(peaks):
                fwhm = molPar[it10, 1]
                fL = 2*fwhm/3.6013
                fG = 2.355*fwhm/3.6013
                voi = Voigt1D(x_0 = molPar[it10, 0], fwhm_L = fL, fwhm_G = fG)
                ynew = ynew - 0.01*voi(x)/voi(molPar[it10, 0])
            
            data = np.array([x, ynew]).transpose()
            stack = widthNormLinStacker(data, molParam, hide = True, extent = 5)
            dict2 = stackCheck(stack, flatReg = [-2.5, 2.5], retMods = False)
            bic1 = dict2['BCI of voigt']
            bic2 = dict2['BCI of null']
            fails[it8, it7] = fails[it8, it7] + int(bic2 < bic1)
            print('Doing mol ' + str(it7) + ', SNR ' + str(it8) + ', iteration ' + str(it9))
# -

plt.plot(SNRs, fails[:, 5], 'o-')

fails

np.sum(fails, axis = 1)/(len(mols)*indIter)

# +
#check for all SNRs for with spectra with single molecules only

SNRs = [20, 30, 40, 50]
mols = ['2MethylNaphthalene', 'Acenaphthene', 'Benzoghiperylene', 'Pentacene', 'Perylene', 'Phenanthrene', 'Pyrene', 'Phenalenyl']
indIter = 30

fails2 = np.zeros((len(SNRs), len(mols)))

for it7, mol in enumerate(mols):
    parFile = PYTHONDIR + '\\utils\\Harshit\\Lab Spectra Parameters\\' + molName + 'Params.txt'
    molPar = np.loadtxt(parFile)
    peaks = molPar.shape[0]
    for it8, snr in enumerate(SNRs):
        for it9 in range(indIter):
            y = np.random.normal(signal, 1/snr, x.shape)
            ynew = copy.deepcopy(y)
            for it10 in range(peaks):
                fwhm = molPar[it10, 1]
                fL = 2*fwhm/3.6013
                fG = 2.355*fwhm/3.6013
                voi = Voigt1D(x_0 = molPar[it10, 0], fwhm_L = fL, fwhm_G = fG)
                ynew = ynew - 0.01*voi(x)/voi(molPar[it10, 0])
            
            data = np.array([x, ynew]).transpose()
            stack = widthNormLinStacker(data, molParam, hide = True, extent = 5)
            dict2 = stackCheck(stack, flatReg = [-2.5, 2.5], retMods = False)
            bic1 = dict2['BCI of voigt']
            bic2 = dict2['BCI of null']
            fails2[it8, it7] = fails2[it8, it7] + int(bic2 < bic1)
            print('Doing mol ' + str(it7) + ', SNR ' + str(it8) + ', iteration ' + str(it9))
# -

fails2


