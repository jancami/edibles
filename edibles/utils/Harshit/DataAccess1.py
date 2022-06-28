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
from edibles.utils.edibles_oracle import EdiblesOracle
from edibles.utils.edibles_spectrum import EdiblesSpectrum
from edibles.utils.ContinuumFitter import ContinuumFitter
from stackingFunctions import widthNormLinStacker

pythia = EdiblesOracle()
starName = 'HD 183143'
minWave = 4000
maxWave = 4200
rawList = pythia.getFilteredObsList(object = [starName], MergedOnly = True, WaveMin = minWave, WaveMax = maxWave)
fnames = rawList.tolist()

datasRaw = np.empty(shape = len(fnames), dtype = object)

# +
ffig1, faxs1 = plt.subplots(-((-len(fnames))//2), 2, figsize=(12,-5*((-len(fnames))//2)))
plt.tight_layout(rect=[0, 0.03, 1, 0.95])

for it1, file in enumerate(fnames):
    sp = EdiblesSpectrum(file)
    sp.getSpectrum(xmin = minWave, xmax = maxWave)
    datasRaw[it1] = np.array([sp.bary_wave, sp.bary_flux]).transpose()
    
    if(len(fnames) < 3):
        desax1 = faxs1[it1]
    else:
        desax1 = faxs1[it1//2, it1 - 2 * (it1//2)]
    
    desax1.plot(sp.bary_wave, sp.bary_flux, label="Barycentric")
    tit1 = 'Merged spectra of ' + starName + ', observation ' + str(it1+1) + ' (On ' + str(sp.datetime.day) + '/' + str(sp.datetime.month) + '/' + str(sp.datetime.year) + ')'
    desax1.set_title(tit1)
    desax1.set(xlabel = r'Wavelength ($\AA$)', ylabel = 'Flux')
    desax1.legend()

if not ((len(fnames) - 2 * (len(fnames)//2)) == 0):
    ffig1.delaxes(faxs1[len(fnames)//2, 1])

plt.subplots_adjust(hspace=0.3, wspace=0.2)
# -

datas = datasRaw

# +
ffig2, faxs2 = plt.subplots(-((-len(datas))//2), 2, figsize=(12,-5*((-len(datas))//2)))
plt.tight_layout(rect=[0, 0.03, 1, 0.95])

for it2 in range(len(datas)):
    CF2 = ContinuumFitter(datasRaw[it2][:, 0], datasRaw[it2][:, 1])
    CS1, contPoints1  = CF2.SplineManualAnchor()
    datas[it2][:, 1] = datasRaw[it2][:, 1]/CS1(datasRaw[it2][:, 0])
    
    if(len(datas) < 3):
        desax2 = faxs2[it2]
    else:
        desax2 = faxs2[it2//2, it2 - 2 * (it2//2)]
    
    desax2.plot(datas[it2][:, 0], datas[it2][:, 1], label="Barycentric")
    tit2 = 'Continuum divided spectra of ' + starName + ', observation ' + str(it1+1) + ' (On ' + str(sp.datetime.day) + '/' + str(sp.datetime.month) + '/' + str(sp.datetime.year) + ')'
    desax2.set_title(tit2)
    desax2.set(xlabel = r'Wavelength ($\AA$)', ylabel = 'Flux')
    desax2.legend()

if not ((len(datas) - 2 * (len(datas)//2)) == 0):
    ffig2.delaxes(faxs2[len(datas)//2, 1])

plt.subplots_adjust(hspace=0.3, wspace=0.2)
# -

perylene = np.loadtxt(r'C:\Users\hkhan\edibles\edibles\utils\Harshit\Lab Spectra Parameters\PeryleneParams.txt')
print(perylene)

# +
#make code to give ranges for regions of avoidance and and insert them in stacker
# -

#print(datas[0].shape)
#print(np.logical_and(datas[:, 0] >= 4085, datas[:, 0] <= 4115).shape)
stackedData1 = widthNormLinStacker(np.delete(datas[0], 
                                             np.logical_or(np.logical_and(datas[0][:, 0] >= 4100, datas[0][:, 0] <= 4170),
                                                           np.logical_and(datas[0][:, 0] >= 4007, datas[0][:, 0] <= 4030))
                                             , 0), perylene)

stackedData2 = widthNormLinStacker(datas[1], perylene)

plt.plot(stacked)
