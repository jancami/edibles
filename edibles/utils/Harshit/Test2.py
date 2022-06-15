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
from peakBasedFunctions import voigtMultiPeakRanges
from peakBasedFunctions import voigtNUniPeak
from stackingFunctions import widthNormLinStacker
# -

data2metNM = np.loadtxt('C:/Users/hkhan/edibles/edibles/data/Labdata/CRDS/2MethylNaphthalene.dat', skiprows = 1)
print(data2metNM)
plt.plot(data2metNM[:, 0], data2metNM[:, 1])

from edibles.utils.functions import vac2air_ciddor

print(np.mean(data2metNM[np.logical_and(data2metNM[:, 0] >= 312.5, data2metNM[:, 0] <= 315), 1]))

data2met = data2metNM
#print(np.mean(data2metNM[np.logical_and(data2metNM[:, 0] >= 312.5, data2metNM[:, 0] <= 315.0), 1]))
data2met[:, 0] = vac2air_ciddor(10.0*data2metNM[:, 0])
data2met[:, 1] = 1 - 0.01*((data2metNM[:, 1] - (18.92082042952208))/(np.max(data2metNM[:, 1]) - - (18.92082042952208)))

# +
# #%matplotlib inline

plt.plot(data2met[:, 0], data2met[:, 1])
# -

voigtMultiPeakRanges(data2met, 5, np.array([[310.4, 310.7],
                                            [310.8, 311.1],
                                            [311.3, 311.6],
                                            [315.1, 315.34],
                                            [315.34, 315.43]]), 0.25)

print(np.std(data2met[np.logical_and(data2met[:, 0] >= 3116, data2met[:, 0] <= 3151), 1]))

# +
# %matplotlib inline

plt.plot(data2met[:, 0], data2met[:, 1], label = 'Data')
raw2met = voigtNUniPeak(data2met, np.array([[3104, 3107],
                                            [3108, 3111],
                                            [3113, 3116],
                                            [3151, 3153.4],
                                            [3153.4, 3154.3]]), 8.481437933819905e-05)
plt.legend()
# -

print(raw2met)

# +
tPar = np.zeros((5,2))

tPar[0, 0] = raw2met['Centre1']
tPar[1, 0] = raw2met['Centre2']
tPar[2, 0] = raw2met['Centre3']
tPar[3, 0] = raw2met['Centre4']
tPar[4, 0] = raw2met['Centre5']

tPar[0, 1] = raw2met['FWHM1']
tPar[1, 1] = raw2met['FWHM2']
tPar[2, 1] = raw2met['FWHM3']
tPar[3, 1] = raw2met['FWHM4']
tPar[4, 1] = raw2met['FWHM5']
# -

spData = widthNormLinStacker(data2met, tPar)


