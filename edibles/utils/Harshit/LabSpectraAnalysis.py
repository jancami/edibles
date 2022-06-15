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
#importing necessary libraries, just run this part

import numpy as np
import matplotlib.pyplot as plt
from peakBasedFunctions import voigtNUniPeak
from stackingFunctions import widthNormLinStacker
from edibles.utils.functions import vac2air_ciddor
import os.path

# +
#loading raw data, change the file name for respective molecule and run this part

#put address for file
dataRaw = np.loadtxt('C:/Users/hkhan/edibles/edibles/data/Labdata/CRDS/2MethylNaphthalene.dat', skiprows = 1)
plt.plot(dataRaw[:, 0], dataRaw[:, 1], label = 'Raw Data')
plt.legend()

# +
#workup of 1st column of data, change the option and run this part

#If in first column of file - wavelength is given in nm set option = 0
#                           - wavelength is given in angstrom set option = 1
#                           - wave number if given in 1/cm set option = 2
option = 0

data = np.zeros(dataRaw.shape)

if option == 0:
    data[:, 0] = vac2air_ciddor(dataRaw[:, 0]*10.0)
elif option == 1:
    data[:, 0] = vac2air_ciddor(dataRaw[:, 0])
elif option == 2:
    data[:, 0] = vac2air_ciddor((1/dataRaw[:, 0])*1e8)
    
data[:, 1] = dataRaw[:, 1]
    
plt.plot(data[:, 0], data[:, 1], label = '1st column worked up')
plt.legend()

# +
#input flat range here and then run this part

flatRange = [3120, 3150]

# +
#workup of 2nd column of data, just run this part
#donot run it twice in sequence!!! (run from loading of raw data if you want to run this part again)

# %matplotlib notebook

data[:, 1] = 1 - 0.01*(dataRaw[:, 1] - np.mean(dataRaw[np.logical_and(data[:, 0] >= flatRange[0], data[:, 0] <= flatRange[1]), 1]))/(np.max(dataRaw[:, 1]) - np.mean(dataRaw[np.logical_and(data[:, 0] >= flatRange[0], data[:, 0] <= flatRange[1]), 1]))
plt.plot(data[:, 0], data[:, 1], label = 'Final data')
plt.legend()

# +
#input rough peak ranges (min, max) here from interacting with above graph and then run this part

peakRanges = np.array([[3104, 3107],
                       [3108, 3111],
                       [3113, 3116],
                       [3151, 3153.4],
                       [3153.4, 3154.3]])

# +
#calculation of sd here, just run this part

# %matplotlib inline

sd = np.std(data[np.logical_and(data[:, 0] >= flatRange[0], data[:, 0] <= flatRange[1]), 1])
print(sd)

# +
#fitting part, just run this part

plt.plot(data[:, 0], data[:, 1], label = 'Data')
rawParams = voigtNUniPeak(data, peakRanges, sd)
plt.legend()
print(rawParams)

# +
#load raw parameters into array passable in stacker, change the file name according to molecule and run this part

params = np.zeros((peakRanges.shape[0], 2))

for it1 in range(peakRanges.shape[0]):
    lab1 = 'Centre' + str(it1+1)
    lab2 = 'FWHM' + str(it1+1)
    params[it1, 0] = rawParams[lab1]
    params[it1, 1] = rawParams[lab2]

print(params)

#change fileName here according to molecule
fileName = '2MethylNaphthaleneParams.txt'

if not os.path.exists(fileName):
    np.savetxt(fileName, params)

# +
#stacking for checking, just run this part

stackedData = widthNormLinStacker(data, params)
# -


