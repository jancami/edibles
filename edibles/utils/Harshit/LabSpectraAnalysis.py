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
from edibles.utils.ContinuumFitter import ContinuumFitter

# +
#loading raw data, change the file name for respective molecule and run this part

#put address for file
dataRaw = np.loadtxt(r'C:\Users\hkhan\edibles\edibles\data\Labdata\SchmidtLab\UNSWphenalenyl_D1.DAT', skiprows = 1)
plt.plot(dataRaw[:, 0], dataRaw[:, 1], label = 'Raw Data')
plt.legend()

# +
#workup of 1st column of data, change the option and run this part

#If in first column of file - wavelength is given in nm, set option = 0
#                           - wavelength is given in angstrom, set option = 1
#                           - wave number if given in 1/cm, set option = 2
option = 2

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
#run this part to select points to find continuum
#select points (in strict increasing wavelength order) by left clicking
#once done selecting, press centre mouse key to end selecting

CF1 = ContinuumFitter(data[:, 0], data[:, 1])
CS, contPoints  = CF1.SplineManualAnchor()


# +
#workup of 2nd column of data, just run this part
#donot run it twice in sequence!!! (run from loading of raw data if you want to run this part again)

data[:, 1] = 1 - 0.01*(dataRaw[:, 1] - CS(data[:, 0]))/(np.max(dataRaw[:, 1]) - CS(data[dataRaw[:, 1] == np.max(dataRaw[:, 1]), 0]))
plt.plot(data[:, 0], data[:, 1], label = 'Final data')
plt.legend()

# +
#run this part to select start and end points of peak (in strict increasing wavelength order)
#select only start and end points of peaks, nothing else
#make sure no. of points selected is two times no. of peaks

CF2 = ContinuumFitter(data[:, 0], data[:, 1])
wvs1 = CF2.SelectPoints(n=100, y_message = 'Select peak start and end points', vetoTimeout = True)[:, 0]
peakRanges = np.reshape(wvs1, (int(wvs1.size/2), 2))

# +
#just run this part to check peak ranges

#peakRanges = np.delete(peakRanges, 3, 0) #use to delete some peaks that just aren't worth changing
print(peakRanges)

# +
#calculation of sd here, just run this part

sdArr = data

for it2 in range(peakRanges.shape[0]):
    #print(np.logical_and(sdArr[:, 0] >= peakRanges[it2, 0], sdArr[:, 0] <= peakRanges[it2, 1]).shape)
    sdArr = np.delete(sdArr, np.logical_and(sdArr[:, 0] >= peakRanges[it2, 0], sdArr[:, 0] <= peakRanges[it2, 1]), 0)

sd = np.std(sdArr[:, 1])
print(sd)

# +
#fitting part, just run this part

plt.plot(data[:, 0], data[:, 1], label = 'Data')
rawParams = voigtNUniPeak(data, peakRanges, sd)
plt.legend()
print(rawParams)

# +
#to load raw parameters into array passable in stacker, change the file name according to molecule and run this part
#if you want to force edit, change fEdit to 1, otherwise keep it to 0

fEdit = 0

params = np.zeros((peakRanges.shape[0], 2))

for it1 in range(peakRanges.shape[0]):
    lab1 = 'Centre' + str(it1+1)
    lab2 = 'FWHM' + str(it1+1)
    params[it1, 0] = rawParams[lab1]
    params[it1, 1] = rawParams[lab2]

print(params)

#change fileName here according to molecule
fileName = 'Lab Spectra Parameters/PhenalenylParams.txt'

if (not os.path.exists(fileName)) or fEdit == 1:
    np.savetxt(fileName, params)

# +
#stacking for checking, just run this part

stackedData = widthNormLinStacker(data, params)
# -


