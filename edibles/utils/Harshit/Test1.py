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
from peakBasedFunctions import voigtMultiPeakNG
from stackingFunctions import widthNormLinStacker

# +
#loading the pentacene lab data and simulated spectrum

data5000 = np.loadtxt('Pentacene_air_snr5000.txt')
data2000 = np.loadtxt('Pentacene_air_snr2000.txt')
data1000 = np.loadtxt('Pentacene_air_snr1000.txt')
data500 = np.loadtxt('Pentacene_air_snr500.txt')
data100 = np.loadtxt('Pentacene_air_snr100.txt')

# +
pParRaw = voigtMultiPeakNG(data5000, 5, 0.0002402523653397399)

pPar = np.zeros((5,2))

pPar[0, 0] = pParRaw['Centre1']
pPar[1, 0] = pParRaw['Centre2']
pPar[2, 0] = pParRaw['Centre3']
pPar[3, 0] = pParRaw['Centre4']
pPar[4, 0] = pParRaw['Centre5']

pPar[0, 1] = pParRaw['FWHM1']
pPar[1, 1] = pParRaw['FWHM2']
pPar[2, 1] = pParRaw['FWHM3']
pPar[3, 1] = pParRaw['FWHM4']
pPar[4, 1] = pParRaw['FWHM5']
# -

spData = widthNormLinStacker(data5000, pPar)

print(spData.shape)

voigtMultiPeakNG(spData, 1, 0.0002402523653397399/np.sqrt(5))


