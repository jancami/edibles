from edibles import DATADIR
from edibles import PYTHONDIR
from edibles.utils.edibles_spectrum import EdiblesSpectrum
from edibles.utils.edibles_oracle import EdiblesOracle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from PyAstronomy import pyasl
import astropy.constants as cst

# 
labfilename = '/data/python/devel/141_REMPI_472-438.txt'
# The file contains wavenumbers (in vacuum) and intensities. Read those in as pandas dataframe. 
labspec = pd.read_csv(labfilename, delimiter=',')
# Add a column to contain the (air) wavelength in AA. 
labspec['wavelength'] = pyasl.vactoair2(1e1 * labspec['wave'],mode='ciddor')
normint = labspec.int - np.median(labspec.int)
labspec['norm'] = np.exp(-0.02 * normint / np.max(normint))

objects = ['HD 147889', 'HD 183143', 'HD 169454']

# Entire wavelength range is roughly 4350-4720AA, but strongest band near 4689 DIB. 
# Match from Tim is in HD 147889, so let's try to reproduce that first.... 
plotrange = [4670,4700]
pythia = EdiblesOracle()
#List = pythia.getFilteredObsList(object=["HD 147889"], OrdersOnly=True, Wave=4689)
List = pythia.getFilteredObsList(object=["HD 183143"], OrdersOnly=True, Wave=4689)
List = pythia.getFilteredObsList(object=["HD 169454"], OrdersOnly=True, Wave=4689)
for filename in List: 
    sp = EdiblesSpectrum(filename)
    sp.getSpectrum(np.min(sp.raw_wave)+1,np.max(sp.raw_wave)-1)
    wave = sp.bary_wave
    # Need to apply radial velocity correction to get it to rest frame. 
    v_rad = -7.8 # in km/s for HD 147889 (from Haoyu)
    v_rad = -10.4 # in km/s for HD 183143 (from Haoyu; one of 2 components)
    v_rad = -8.7 # in km/s for HD 183143 (from Haoyu; one of 2 components)
    wave_rest = wave / (1+v_rad/cst.c.to("km/s").value)
    flux = np.clip(sp.bary_flux, 0, None) 
    bool_keep = (wave_rest > plotrange[0]) & (wave < plotrange[1])
    plotwave = wave_rest[bool_keep]
    plotflux = flux[bool_keep]
    #print(plotwave)
    normflux = plotflux / np.median(plotflux)
    plt.figure()
    plt.xlim(plotrange)
    #print(normflux)
    ylim = [np.min(normflux), np.max(normflux)]
    plt.ylim(ylim)
    plt.plot(plotwave, normflux)
    # Rescale lab spectrum to plot range
    #dynrange = ylim[1]-ylim[0]
    plt.plot(labspec.wavelength, labspec.norm)
    plt.show()

#plt.plot(labspec.wavelength,labspec.int)
#plt.show()