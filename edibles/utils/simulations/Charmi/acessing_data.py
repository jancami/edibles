# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 13:54:26 2022

@author: Charmi Bhatt
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
from edibles.utils.edibles_oracle import EdiblesOracle
from edibles.utils.edibles_spectrum import EdiblesSpectrum
import warnings
from astropy.modeling import models
from astropy import units as u
from specutils.spectra import Spectrum1D
from specutils.fitting import fit_generic_continuum
from edibles.utils.edibles_oracle import EdiblesOracle
from edibles.utils.edibles_spectrum import EdiblesSpectrum
import warnings
from astropy.modeling import models
from astropy import units as u
from specutils.spectra import Spectrum1D
from specutils.fitting import fit_generic_continuum

starName = 'HD 23180'  #166937
#put lower range of wavelengths to extract from edibles data
minWave = 6614

#put upper range of wavelengths to extract from edibles data
maxWave = 6616

pythia = EdiblesOracle()
rawList = pythia.getFilteredObsList(object = [starName], MergedOnly = True, WaveMin = minWave, WaveMax = maxWave)
fnames = rawList.tolist()
obs = len(fnames)


sp = EdiblesSpectrum(fnames[0])
    
sp.getSpectrum(xmin = max(minWave, np.min(sp.raw_wave)+1)
                   , xmax = min(maxWave, np.max(sp.raw_wave)-1))

                   
#data = np.array([sp.bary_wave, sp.bary_flux]).transpose()
leftEdge = 0
rightEdge = 0
    
if minWave <= np.min(sp.raw_wave):
    leftEdge = 1
    #print('Left edge detected')
if maxWave >= np.max(sp.raw_wave):
    rightEdge = 1

data = np.delete(np.array([sp.bary_wave, sp.bary_flux]).transpose(), 
                           np.logical_or(sp.bary_wave <= np.min(sp.bary_wave) + 40.0*leftEdge, 
                                         sp.bary_wave >= np.max(sp.bary_wave) - 40.0*rightEdge), 0)

v = -6.5

data[:, 0] = data[:, 0]*(1+v/299792.458)

x1 = data[:,0]
y1= data[:,1]

spectrum1 = Spectrum1D(flux = y1*u.dimensionless_unscaled, spectral_axis = x1*u.angstrom)

with warnings.catch_warnings():  # Ignore warnings
    warnings.simplefilter('ignore')
    g1_fit = fit_generic_continuum(spectrum1, model = models.Legendre1D(degree = 5))


data[:,1] = y1/g1_fit(x1*u.angstrom)

data[:, 0] = 1/(data[:, 0])*1e+8


plt.plot(data[:, 0], data[:, 1])
#plt.xlim(6612, 6615)
plt.ylim(0.9, 1.05)

#print(data[:,0])

for i in range(len(data[:,0])):
    print((data[:, 0][i]- data[:, 0][i+1])/data[:,0])






# '''works:'''

# starName = 'HD 166937'
#     #put lower range of wavelengths to extract from edibles data
#     minWave = 6612
    
#     #put upper range of wavelengths to extract from edibles data
#     maxWave = 6616
    
#     pythia = EdiblesOracle()
#     rawList = pythia.getFilteredObsList(object = [starName], MergedOnly = True, WaveMin = minWave, WaveMax = maxWave)
#     fnames = rawList.tolist()
#     obs = len(fnames)
    
    
#     sp = EdiblesSpectrum(fnames[0])
        
#     sp.getSpectrum(xmin = max(minWave, np.min(sp.raw_wave)+1)
#                         , xmax = min(maxWave, np.max(sp.raw_wave)-1))
    
                       
#     #data = np.array([sp.bary_wave, sp.bary_flux]).transpose()
#     leftEdge = 0
#     rightEdge = 0
        
#     if minWave <= np.min(sp.raw_wave):
#         leftEdge = 1
#         #print('Left edge detected')
#     if maxWave >= np.max(sp.raw_wave):
#         rightEdge = 1
    
#     data = np.delete(np.array([sp.bary_wave, sp.bary_flux]).transpose(), 
#                                 np.logical_or(sp.bary_wave <= np.min(sp.bary_wave) + 40.0*leftEdge, 
#                                               sp.bary_wave >= np.max(sp.bary_wave) - 40.0*rightEdge), 0)
    
#     v = -6.5
    
#     data[:, 0] = data[:, 0]*(1+v/299792.458)
    
#     x1 = data[:,0]
#     y1= data[:,1]
    
#     spectrum1 = Spectrum1D(flux = y1*u.dimensionless_unscaled, spectral_axis = x1*u.angstrom)
    
#     with warnings.catch_warnings():  # Ignore warnings
#         warnings.simplefilter('ignore')
#         g1_fit = fit_generic_continuum(spectrum1, model = models.Legendre1D(degree = 5))
    
    
#     data[:,1] = y1/g1_fit(x1*u.angstrom)
    
#     plt.figure(figsize=(20,6))
#     #plt.plot(data[:, 0] + 0.5, data[:, 1]/max(data[:, 1]))
    
    