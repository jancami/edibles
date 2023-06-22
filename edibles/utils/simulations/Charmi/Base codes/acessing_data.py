# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 13:54:26 2022

@author: Charmi Bhatt
"""

import numpy as np
import matplotlib.pyplot as plt
from edibles.utils.edibles_oracle import EdiblesOracle
from edibles.utils.edibles_spectrum import EdiblesSpectrum
import warnings
from astropy.modeling import models
from astropy import units as u
from specutils.spectra import Spectrum1D
from specutils.fitting import fit_generic_continuum
from astropy import constants as const



'''works:'''

starName = 'HD 185418'
#put lower range of wavelengths to extract from edibles data
minWave = 6612

#put upper range of wavelengths to extract from edibles data
maxWave = 6615

pythia = EdiblesOracle()
rawList = pythia.getFilteredObsList(object = [starName], MergedOnly = True, WaveMin = minWave, WaveMax = maxWave)
fnames = rawList.tolist()
obs = len(fnames)


sp = EdiblesSpectrum(fnames[0])

# print(sp.target)
# print(sp.datetime.date())
# print(sp.raw_wave)

    
sp.getSpectrum(xmin = max(minWave, np.min(sp.raw_wave)+1)
                    , xmax = min(maxWave, np.max(sp.raw_wave)-1))

                   
data = np.array([sp.bary_wave, sp.bary_flux]).transpose()

wave_data = data[:,0]
int_data = data[:,1]

c = 299792458  # in m/s speed of light
v = -6.5 #velocity_of_cloud

wave_data = wave_data*(1+ (v/c)) # doppler shift


plt.figure(figsize=(12,6))
plt.plot(wave_data, (int_data/max(int_data)), color = 'black', label = 'getSpectrum')


spectrum1 = Spectrum1D(flux = int_data*u.dimensionless_unscaled, spectral_axis = wave_data*u.angstrom)

g1_fit = fit_generic_continuum(spectrum1, model = models.Legendre1D(degree = 1))

int_data = int_data/g1_fit(wave_data*u.angstrom)

plt.plot(wave_data, int_data/max(int_data), label='fit_generic_continuum')
plt.legend()

print(len(wave_data))
print(wave_data)
print(len(int_data))



'''Previously Used Code'''

# leftEdge = 0
# rightEdge = 0
    
# if minWave <= np.min(sp.raw_wave):
#     leftEdge = 1
# if maxWave >= np.max(sp.raw_wave):
#     rightEdge = 1

# data = np.delete(np.array([sp.bary_wave, sp.bary_flux]).transpose(), 
#                             np.logical_or(sp.bary_wave <= np.min(sp.bary_wave) + 40.0*leftEdge, 
#                                           sp.bary_wave >= np.max(sp.bary_wave) - 40.0*rightEdge), 0)




# spectrum1 = Spectrum1D(flux = int_data*u.dimensionless_unscaled, spectral_axis = wave_data*u.angstrom)

# with warnings.catch_warnings():  # Ignore warnings
#     warnings.simplefilter('ignore')
#     g1_fit = fit_generic_continuum(spectrum1, model = models.Legendre1D(degree = 5))


# int_data = int_data/g1_fit(wave_data*u.angstrom)
