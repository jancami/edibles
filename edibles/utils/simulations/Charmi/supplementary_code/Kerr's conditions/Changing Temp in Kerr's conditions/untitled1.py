# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 13:35:15 2022

@author: Charmi Bhatt
"""

import numpy as np
import pandas as pd
import astropy.constants as const
import matplotlib.pyplot as plt
import timeit
import scipy.stats as ss
from edibles.utils.edibles_oracle import EdiblesOracle
from edibles.utils.edibles_spectrum import EdiblesSpectrum
import warnings
from astropy.modeling import models
from astropy import units as u
from specutils.spectra import Spectrum1D
from specutils.fitting import fit_generic_continuum
import seaborn as sns
from scipy.signal import argrelextrema

plt.figure(figsize=(15,8))

HD166937 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Heather's data\HD166937_avg_spectra.txt", sep=',')

plt.plot(((1/HD166937['Wavelength'])- 0.0001512079)*1e8, (1-(1-HD166937['Flux'])/(1-0.903)*0.1) + 0.08 ,  label = '$\lambda$6614')
#plt.xlim(-2.5, 2.5)

HD_5797 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Heather's data\5797_HD_166937.txt", sep=',')



plt.plot(((1/HD_5797['Wavelength'])- 0.0001725111)*1e8 + 0.4, (1-(1-HD_5797['Flux'])/(1-0.909)*0.1 ) + 0.04,  label = '$\lambda$5797')
#plt.xlim(-2.5, 2.5)

HD_6379 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Heather's data\6379_HD_166937.txt", sep=',')

plt.plot(((1/HD_6379['Wavelength'])- 0.0001567584)*1e8 - 0.3, 1-(1-HD_6379['Flux'])/(1-0.915)*0.1 , label = '$\lambda$6379')
plt.xlim(-5,5)

print(min(HD_6379['Flux']))

# plt.xticks(fontsize=15) #, rotation=90)
# plt.yticks(fontsize=15) #, rotation=90)
plt.legend(title = ''r'$\zeta$', fontsize=15)
#plt.legend(fontsize = 15)
        
        
plt.gca().xaxis.set_major_locator(plt.MultipleLocator(1))
plt.gca().xaxis.set_minor_locator(plt.MultipleLocator(0.5))

plt.gca().yaxis.set_minor_locator(plt.MultipleLocator(0.01))
plt.xlabel('Wavenumber (in cm$^{-1}$)', fontsize = 15)
plt.ylabel('Normalized Intensity', fontsize = 15)



plt.legend()