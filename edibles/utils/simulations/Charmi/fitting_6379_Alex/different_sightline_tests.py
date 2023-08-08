# -*- coding: utf-8 -*-
"""
Created on Fri Aug  4 12:24:43 2023

@author: alexr
"""

#%% Parameters

import functions as fn
import numpy as np
import matplotlib.pyplot as plt

sightlines = ['24398','144470','147165','147683','149757','166937',
              '170740','184915','185418','203532']#,'185859']


sightline = sightlines[1]
# sightline = '185859'
Jmax = 600

#These parameters are given to 3sf however for their uncertainties they may be meant to be given to less

#Parallel best fit parameters:
B = 0.00129
delta_B = 0.594
zeta = -0.969
T = 49.7
sigma = 0.198
origin = -0.0793

# #Perpendicular best fit parameters:
# B = 0.00179
# delta_B = 0.717
# zeta = -1.10
# T = 28.6
# sigma = 0.206
# origin = -0.157

#%% Getting data

#Observational data:
Obs_data, xs, y_data, std = fn.obs_curve_to_fit(sightline)

#Computing model:
combinations = fn.allowed_parallel_transitions(Jmax)
x_mod, y_mod = fn.get_rotational_spectrum(B, delta_B, zeta, T, sigma, origin, combinations, transition = 'parallel', bell = False)
y_mod_interp = np.interp(xs, x_mod, y_mod)

#%% Chi squared value

chsq = fn.chi_sq(y_mod_interp, y_data, std)
print(chsq)
#%% Plotting

fig, ax = plt.subplots()#figsize = (15,10))
ax.plot(xs, y_data, label = 'HD{}'.format(sightline)) # Observational data plot
# for star in sightlines: # For plotting multiple stars 
#     Obs_data, x, y, std = fn.obs_curve_to_fit(star)
#     ax.plot(x, y, label = 'HD{}'.format(star), linewidth = 1.3)#, linestyle = '-.')
# Obs_data, x, y, std = fn.obs_curve_to_fit(sightline)
# ax.plot(x, y, label = 'HD{}'.format(sightline), linewidth = 3)
ax.plot(xs, y_mod_interp, label = 'HD185859 Model', linewidth = 4) # Model data plot
ax.set_title('Chi Squared = {}'.format(chsq))
ax.xaxis.set_major_locator(plt.MultipleLocator(1))
ax.xaxis.set_minor_locator(plt.MultipleLocator(0.5))
ax.set_xlabel('Wavenumber / $cm^{-1}$')
ax.set_ylabel('Flux')
ax.legend()

plt.show()
