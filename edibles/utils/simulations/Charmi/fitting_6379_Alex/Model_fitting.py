# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 11:32:49 2023

@author: alexr
"""

from functions import obs_curve_to_plot, obs_curve_to_fit, allowed_perperndicular_transitions,\
    get_rotational_spectrum, model_curve_to_fit, fit_model, curve_to_fit_wavenos, allowed_parallel_transitions
import matplotlib.pyplot as plt
import numpy as np

Jmax = 600

# curve_to_fit_wavenos()
#%%
sightline = '185859'


B = 0.0016 # cm^-1
delta_B = 0.6 # %
zeta = -0.55 # cm^-1
T = 50 # K
sigma = 0.1953  
origin = 0.05
# combinations  = allowed_perperndicular_transitions(Jmax)
combinations  = allowed_parallel_transitions(Jmax)


result = fit_model(B, delta_B, zeta, T, sigma, origin, combinations, sightline, transition = 'parallel', Jmax = Jmax)


#%%
Obs_data, x_equal_spacing, y_data_fit, std_dev = obs_curve_to_fit(sightline)

def residual(y_obs_data, y_model_data, sigma):
    y_obs_data = np.array(y_obs_data)
    y_model_data = np.array(y_model_data)
    r = (y_obs_data - y_model_data)/sigma
    return r

r = residual(y_data_fit, result.best_fit, std_dev)
print(r)

print(len(r))
print(type(r))

#%%



combinations  = allowed_perperndicular_transitions(300)

plt.figure(1).add_axes((0,0.2,0.6,0.5))
plt.scatter(x_equal_spacing, y_data_fit, label='Observations')
plt.plot(x_equal_spacing, result.best_fit, 'r-', label='Best Fit')
plt.legend()
plt.ylabel('Flux')

plt.figure(1).add_axes((0,0,0.6,0.2)) #residual plot
plt.plot(x_equal_spacing, r, linestyle = ':')
plt.ylabel('Residual')
plt.xlabel('Wavenumber / cm$^{-1}$')
plt.show()
# 
