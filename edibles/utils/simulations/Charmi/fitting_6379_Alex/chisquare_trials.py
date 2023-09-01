# -*- coding: utf-8 -*-
"""
Created on Fri Aug  4 10:37:16 2023

@author: alexr
"""

import scipy.stats as sc
import functions as fn
import numpy as np

sightline = '185859'
Jmax = 600
B = 0.0013
delta_B = 0.6
zeta = -0.97
T = 50
sigma = 0.20
origin = -0.08

Obs_data, xs, y_data, std = fn.obs_curve_to_fit(sightline)
# print(Obs_data)
# print(len(xs), len(y_data))
# print(std)


combinations = fn.allowed_parallel_transitions(Jmax)
x_mod, y_mod = fn.get_rotational_spectrum(B, delta_B, zeta, T, sigma, origin, combinations, transition = 'parallel', bell = False)
y_mod_interp = np.interp(xs, x_mod, y_mod)

def chi_sq(model_ys, obs_ys, std):
    '''
    Function to calcualate redued chi-squared value of model compared to observations. 
    
    model_ys and obs_ys must be of the same dimensions, this is checked and will raise an AssertionError if not.

    Parameters
    ----------
    model_ys : array-like
        Generated model fluxes, after having been passed through the np.interp process. 
    obs_ys : array-like
        Observed data.
    std : float
        Standard deviation of observations.

    Returns
    -------
    red_chi_sq : float
        Reduced chi-square value.

    '''
    assert len(model_ys) == len(obs_ys), 'Ensure the model and observed data are of the same dimension!'
    
    chi_sq = 0
    for i in range(len(obs_ys)):
        difference = (model_ys[i] - obs_ys[i])**2
        chi_sq += difference
    
    dof  = len(obs_ys) - 6
    chi_sq = chi_sq/(std**2)
    red_chi_sq = chi_sq/dof
    return red_chi_sq

print(chi_sq([12,8]*50, y_data, std))

scipy = sc.chisquare(y_data, y_mod_interp, 94)
print(scipy)






