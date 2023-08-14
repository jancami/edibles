# -*- coding: utf-8 -*-
"""
Created on Sat Aug 12 19:12:45 2023

@author: alexr

Making models for different sightlines based on certain initial parameters
"""

import functions as fn
import numpy as np

sightlines = ['24398','144470','147165','147683','149757','166937',
              '170740','184915','185418','203532','185859']
sightline = sightlines[-1]
print(sightline)

#%% Perpendicular best fit parameters

Jmax = 600

combinations = fn.allowed_parallel_transitions(Jmax)

# Set 1:

# B = 0.0016
# delta_B = 0.6
# zeta = -0.55
# T = 70
# sigma = 0.1953
# origin = 0

# Set 2:
    
B = 0.0013
delta_B = 0.59
zeta = -0.97
T = 50
sigma = 0.2
origin = -0.08

# # Set 3:
    
# B = 0.005
# delta_B = 0.6
# zeta = -0.35
# T = 30
# sigma = 0.1953
# origin = 0.22

#%% Getting data

#Observational data:
Obs_data, xs, y_data, std = fn.obs_curve_to_fit(sightline)

#Computing model:
result = fn.fit_model(B, delta_B, zeta, T, sigma, origin, combinations, sightline, transition = 'parallel', Jmax = Jmax)

