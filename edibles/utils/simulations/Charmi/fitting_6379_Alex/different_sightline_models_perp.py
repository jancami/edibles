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
sightline = sightlines[2]
print(sightline)

#%% Perpendicular best fit parameters

Jmax = 800

combinations = fn.allowed_perperndicular_transitions(Jmax)

# Set 1:

# B = 0.0016
# delta_B = 0.35
# zeta = -0.4
# T = 70
# sigma = 0.1953
# origin = 0.22

# Set 2:
    
# B = 0.0018
# delta_B = 0.72
# zeta = -1.1
# T = 28.5
# sigma = 0.2
# origin = 0.16

# # Set 3:
    
# B = 0.005
# delta_B = 0.6
# zeta = -0.35
# T = 30
# sigma = 0.1953
# origin = 0.22

# # Set 4

# B = 0.0010 # cm^-1
# delta_B = 0.72 # %
# zeta = -1.1 # cm^-1
# T = 160 # K
# sigma = 0.2 
# origin = 0.08

# # Set 5 - best fit from set 2 on HD185859

B = 0.0014747
delta_B = 0.1460773
zeta = 0.188564
# T = 299.91933
T = 80
sigma = 0.1934417
# origin = -0.258784
origin = -0.4

#%% Getting data

#Observational data:
Obs_data, xs, y_data, std = fn.obs_curve_to_fit(sightline)

#Computing model:
result = fn.fit_model(B, delta_B, zeta, T, sigma, origin, combinations, sightline, transition = 'perpendicular', Jmax = Jmax)

