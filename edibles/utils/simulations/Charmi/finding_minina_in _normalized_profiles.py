#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 11:46:04 2023

@author: charmibhatt
"""
import numpy as np
import pandas as pd
import astropy.constants as const
import matplotlib.pyplot as plt
import timeit
import scipy
import scipy.stats as ss
from scipy.signal import argrelextrema


Obs_data = pd.read_csv(r"/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/Heather's_data/6614_HD_166937.txt", sep = ',')
print(Obs_data.shape)
y_obs_data =  np.array(Obs_data['Flux'])
x_obs_data = np.array(Obs_data['Wavelength'])

# plt.Figure(figsize=(15,8))
# plt.plot(x_obs_data, y_obs_data)

# minima = [argrelextrema(y_obs_data, np.less)]
# minima_ind = minima[0]

# minimas = np.array([minima_ind, x_obs_data[minima_ind]], y_obs_data[minima_ind]).transpose()
# print(minimas)

Obs_data = Obs_data[Obs_data['Wavelength'] >= 6615.5]

Obs_data = Obs_data.drop(columns= ['index'], axis = 1)

plt.plot(Obs_data.iloc[:,0], Obs_data.iloc[:,1])

print(Obs_data)
print(np.mean(Obs_data['Flux']))
std_dev = np.sqrt((np.sum((Obs_data['Flux'] - np.mean(Obs_data['Flux']))**2))/(len(Obs_data['Flux']) - 1))
print(std_dev)
print(np.std(Obs_data, axis = 0))