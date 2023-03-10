#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 18:24:12 2023

@author: charmibhatt
"""

import numpy as np
import pandas as pd
import astropy.constants as const
import matplotlib.pyplot as plt
import timeit
import scipy.stats as ss
from scipy.signal import argrelextrema

# chi_plane_cordinates = pd.read_csv(r'/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/chi_plane_cordinates.txt', delim_whitespace=(True))
# #ax.scatter3D(chi_plane_cordinates.iloc[:,0], chi_plane_cordinates.iloc[:,1], chi_plane_cordinates.iloc[:,2]) #, c = chi_plane_cordinates.iloc[:,2]) #rstride=1, cstride=1
# chi_plane_cordinates = chi_plane_cordinates.sort_values(by=['chi2'])
# print(chi_plane_cordinates.to_string())


BB = pd.read_csv(r'/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/fitting methods/BB_Bmin_0.0005_Bmax_0.0485_stepsize_0.002_.txt', delim_whitespace=(True), header = None)
TT = pd.read_csv(r'/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/fitting methods/Bmax_0.05_TT_Tmin_5_Tmax_95_stepsize_5_.txt', delim_whitespace=(True), header = None)
red_chi = pd.read_csv(r'/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/fitting methods/red_chi_coarse.txt', delim_whitespace=(True), header = None)

#upto 0.01
BB = pd.read_csv(r'/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/fitting methods/BB_Bmin_0.0005_Bmax_0.009500000000000001_stepsize_0.001_.txt', delim_whitespace=(True), header = None)
TT = pd.read_csv(r'/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/fitting methods/TT_Tmin_5_Tmax_95_stepsize_5_.txt', delim_whitespace=(True), header = None)
red_chi = pd.read_csv(r'/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/fitting methods/red_chi_coarse_Bmax_0.01.txt', delim_whitespace=(True), header = None)

print(BB.shape)
print(TT.shape)
print(red_chi.shape)

#print(red_chi)
#BB = BB.pop(BB.columns[0])

print((BB))
print(BB.shape)

fig = plt.figure(figsize=(12,10))
ax = plt.axes(projection='3d')
ax.tick_params(axis='both', which='major', labelsize=20)
ax.tick_params(axis='both', which='minor', labelsize=8)
ax.plot_surface(BB, TT, red_chi, rstride=1, cstride=1, cmap='summer', edgecolor='none', alpha = 0.7)
ax.view_init(20, 100)

ax.set_xlabel('B', size = 20, labelpad = 20)
ax.set_ylabel('T', size = 20, labelpad = 20)
ax.set_zlabel(r'$\chi^{2}$',  size = 20, labelpad = 20, rotation   = 90)


import scipy.ndimage.filters as filters
import scipy.ndimage.morphology as morphology

def detect_local_minima(arr):
    # https://stackoverflow.com/questions/3684484/peak-detection-in-a-2d-array/3689710#3689710
    """
    Takes an array and detects the troughs using the local maximum filter.
    Returns a boolean mask of the troughs (i.e. 1 when
    the pixel's value is the neighborhood maximum, 0 otherwise)
    """
    # define an connected neighborhood
    # http://www.scipy.org/doc/api_docs/SciPy.ndimage.morphology.html#generate_binary_structure
    neighborhood = morphology.generate_binary_structure(len(arr.shape),2)
    # apply the local minimum filter; all locations of minimum value 
    # in their neighborhood are set to 1
    # http://www.scipy.org/doc/api_docs/SciPy.ndimage.filters.html#minimum_filter
    local_min = (filters.minimum_filter(arr, footprint=neighborhood)==arr)
    # local_min is a mask that contains the peaks we are 
    # looking for, but also the background.
    # In order to isolate the peaks we must remove the background from the mask.
    # 
    # we create the mask of the background
    background = (arr==0)
    # 
    # a little technicality: we must erode the background in order to 
    # successfully subtract it from local_min, otherwise a line will 
    # appear along the background border (artifact of the local minimum filter)
    # http://www.scipy.org/doc/api_docs/SciPy.ndimage.morphology.html#binary_erosion
    eroded_background = morphology.binary_erosion(
        background, structure=neighborhood, border_value=1)
    # 
    # we obtain the final mask, containing only peaks, 
    # by removing the background from the local_min mask
    detected_minima = local_min ^ eroded_background
    return np.where(detected_minima)  

    
local_minima_locations = detect_local_minima(red_chi)
print(local_minima_locations)

lml_i = local_minima_locations[0]
lml_j= local_minima_locations[1]

print('------------')
for i,j in zip(lml_i, lml_j):
    ax.scatter3D(BB.iloc[i][j], TT.iloc[i][j], red_chi.iloc[i][j], marker = 'o', c = 'black')
    print(BB.iloc[i][j])
    print(TT.iloc[i][j])
    print(red_chi.iloc[i][j])
    print('------------')









#ax.contour3D(BB, TT, red_chi, cmap='binary')

#ax.plot_surface(BB, TT, red_chi, rstride=1, cstride=1, cmap='viridis', edgecolor='none', alpha = 0.5)

# ax.contour(BB, TT, red_chi, zdir='z', offset=500, cmap='coolwarm')
# ax.contour(BB, TT, red_chi, zdir='z', offset=400, cmap='coolwarm')
# ax.contour(BB, TT, red_chi, zdir='z', offset=40, cmap='coolwarm')