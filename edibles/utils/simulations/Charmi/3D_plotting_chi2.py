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
from matplotlib import cm

Obs_data = pd.read_csv(r"/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/Heather's_data/6614_HD_166937.txt", sep = ',')

y_obs_data =  np.array(Obs_data['Flux'])
x_obs_data = np.array(Obs_data['Wavelength'])


# chi_plane_cordinates = pd.read_csv(r'/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/chi_plane_cordinates.txt', delim_whitespace=(True))
# #ax.scatter3D(chi_plane_cordinates.iloc[:,0], chi_plane_cordinates.iloc[:,1], chi_plane_cordinates.iloc[:,2]) #, c = chi_plane_cordinates.iloc[:,2]) #rstride=1, cstride=1
# chi_plane_cordinates = chi_plane_cordinates.sort_values(by=['chi2'])
# print(chi_plane_cordinates.to_string())


BB = pd.read_csv(r'/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/fitting methods/BB_Bmin_0.0005_Bmax_0.0485_stepsize_0.002_.txt', delim_whitespace=(True), header = None)
TT = pd.read_csv(r'/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/fitting methods/Bmax_0.05_TT_Tmin_5_Tmax_95_stepsize_5_.txt', delim_whitespace=(True), header = None)
red_chi = pd.read_csv(r'/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/fitting methods/red_chi_coarse.txt', delim_whitespace=(True), header = None)

#coarse grid
#upto 0.01
BB = pd.read_csv(r'/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/fitting methods/BB_Bmin_0.0005_Bmax_0.009500000000000001_stepsize_0.001_.txt', delim_whitespace=(True), header = None)
TT = pd.read_csv(r'/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/fitting methods/TT_Tmin_5_Tmax_95_stepsize_5_.txt', delim_whitespace=(True), header = None)
red_chi = pd.read_csv(r'/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/fitting methods/red_chi_coarse_Bmax_0.01.txt', delim_whitespace=(True), header = None)

print((np.argmin(red_chi)))

def find_min_2d_array(arr):
    flat_arr = np.array(arr).flatten()  # flatten the 2D array into a 1D array
    min_index = np.argmin(flat_arr)  # find the index of the minimum value
    min_val = flat_arr[min_index]  # get the minimum value using the index
    return min_val

min_val = find_min_2d_array(red_chi)
print(min_val) 

#fine grid
BB = pd.read_csv(r'/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/fitting methods/BB_Bmin_0.002_Bmax_0.004899999999999995_stepsize_0.0001_.txt', delim_whitespace=(True), header = None)
TT = pd.read_csv(r'/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/fitting methods/Bmax_0.05_TT_Tmin_35_Tmax_89_stepsize_1_.txt', delim_whitespace=(True), header = None)
red_chi = pd.read_csv(r'/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/fitting methods/red_chi_fine_B_0.002_to_0.005.txt', delim_whitespace=(True), header = None)

#finer grid
BB = pd.read_csv(r'//Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/fitting methods/BB_Bmin_0.0028_Bmax_0.0034900000000000018_stepsize_1e-05_.txt', delim_whitespace=(True), header = None)
TT = pd.read_csv(r'/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/fitting methods/Bmax_0.05_TT_Tmin_60.0_Tmax_67.5_stepsize_0.5_.txt', delim_whitespace=(True), header = None)
red_chi = pd.read_csv(r'/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/fitting methods/red_chi_finer_B_0.0028_to_0.0035.txt', delim_whitespace=(True), header = None)

#correct origin 166937
BB = pd.read_csv(r'/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/fitting methods/166937_BB_Bmin_0.0028_Bmax_0.004899999999999996_stepsize_0.0001_.txt', delim_whitespace=(True), header = None)
TT = pd.read_csv(r'/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/fitting methods/166937_Bmax_0.005_TT_Tmin_40_Tmax_69_stepsize_1_.txt', delim_whitespace=(True), header = None)
red_chi = pd.read_csv(r'/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/fitting methods/166937_red_chi_Bmin_0.0028_Bmax_0.004899999999999996_Tmin_40_Tmax_69_.txt', delim_whitespace=(True), header = None)

min_val = find_min_2d_array(red_chi)
print(min_val) 


# print(BB.to_string())
# print(TT.to_string())
# print(red_chi.to_string())

#print(red_chi)
#BB = BB.pop(BB.columns[0])


chi = 179 * red_chi
fig = plt.figure(figsize=(15,10))
ax = plt.axes(projection='3d', computed_zorder=False)
ax.tick_params(axis='x', which='major', labelsize=15, rotation = 60)
ax.tick_params(axis='both', which='minor', labelsize=8)
ax.plot_surface(BB, TT, red_chi, rstride=1, cstride=1, cmap='summer', edgecolor='none', alpha = 1)
ax.view_init(20, 100)
#ax.set_xlim(0.0028,0.0032)

ax.set_xlabel('B', size = 20, labelpad = 50)
ax.set_ylabel('T', size = 20, labelpad = 30)
ax.set_zlabel(r'$\chi^{2}$',  size = 20, labelpad = 20, rotation   = 90, zorder = 1, alpha = 0.8)


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


#value = 214.26
chi = 179 * red_chi 

local_minima_locations = detect_local_minima(red_chi)
print(local_minima_locations)

lml_i = local_minima_locations[0]
lml_j= local_minima_locations[1]

#print(TT.iloc[26][11])

# i = 26
# j= 11
# ax.scatter3D(BB.iloc[i][j], TT.iloc[i][j], red_chi.iloc[i][j], marker = 'o', c = 'black', zorder = 2)

print('------------')
for i,j in zip(lml_i, lml_j):
    ax.scatter3D(BB.iloc[i][j], TT.iloc[i][j], red_chi.iloc[i][j], marker = 'o', c = 'black', zorder = 2)
    print(BB.iloc[i][j])
    #print(BB[i][j])
    print(TT.iloc[i][j])
    print(red_chi.iloc[i][j])
    print(chi.iloc[i][j])
    print('----///////--------')



'''chi2 + 1'''

chi = 179 * red_chi 
indices = []   
lower_value = 197.6
upper_value = 197.7

for i in range(len(chi)):
        for j in range(len(chi.iloc[i])):
            if lower_value < chi.iloc[i][j] < upper_value:
                indices.append((i, j))
                
local_chi = []
local_red_chi = []
local_BB = []
local_TT = []

for i in range(len(indices)):
    index = indices[i]
    index = indices[i]
    ii = index[0]
    jj = index[1]
    # print(ii)
    # print(jj)
    local_chi.append(chi[ii][jj])
    local_red_chi.append(red_chi[ii][jj])
    local_BB.append(BB[ii][jj])
    local_TT.append(TT[ii][jj])
    

# chi_plus_one_coordinates = np.array([local_BB, local_TT, local_chi, local_red_chi]).transpose()
# np.set_printoptions(suppress=True, precision=4)

# print((chi_plus_one_coordinates))
            
# print(red_chi.iloc[3:7][23:27])

# # Specify the portion to be printed
# start_row = 1
# end_row = 9
# start_col = 19
# end_col = 31

# np.set_printoptions(suppress=True, precision=6)
# # Loop through the specified portion and print each element
# for i in range(start_row, end_row):
#     for j in range(start_col, end_col):
#         print(chi.iloc[i][j], end=" ")
#     print()
    
# B_uncertainty_in_chi= np.array(chi.iloc[start_row:end_row, start_col:end_col])
# B_uncertainty_in_BB= np.array(BB.iloc[start_row:end_row, start_col:end_col])
# B_uncertainty_in_TT= np.array(TT.iloc[start_row:end_row, start_col:end_col])

# ax.plot_surface(B_uncertainty_in_BB, B_uncertainty_in_TT, B_uncertainty_in_chi,rstride=1, cstride=1, cmap='autumn', edgecolor='none', alpha = 1)

#print(B_uncertainty)

# Save the new array into a text file
# np.savetxt("B_uncertainty.txt", B_uncertainty)



print(red_chi)

value = chi.iloc[i][j] 

print(value)
value_plus_one = value + 1
print(value_plus_one)

arr = chi
lower_value = value
upper_value = value_plus_one

def modify_array_between_values(arr, lower_value, upper_value):
    # Find the indices of all elements that are between lower_value and upper_value
    indices = []
    for i in range(len(arr)):
        for j in range(len(arr.iloc[i])):
            if lower_value < arr.iloc[i][j] < upper_value:
                indices.append((i, j))
    
def modify_array_to_keep_shape(arr, indices):
# Create a new 2D array of None values with the same shape as the original array
    new_arr = [[None for j in range(len(arr.iloc[i]))] for i in range(len(arr))]

    # Copy the values from the original array to the new array at the specified indices
    for i, j in indices:
        new_arr[i][j] = arr[i][j]
    
    return new_arr
            

arr = chi
lower_value = value
upper_value = value_plus_one
print(modify_array_between_values(arr, lower_value, upper_value))



print(modify_array_to_keep_shape(chi, indices))




#ax.contour3D(BB, TT, red_chi, cmap='binary')

ax.plot_surface(BB, TT, red_chi, rstride=1, cstride=1, cmap='viridis', edgecolor='none', alpha = 0.5)

# ax.contour(BB, TT, red_chi, zdir='y', offset=50, cmap='coolwarm')
# ax.contour(BB, TT, red_chi, zdir='y', offset=40, cmap='coolwarm')
# ax.contour(BB, TT, red_chi, zdir='z', offset=40, cmap='coolwarm')