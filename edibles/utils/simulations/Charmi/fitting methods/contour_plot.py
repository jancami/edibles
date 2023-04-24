#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 10:12:06 2023

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
from scipy import interpolate
import scipy.ndimage.filters as filters
import scipy.ndimage.morphology as morphology
from scipy.interpolate import griddata


#grid of original data

Bmin = 0.0028
Bmax = 0.005
stepsize_B = 0.0001
Bs  = np.arange(Bmin, Bmax, stepsize_B)
print(Bs)
print(len(Bs))

Tmin = 40
Tmax = 70
stepsize_T= 1
Ts =  np.arange(Tmin, Tmax, stepsize_T)
print(Ts)
print(len(Ts))

BB, TT = np.meshgrid(Bs, Ts)
points = np.empty((660, 2))
points[:, 0] = BB.flatten()
points[:, 1] = TT.flatten()


plt.figure(figsize = (15,8))

BB = pd.read_csv(r'/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/fitting methods/BB_Bmin_0.0028_Bmax_0.003200000000000001_stepsize_1e-05_.txt', delim_whitespace=(True), header = None)
TT = pd.read_csv(r'/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/fitting methods/Bmax_0.0032_TT_Tmin_57.0_Tmax_67.5_stepsize_0.5_.txt', delim_whitespace=(True), header = None)
red_chi = pd.read_csv(r'/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/fitting methods/red_chi_finer_B_0.0028_to_0.0032_57_to_68.txt', delim_whitespace=(True), header = None)

#correct origin 166937
BB = pd.read_csv(r'/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/fitting methods/166937_BB_Bmin_0.0028_Bmax_0.004899999999999996_stepsize_0.0001_.txt', delim_whitespace=(True), header = None)
TT = pd.read_csv(r'/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/fitting methods/166937_Bmax_0.005_TT_Tmin_40_Tmax_69_stepsize_1_.txt', delim_whitespace=(True), header = None)
red_chi = pd.read_csv(r'/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/fitting methods/166937_red_chi_Bmin_0.0028_Bmax_0.004899999999999996_Tmin_40_Tmax_69_.txt', delim_whitespace=(True), header = None)


print(BB.shape)
print(TT.shape)
print(red_chi.shape)

chi = 179 * red_chi


#plt.contourf(BB, TT, red_chi, 20, c=chi, cmap='Oranges_r')
# cbar = plt.colorbar()
# cbar.ax.tick_params(labelsize=15 )  # adjust colorbar tick label size
# cbar.set_label(r'Reduced $\chi^{2}$', fontsize=20, labelpad = 15)  # set colorbar title
# plt.xlabel('B', size=15)
# plt.ylabel('T', size=15)
# plt.xticks(size=15)
# plt.yticks(size=14)

#print(min(red_chi.min()))

min_idx = np.argmin(red_chi)
min_x, min_y = np.unravel_index(min_idx, red_chi.shape)
b = (BB.iloc[1,:])
t = (TT.iloc [:,1])
print('Global minimum point:', (b[min_x], t[min_y]))
print('Global minimum value:', red_chi[min_x][ min_y])
plt.scatter(b[min_x], t[min_y], marker='o', c = 'black')



#interpolating and finding 1 sigma uncertainty

Bmin = 0.0028
Bmax = 0.0048
stepsize_B = 0.000001
Bs  = np.arange(Bmin, Bmax, stepsize_B)
# print(Bs)
print(len(Bs))

Tmin = 40
Tmax = 69
stepsize_T= 0.1
Ts =  np.arange(Tmin, Tmax, stepsize_T)
# print(Ts)
print(len(Ts))
BBM, TTM = np.meshgrid(Bs, Ts)
print(BBM.shape)
print(TTM.shape)
#print(red_chinew.shape)
red_chi = np.array(red_chi)
red_chinew = griddata(points, red_chi.flatten(), (BBM, TTM), method='linear')
print(red_chinew.shape)

chi = 179* red_chinew
print(chi)


plt.contourf(BBM, TTM, chi, 20, c=chi, cmap='Oranges_r')
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=15 )  # adjust colorbar tick label size
cbar.set_label(r' $\chi^{2}$', fontsize=20, labelpad = 15)  # set colorbar title
plt.xlabel('B', size=15)
plt.ylabel('T', size=15)
plt.xticks(size=15)
plt.yticks(size=14)

print(chi.min() + 1)

chi[chi >= chi.min() + 1] = None
chi_local = chi
indi = np.argwhere(~np.isnan(chi_local))


# for i in indi:
#     ii = i[0]
#     jj = i[1]
#     BBM[ii][jj] = None
#     TTM[ii][jj] = None
    
# print(BBM)
    
plt.contourf(BBM, TTM, chi_local, 10, c=chi_local, cmap='Greens_r')
# cbar = plt.colorbar()
# cbar.ax.tick_params(labelsize=15 )  # adjust colorbar tick label size
# cbar.set_label(r'Reduced $\chi^{2}$', fontsize=20, labelpad = 15)  # set colorbar title
# plt.xlabel('B', size=15)
# plt.ylabel('T', size=15)
# plt.xticks(size=15)
# plt.yticks(size=14)

plt.vlines(b[min_x], ymin = 40, ymax = t[min_y], color = 'white')
plt.vlines(0.00418, ymin = 40, ymax = 48, color = 'white', linestyles = 'dotted')
plt.vlines(0.00373, ymin = 40, ymax = 54.5, color = 'white', linestyles = 'dotted')

plt.hlines(t[min_y], xmin = 0.0028, xmax = b[min_x], color = 'white')
plt.hlines(48, xmin = 0.0028, xmax = 0.00418, color = 'white', linestyles = 'dotted')
plt.hlines(54.5, xmin = 0.0028, xmax = 0.00373, color = 'white', linestyles = 'dotted')

plt.scatter(b[min_x], t[min_y], zorder = 3, marker='o', c = 'black')

print(b[min_x])
print(t[min_y])
