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

plt.figure(figsize = (15,8))

BB = pd.read_csv(r'/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/fitting methods/BB_Bmin_0.0028_Bmax_0.003200000000000001_stepsize_1e-05_.txt', delim_whitespace=(True), header = None)
TT = pd.read_csv(r'/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/fitting methods/Bmax_0.0032_TT_Tmin_57.0_Tmax_67.5_stepsize_0.5_.txt', delim_whitespace=(True), header = None)
red_chi = pd.read_csv(r'/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/fitting methods/red_chi_finer_B_0.0028_to_0.0032_57_to_68.txt', delim_whitespace=(True), header = None)
print(BB.shape)
print(TT.shape)
print(red_chi.shape)

chi = 179 * red_chi

#plt.imshow(red_chi, extent=[0.0028, 0.0032, 57, 68], origin='lower',
          # cmap='RdGy', alpha=0.5)
# plt.contourf(BB, TT, chi, 20, c=chi, cmap='Oranges')

# plt.xlim(0.00275, 0.00325)
# plt.ylim(56.5, 68)


# plt.xlabel('B', size = 20)
# plt.ylabel('T', size = 20)
# plt.xticks(size = 20)
# plt.yticks(size = 20)

# #plt.colorbar()
# plt.colorbar().set_label((r'$\chi^{2}$'), fontsize = 10)

plt.contourf(BB, TT, red_chi, 20, c=chi, cmap='Oranges_r')
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=15 )  # adjust colorbar tick label size
cbar.set_label(r'$\chi^{2}$', fontsize=20, labelpad = 15)  # set colorbar title

# plt.xlim(0.00275, 0.00325)
# plt.ylim(56.5, 68)

plt.xlabel('B', size=15)
plt.ylabel('T', size=15)
plt.xticks(size=15)
plt.yticks(size=14)

print(min(red_chi.min()))

