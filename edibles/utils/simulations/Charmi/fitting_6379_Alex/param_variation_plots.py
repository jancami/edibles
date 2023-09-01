# -*- coding: utf-8 -*-
"""
Created on Mon Aug 28 10:43:48 2023

@author: alexr
"""

import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt#
import functions as fn
import numpy as np



path = Path(r'C:/Users/alexr/OneDrive - Durham University/GRI Mitacs Summer 23/Project')
file = 'bluewing_results.csv'

params_df = pd.read_csv(path / file, sep = ',')
params_df = params_df.iloc[:4]

print(params_df)

sightlines = params_df['Sightline']
chisqs = params_df['chisq']
Bs = params_df['B']
B_errs = params_df['B_err']
delta_Bs = params_df['delta_B']
delta_B_errs = params_df['delta_B_err']
Ts = params_df['T']
T_errs = params_df['T_err']
sigmas = params_df['sigma']
sigma_errs = params_df['sigma_err']
origins = params_df['origin']
origin_errs = params_df['origin_err']

def residual(y_obs_data, y_model_data, sigma):
    y_obs_data = np.array(y_obs_data)
    y_model_data = np.array(y_model_data)
    r = (y_obs_data - y_model_data)/sigma
    return r

x_indices = range(len(sightlines))

#%%
# plt.plot()
print(Bs)
print(B_errs)
fig, ax = plt.subplots(facecolor = 'none')
ax.errorbar(x_indices, Bs, yerr = B_errs, xerr = None, 
             ecolor = 'gray', elinewidth = 2, fmt = 'o', capsize = 3, 
             ms = 4, mec = 'black',mfc = 'black')
ax.set_xticks(x_indices, sightlines)
ax.set_xticklabels(sightlines, rotation = 45)
ax.set_ylim((-0.001,0.006))

ax.set_xlabel('Sightlines', fontsize = 'x-large')
ax.set_ylabel('B / $cm^{-1}$', fontsize = 'x-large')

plt.savefig("C:\\Users\\alexr\\OneDrive - Durham University\\GRI Mitacs Summer 23\\Project\\Report\\Figure\\B_var.png", bbox_inches = 'tight')

plt.show()

#%% T plot

fig, ax = plt.subplots(facecolor = 'none')
ax.errorbar(x_indices, Ts, yerr = T_errs, xerr = None, 
             ecolor = 'gray', elinewidth = 2, fmt = 'o', capsize = 3, 
             ms = 4, mec = 'black',mfc = 'black')
ax.set_xticks(x_indices, sightlines)
ax.set_xticklabels(sightlines, rotation = 45)
# ax.set_ylim((-0.001,0.006))

ax.set_xlabel('Sightlines', fontsize = 'x-large')
ax.set_ylabel('T / K', fontsize = 'x-large')

plt.savefig("C:\\Users\\alexr\\OneDrive - Durham University\\GRI Mitacs Summer 23\\Project\\Report\\Figure\\T_var.png", bbox_inches = 'tight')

plt.show()
#%% sigma plot

fig, ax = plt.subplots(facecolor = 'none')
ax.errorbar(x_indices, sigmas, yerr = sigma_errs, xerr = None, 
             ecolor = 'gray', elinewidth = 2, fmt = 'o', capsize = 3, 
             ms = 4, mec = 'black',mfc = 'black')
ax.set_xticks(x_indices, sightlines)
ax.set_xticklabels(sightlines, rotation = 45)
# ax.set_ylim((-0.001,0.006))

ax.set_xlabel('Sightlines', fontsize = 'x-large')
ax.set_ylabel('$\sigma$ / $cm^{-1}$', fontsize = 'x-large')

plt.savefig("C:\\Users\\alexr\\OneDrive - Durham University\\GRI Mitacs Summer 23\\Project\\Report\\Figure\\sigma_var.png", bbox_inches = 'tight')

plt.show()


#%% Results plot
sightline_data = params_df.iloc[3]
sightline = sightline_data['Sightline'][2:]
print(sightline)

B = sightline_data['B']
delta_B = sightline_data['delta_B']
zeta = sightline_data['zeta']
T = sightline_data['T']
sigma = sightline_data['sigma']
origin = sightline_data['origin']
Jmax = 800

combinations = fn.allowed_perperndicular_transitions(Jmax)

xs, ys = fn.get_rotational_spectrum(B, delta_B, zeta, T, sigma, origin, combinations, Jmax = Jmax, transition = 'perpendicular')
Obs_data, x_equal_spacing, y_data_fit, std_dev = fn.obs_curve_to_fit(sightline, fitrange='bluewing')
#%%

l = 40
u = 93
x_obs_cut = Obs_data['Wavelength'][l:u]
y_obs_cut = Obs_data['Flux'][l:u]

plt.plot(x_obs_cut, y_obs_cut, label = 'HD{} Observations'.format(sightline), color = 'dimgray')
plt.plot(xs, ys, label = 'Model', color = 'black')
plt.plot(x_equal_spacing, y_data_fit, label = 'Fitting range', linewidth = 2, color = 'r')
# ax.xaxis.set_major_locator(plt.MultipleLocator(1))
# ax.xaxis.set_minor_locator(plt.MultipleLocator(0.5))
plt.xlabel('Wavenumber/$cm^{-1}$')
plt.ylabel('Flux')

plt.legend(loc='lower left', fontsize = 'small')

plt.savefig("C:\\Users\\alexr\\OneDrive - Durham University\\GRI Mitacs Summer 23\\Project\\Report\\Figure\\HD{}_Model_plot.png".format(sightline), bbox_inches = 'tight')

plt.show()

