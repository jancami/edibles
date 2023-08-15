# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 12:34:14 2023

@author: alexr
"""

# from functions import obs_curve_to_plot, allowed_perperndicular_transitions,\
#     get_rotational_spectrum, obs_curve_to_fit
import functions
import matplotlib.pyplot as plt
import beepy as bp
import numpy as np
import pandas as pd

sightlines = ['24398','144470','147165','147683','149757','166937',
              '170740','184915','185418','185859','203532']

# print(len(sightlines))



#%% Plotting observational data of DIBs

fig, ax = plt.subplots(figsize = (9,6), facecolor = 'none')
xs, ys = [], []
std = 0

l = 35
h = 85
for star in sightlines: 
    xs, ys, std, x_axis = functions.obs_curve_to_plot(star, wavenos = 0, scaled = 0, zero = '6379')
    ax.plot(xs[l:h], ys[l:h], label = 'HD{}'.format(star), linewidth = 2)

ax.xaxis.set_major_locator(plt.MultipleLocator(0.5))
ax.xaxis.set_minor_locator(plt.MultipleLocator(0.1))
ax.set_xlabel(x_axis)
ax.set_ylabel('Flux')
# ax.set_facecolor('none')
ax.legend()

# plt.savefig("C:\\Users\\alexr\\OneDrive - Durham University\\GRI Mitacs Summer 23\\Project\\Presentation\\wavelength_plot_partial_transparent.png", bbox_inches = 'tight')
# plt.savefig('fig.png')

plt.show()
# plot_colortable(mcolors.TABLEAU_COLORS, ncols=2, sort_colors=False)

#%% Simulations

Jmax = 300

# combinations = functions.allowed_parallel_transitions(Jmax)
combinations = functions.allowed_perperndicular_transitions(Jmax)
# print(combinations)

B = 0.0016 # cm^-1
delta_B = -0.3 # %
zeta = -0.35 # cm^-1
T = 70 # K
sigma = 0.1953 
origin = 0.1
xs, ys = functions.get_rotational_spectrum(B, delta_B, zeta, T, sigma, origin, combinations, Jmax = Jmax, transition = 'perpendicular')


# B1 = 0.0016 # cm^-1
# delta_B1 = 0.6 # %
# zeta1 = -0.55 # cm^-1
# T1 = 100 # K
# sigma1 = 0.1953  
# origin1 = 0

# xs1, ys1 = functions.get_rotational_spectrum(B1, delta_B1, zeta1, T1, sigma1, origin1, combinations, Jmax = Jmax, transition = 'parallel')

# bp.beep(sound='wilhelm')
# plt.plot(xs, ys)
#%% Plot obs data and model

sightline = '185859'

x_obs, y_obs, std, x_label = functions.obs_curve_to_plot(sightline, wavenos=True, zero = 'min')

fig, ax = plt.subplots()#facecolor = 'none')
ax.plot(x_obs, y_obs, label = 'HD{}'.format(sightline))
ax.plot(xs, ys, label = 'Model')
# ax.plot(xs1, ys1, label = 'T = 100K')
ax.set_title('B = {}, $\Delta B =${}, $\zeta = {}$,\n $\sigma$ = {}, T = {}, origin = {}'.format(str(B), str(delta_B), str(zeta), str(sigma), str(T), str(origin)))
ax.xaxis.set_major_locator(plt.MultipleLocator(1))
ax.xaxis.set_minor_locator(plt.MultipleLocator(0.5))
ax.set_xlabel(x_label)
ax.set_ylabel('Flux')
ax.legend()

# plt.savefig('B={}_DeltaB={}_zeta={}_T={}_sigma={}origin={}_HD{}.png'.format(str(B), str(delta_B), str(zeta), str(T), str(sigma), str(origin),sightline), bbox_inches = 'tight', format = 'png')
# plt.savefig("C:\\Users\\alexr\\OneDrive - Durham University\\GRI Mitacs Summer 23\\Project\\Presentation\\Figures\\bluewing_demo.png", bbox_inches = 'tight')

plt.show()
#%% lmfit model plot creation function

sightline = '185859'


Obs_data, x_equal_spacing, y_data_fit, std_dev = functions.obs_curve_to_fit(sightline)

combinations  = functions.allowed_perperndicular_transitions(300)
x_model_data, y_model_data = functions.get_rotational_spectrum(B, delta_B, zeta, T, sigma, origin, combinations, Jmax = 300, bell = True)

plt.figure(1).add_axes((0,0.2,0.6,0.5))
plt.scatter(x_equal_spacing, y_data_fit, label='Observations')
plt.plot(x_model_data, y_model_data, 'r-', label='Best Fit')

plt.ylabel('y')

plt.figure(1).add_axes((0,0,0.6,0.2)) #residual plot
plt.legend()
plt.show()
plt.xlabel('x')



#%% Checking interpolations

sightline ='166937'

x_obs, y_obs, std1, x_label = functions.obs_curve_to_plot(sightline)
Obs_data, x_equal_spacing, y_interp, std2 = functions.obs_curve_to_fit(sightline)

assert std1 == std2

# print(len(x_obs))
l = 65
u = 82
x_obs_cut = x_obs[l:u]
y_obs_cut = y_obs[l:u]

# plt.plot(x_obs, y_obs, label = 'Non interpolated')
plt.plot(x_obs_cut, y_obs_cut, label = 'Non interpolated')
plt.plot(x_equal_spacing, y_interp, label = 'Interpolated')

plt.xlabel(x_label)
plt.ylabel('Flux')
plt.legend()

plt.show()

#%%

sightline = '166937'

x_obs, y_obs, std1, x_label = functions.obs_curve_to_plot(sightline, wavenos = False, scaled = False)

plt.plot(x_obs, y_obs)


plt.xlabel(x_label)
plt.ylabel('Flux')
# plt.legend()

plt.show()
#%% Plotting calculated model

B = 0.00147468
delta_B = 0.14607728
zeta = 0.18856398
T = 299.919328
sigma = 0.19344172
origin = -0.25878401
Jmax = 300

combinations = functions.allowed_perperndicular_transitions(300)
xs, ys = functions.get_rotational_spectrum(B, delta_B, zeta, T, sigma, origin, combinations, Jmax = Jmax, transition = 'perpendicular')

Obs_data, x_equal_spacing, y_data_fit, std_dev = functions.obs_curve_to_fit('185859')
# print(Obs_data)
# print(len(xs))
# print(len(ys))
# print(len(Obs_data))
# print(len(x_equal_spacing))
# print(len(y_data_fit))
#%%

obs = {'xs': xs, 'ys': ys}
df = pd.DataFrame(obs)
print(df)
df = df[df['ys']<=0.95]
print(df)
y_obs_interp = np.interp(x_equal_spacing, df['xs'],df['ys'])

#%%
def residual(y_obs_data, y_model_data, sigma):
    y_obs_data = np.array(y_obs_data)
    y_model_data = np.array(y_model_data)
    r = (y_obs_data - y_model_data)/sigma
    return r
r = residual(y_data_fit, y_obs_interp, std_dev)

# Plot best fit model, observation and residuals on same axis
def plot_best_fit(result, x_equal_spacing, y_obs_data):
    plt.figure(1, facecolor='none').add_axes((0,0.2,0.6,0.5))
    plt.title('B = {}, $\Delta B =${}, $\zeta = {}$,\n $\sigma$ = {}, T = {}, origin = {}'.format(str(B), str(delta_B), str(zeta), str(sigma), str(T), str(origin)))

    plt.scatter(x_equal_spacing, y_data_fit, label='Observations')
    plt.plot(x_equal_spacing, result, 'r-', label='Best Fit')
    plt.legend()
    plt.ylabel('Flux')

    plt.figure(1).add_axes((0,0,0.6,0.2)) #residual plot
    plt.plot(x_equal_spacing, r, linestyle = ':')
    plt.ylabel('Residual')
    plt.xlabel('Wavenumber / cm$^{-1}$')
    plt.savefig("C:\\Users\\alexr\\OneDrive - Durham University\\GRI Mitacs Summer 23\\Project\\Presentation\\Figures\\perp_model.png", bbox_inches = 'tight')

    plt.show()
plot_best_fit(y_obs_interp, x_equal_spacing, y_data_fit)

