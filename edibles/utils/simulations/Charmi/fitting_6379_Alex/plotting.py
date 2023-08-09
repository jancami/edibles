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

sightlines = ['24398','144470','147165','147683','149757','166937',
              '170740','184915','185418','185859','203532']

# print(len(sightlines))


#%% Plotting observational data of DIBs

fig, ax = plt.subplots()
xs, ys = [], []
std = 0
for star in sightlines: 
    xs, ys, std, x_axis = functions.obs_curve_to_plot(star, wavenos = 1, scaled = False)
    ax.plot(xs, ys, label = 'HD{}'.format(star))

ax.xaxis.set_major_locator(plt.MultipleLocator(1))
ax.xaxis.set_minor_locator(plt.MultipleLocator(0.5))
ax.set_xlabel(x_axis)
ax.set_ylabel('Flux')
ax.legend()
# plt.savefig('6379_obs_waveno.png')

# plot_colortable(mcolors.TABLEAU_COLORS, ncols=2, sort_colors=False)

#%% Simulations

Jmax = 600

combinations = functions.allowed_parallel_transitions(Jmax)
# print(combinations)

B = 0.0016 # cm^-1
delta_B = 0.6 # %
zeta = -0.55 # cm^-1
T = 50 # K
sigma = 0.1953  
origin = 0
xs, ys = functions.get_rotational_spectrum(B, delta_B, zeta, T, sigma, origin, combinations, Jmax = Jmax, transition = 'parallel')


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

x_obs, y_obs, std, x_label = functions.obs_curve_to_plot(sightline)

fig, ax = plt.subplots()
ax.plot(x_obs, y_obs, label = 'HD{}'.format(sightline))
ax.plot(xs, ys, label = 'T = {}K'.format(str(T)))
# ax.plot(xs1, ys1, label = 'T = 100K')
ax.set_title('B = {}, $\Delta B =${}, $\zeta = {}$,\n $\sigma$ = {}, origin = {}, Jmax = {}'.format(str(B), str(delta_B), str(zeta), str(sigma), str(origin), str(Jmax)))
ax.xaxis.set_major_locator(plt.MultipleLocator(1))
ax.xaxis.set_minor_locator(plt.MultipleLocator(0.5))
ax.set_xlabel(x_label)
ax.set_ylabel('Flux')
ax.legend()

plt.show()
# plt.savefig('B={}_DeltaB={}_zeta={}_T={}_sigma={}origin={}_HD{}.png'.format(str(B), str(delta_B), str(zeta), str(T), str(sigma), str(origin),sightline), bbox_inches = 'tight', format = 'png')

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
