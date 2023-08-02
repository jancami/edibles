# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 14:55:49 2023

@author: alexr
"""

import astropy.constants as const
import numpy as np
from functions import allowed_perperndicular_transitions
import timeit
import numba as nb
import matplotlib.pyplot as plt

#%%

def get_rotational_spectrum(B, delta_B, zeta, T, sigma, origin, Jmax = 300, bell = True):
    
    '''
    Generates a model of the rotational spectrum of a molecule based on input parameters
    - Requires allowed_perperndicular_transitions to be defined.
    - Prints the time taken to compute the model
    
    Args (all floats): 
        B (cm^-1): Ground rotational constant of the molecule. Assume 3D oblate symmetric top A = B = 2C.
        delta_B (%): Change in B to excited state. Assume delta_B = delta_C.
        zeta (cm^-1): Coriolis constant.
        T (K): Rotational Temperature.
        sigma: std of Gaussian smoothing.
        origin: allows for a shift along the x-axis.
        Jmax: Default 300, for passing to allowed_perperndicular_transitions function.
        bell (bool): default True, sounds a notification bell once the computation is complete.
    
    Returns:
        x_model_data, y_model_data (numpy array): calculated wavenumber and flux values for the model
    '''
    
    startg = timeit.default_timer()

    combinations  = allowed_perperndicular_transitions(Jmax)
    # rotational constants in cm-1
    ground_B = B
    ground_C = ground_B / 2
    delta_C = delta_B
    excited_B = ground_B + ((delta_B / 100) * ground_B)
    excited_C = ground_C + ((delta_C / 100) * ground_C)
    
    ground_Js = combinations['ground_J']
    excited_Js = combinations['excited_J']
    ground_Ks = combinations['ground_K']
    excited_Ks = combinations['excited_K']

    linelist = combinations

    delta_J = linelist['excited_J'] - linelist['ground_J']
    delta_K = linelist['excited_K'] - linelist['ground_K']

    # Calculating Linelist

    ground_Es = [] 
    for J, K in zip(ground_Js, ground_Ks):
        ground_E = ground_B * J * (J + 1) + (ground_C - ground_B) * (K ** 2)
        ground_Es.append(ground_E)

    linelist['ground_Es'] = ground_Es # Es are in units of cm^-1

    excited_Es = []
    for J, K, del_K in zip(excited_Js, excited_Ks, delta_K):
        if del_K == -1:
            excited_E = excited_B * J * (J + 1) + (excited_C - excited_B) * (K ** 2) - (
                (-2 * excited_C * zeta)) * K + excited_C ** 2
        elif del_K == 1:
            excited_E = excited_B * J * (J + 1) + (excited_C - excited_B) * (K ** 2) + (
                (-2 * excited_C * zeta)) * K + excited_C ** 2

        excited_Es.append(excited_E)

    linelist['excited_Es'] = excited_Es

    wavenos = []
    for i in range(len(linelist.index)):
        wavenumber = origin + excited_Es[i] - ground_Es[i]
        wavenos.append(wavenumber)

    linelist['wavenos'] = wavenos

    #Calculating Honl-London factors to determine relative intensities

    HL_factors = []


    for J, K, delta_J, delta_K in zip(ground_Js, ground_Ks, delta_J, delta_K):
        if delta_J == -1 and delta_K == -1:
            HL_factor = ((J - 1 + K) * (J + K)) / (J * ((2 * J) + 1))
        elif delta_J == -1 and delta_K == 1:
            HL_factor = ((J - 1 - K) * (J - K)) / (J * ((2 * J) + 1))
        elif delta_J == 0 and delta_K == -1:
            HL_factor = (J + 1 - K) * (J + K) / (J * (J + 1))
        elif delta_J == 0 and delta_K == 1:
            HL_factor = (J + 1 + K) * (J - K) / (J * (J + 1))
        elif delta_J == 1 and delta_K == -1:
            HL_factor = (J + 2 - K) * (J + 1 - K) / ((J + 1) * ((2 * J) + 1))
        elif delta_J == 1 and delta_K == 1:
            HL_factor = (J + 2 + K) * (J + 1 + K) / ((J + 1) * ((2 * J) + 1))

        HL_factors.append(HL_factor)

    linelist['HL_factors'] = HL_factors

    # Calculate populations of each level with Boltzmann eqn

    BD_factors = []

    h = const.h.cgs.value
    c = const.c.to('cm/s').value
    k = const.k_B.cgs.value

    for J, K, E in zip(ground_Js, ground_Ks, ground_Es):
        if K == 0:
            boltzmann_equation = (2 * ((2 * J) + 1)) * (np.exp((-h * c * E) / (k * T)))
        else:
            boltzmann_equation = (1 * ((2 * J) + 1)) * (np.exp((-h * c * E) / (k * T)))

        BD_factors.append(boltzmann_equation)

    linelist['BD_factors'] = BD_factors

    intensities = []
    for i in range(len(linelist.index)):
        strength = (HL_factors[i] * BD_factors[i])
        intensities.append(strength)

    linelist['intensities'] = intensities
    
    # linelist_time = timeit.default_timer()
    
    # print('>>>> Time taken to compute full linelist '+ str(linelist_time - startg) +' sec')

    # Smoothening the linelist using a Gaussian smoothing with std sigma and the numba decorator
    
    smooth_wavenos = np.linspace(np.min(linelist['wavenos']) - 1, np.max(linelist['wavenos']) + 1, 1000)  # grid_size

    Wavenos_arr = np.array(linelist['wavenos'])
    Intenisty_arr = np.array(linelist['intensities'])

    @nb.njit(parallel=True)
    def calculate_smooth_intensities(wavenos, intensities, smooth_wavenos, sigma):
        smooth_intensities = np.zeros(smooth_wavenos.shape)

        for i in nb.prange(len(smooth_wavenos)):
            wavepoint = smooth_wavenos[i]
            w_int = np.exp(-(wavenos - wavepoint) ** 2 / (2 * sigma ** 2)) * intensities
            smooth_intensities[i] = np.sum(w_int)

        return smooth_intensities

    # call the numba function with input data
    smooth_intensities = calculate_smooth_intensities(Wavenos_arr, Intenisty_arr, smooth_wavenos, sigma)

    smooth_data = np.array([smooth_wavenos, smooth_intensities]).transpose()
    smooth_data = np.delete(smooth_data, np.where(smooth_data[:, 1] <= 0.001 * (max(smooth_data[:, 1]))), axis=0)

    simu_waveno = smooth_data[:, 0]
    simu_intenisty = 1 - 0.1 * (smooth_data[:, 1] / max(smooth_data[:, 1])) # Scale to between 0.9 and 1
    
    model_data = np.array([simu_waveno, simu_intenisty]).transpose()

    # for units in wavelength
    # simu_wavelength = (1/simu_waveno)*1e8
    # model_data = np.array([simu_wavelength, simu_intenisty]).transpose()
    # model_data = model_data[::-1]

    y_model_data = model_data[:,1]
    x_model_data = model_data[:,0]
    endg = timeit.default_timer()

    # print('>>>> Time taken to smooth profile  ' + str(endg - linelist_time) + '  sec')
    print('>>>> Time taken to simulate this profile  ' + str(endg - startg) + '  sec')
    # print(linelist['ground_Es'])
    # print('>>>> Time for H-L factors if statement ' + str(endif - startif) + ' sec')
    print('==========')
    if bell:
        print('\a')
        
    time = endg - startg
    return x_model_data, y_model_data, time

#%% Optimized fn

def get_rotational_spectrum_opt(B, delta_B, zeta, T, sigma, origin, Jmax = 300, bell = True):

    startg = timeit.default_timer()

    combinations  = allowed_perperndicular_transitions(Jmax)
    # rotational constants in cm-1
    ground_B = B
    ground_C = ground_B / 2
    delta_C = delta_B
    excited_B = ground_B + ((delta_B / 100) * ground_B)
    excited_C = ground_C + ((delta_C / 100) * ground_C)
    
    ground_Js = combinations['ground_J']
    excited_Js = combinations['excited_J']
    ground_Ks = combinations['ground_K']
    excited_Ks = combinations['excited_K']

    linelist = combinations

    delta_J = linelist['excited_J'] - linelist['ground_J']
    delta_K = linelist['excited_K'] - linelist['ground_K']

    # Calculating Linelist

    ground_Es = [] 
    for J, K in zip(ground_Js, ground_Ks):
        ground_E = ground_B * J * (J + 1) + (ground_C - ground_B) * (K ** 2)
        ground_Es.append(ground_E)

    linelist['ground_Es'] = ground_Es # Es are in units of cm^-1

    excited_Es = []
    for J, K, del_K in zip(excited_Js, excited_Ks, delta_K):
        if del_K == -1:
            excited_E = excited_B * J * (J + 1) + (excited_C - excited_B) * (K ** 2) - (
                (-2 * excited_C * zeta)) * K + excited_C ** 2
        elif del_K == 1:
            excited_E = excited_B * J * (J + 1) + (excited_C - excited_B) * (K ** 2) + (
                (-2 * excited_C * zeta)) * K + excited_C ** 2

        excited_Es.append(excited_E)

    linelist['excited_Es'] = excited_Es

    wavenos = []
    for i in range(len(linelist.index)):
        wavenumber = origin + excited_Es[i] - ground_Es[i]
        wavenos.append(wavenumber)

    linelist['wavenos'] = wavenos

    #Calculating Honl-London factors to determine relative intensities

    HL_factors = []


    for J, K, delta_J, delta_K in zip(ground_Js, ground_Ks, delta_J, delta_K):
        if delta_J == -1 and delta_K == -1:
            HL_factor = ((J - 1 + K) * (J + K)) / (J * ((2 * J) + 1))
        elif delta_J == -1 and delta_K == 1:
            HL_factor = ((J - 1 - K) * (J - K)) / (J * ((2 * J) + 1))
        elif delta_J == 0 and delta_K == -1:
            HL_factor = (J + 1 - K) * (J + K) / (J * (J + 1))
        elif delta_J == 0 and delta_K == 1:
            HL_factor = (J + 1 + K) * (J - K) / (J * (J + 1))
        elif delta_J == 1 and delta_K == -1:
            HL_factor = (J + 2 - K) * (J + 1 - K) / ((J + 1) * ((2 * J) + 1))
        elif delta_J == 1 and delta_K == 1:
            HL_factor = (J + 2 + K) * (J + 1 + K) / ((J + 1) * ((2 * J) + 1))

        HL_factors.append(HL_factor)

    linelist['HL_factors'] = HL_factors

    # Calculate populations of each level with Boltzmann eqn

    BD_factors = []

    h = const.h.cgs.value
    c = const.c.to('cm/s').value
    k = const.k_B.cgs.value

    for J, K, E in zip(ground_Js, ground_Ks, ground_Es):
        if E/T >= 7:
            boltzmann_equation = 0
        elif K == 0:
            boltzmann_equation = (2 * ((2 * J) + 1)) * (np.exp((-h * c * E) / (k * T)))
        else:
            boltzmann_equation = (1 * ((2 * J) + 1)) * (np.exp((-h * c * E) / (k * T)))

        BD_factors.append(boltzmann_equation)
    
    linelist['BD_factors'] = BD_factors

    intensities = []
    for i in range(len(linelist.index)):
        strength = (HL_factors[i] * BD_factors[i])
        intensities.append(strength)

    linelist['intensities'] = intensities
    
    # linelist_time = timeit.default_timer()
    
    # print('>>>> Time taken to compute full linelist '+ str(linelist_time - startg) +' sec')

    # Smoothening the linelist using a Gaussian smoothing with std sigma and the numba decorator
    
    smooth_wavenos = np.linspace(np.min(linelist['wavenos']) - 1, np.max(linelist['wavenos']) + 1, 1000)  # grid_size

    Wavenos_arr = np.array(linelist['wavenos'])
    Intenisty_arr = np.array(linelist['intensities'])

    @nb.njit(parallel=True)
    def calculate_smooth_intensities(wavenos, intensities, smooth_wavenos, sigma):
        smooth_intensities = np.zeros(smooth_wavenos.shape)

        for i in nb.prange(len(smooth_wavenos)):
            wavepoint = smooth_wavenos[i]
            w_int = np.exp(-(wavenos - wavepoint) ** 2 / (2 * sigma ** 2)) * intensities
            smooth_intensities[i] = np.sum(w_int)

        return smooth_intensities

    # call the numba function with input data
    smooth_intensities = calculate_smooth_intensities(Wavenos_arr, Intenisty_arr, smooth_wavenos, sigma)

    smooth_data = np.array([smooth_wavenos, smooth_intensities]).transpose()
    smooth_data = np.delete(smooth_data, np.where(smooth_data[:, 1] <= 0.001 * (max(smooth_data[:, 1]))), axis=0)

    simu_waveno = smooth_data[:, 0]
    simu_intenisty = 1 - 0.1 * (smooth_data[:, 1] / max(smooth_data[:, 1])) # Scale to between 0.9 and 1
    
    model_data = np.array([simu_waveno, simu_intenisty]).transpose()

    # for units in wavelength
    # simu_wavelength = (1/simu_waveno)*1e8
    # model_data = np.array([simu_wavelength, simu_intenisty]).transpose()
    # model_data = model_data[::-1]

    y_model_data = model_data[:,1]
    x_model_data = model_data[:,0]
    endg = timeit.default_timer()

    # print('>>>> Time taken to smooth profile  ' + str(endg - linelist_time) + '  sec')
    print('>>>> Time taken to simulate this profile  ' + str(endg - startg) + '  sec')
    # print(linelist['ground_Es'])
    # print('>>>> Time for H-L factors if statement ' + str(endif - startif) + ' sec')
    print('==========')
    if bell:
        print('\a')
        
    time = endg - startg
    return x_model_data, y_model_data, time

#%%

B = 0.0016
delta_B = 0.2
zeta = -0.55
T = 100
sigma = 0.1953
origin = 0.22

xs, ys, time = get_rotational_spectrum(B, delta_B, zeta, T, sigma, origin, bell = True)
xs_opt, ys_opt, time_opt = get_rotational_spectrum_opt(B, delta_B, zeta, T, sigma, origin, Jmax = 300, bell = True) 


#%% Plot
fig, ax = plt.subplots()
ax.plot(xs, ys, label = 'Original fn, {}s'.format(time))
ax.plot(xs_opt, ys_opt, label = 'Opt fn, {}s'.format(time_opt))
ax.xaxis.set_major_locator(plt.MultipleLocator(1))
ax.xaxis.set_minor_locator(plt.MultipleLocator(0.5))
ax.set_title('B = {}, $\Delta B =${}, $\zeta = {}$,\n T = {}, $\sigma$ = {}, origin = {}'.format(str(B), str(delta_B), str(zeta), str(T,), str(sigma), str(origin)))
ax.set_xlabel('Wavenumber / cm$^{-1}$')
ax.set_ylabel('Flux')
ax.legend()


#%%

    
    # ground_Es = [] 
    # for J, K in zip(ground_Js, ground_Ks):
    #     ground_E = ground_B * J * (J + 1) + (ground_C - ground_B) * (K ** 2)
    #     ground_Es.append(ground_E)

    # linelist['ground_Es'] = ground_Es # Es are in units of cm^-1
    
    # end2_1 = timeit.default_timer()
    # print('Time to calculate ground Es '+str(end2_1-start2))
    
    # start2_2 = timeit.default_timer()
    # # @nb.njit()
    # excited_Es = []
    # for J, K, del_K in zip(excited_Js, excited_Ks, delta_K):
    #     E_const = excited_B * J * (J + 1) + (excited_C - excited_B) * (K ** 2) + excited_C ** 2
    #     if del_K == -1:
    #         excited_E = E_const - ((-2 * excited_C * zeta)) * K 
    #     elif del_K == +1:
    #         excited_E = E_const + ((-2 * excited_C * zeta)) * K 

    #     excited_Es.append(excited_E)

    # end2_2 = timeit.default_timer()

    # print('Time to calculate excited Es '+str(end2_2 - start2_2))

    # excited_Es = []
    # for J, K, del_K in zip(excited_Js, excited_Ks, delta_K):
    #     if del_K == -1:
    #         excited_E = excited_B * J * (J + 1) + (excited_C - excited_B) * (K ** 2) - (
    #             (-2 * excited_C * zeta)) * K + excited_C ** 2
    #     elif del_K == 1:
    #         excited_E = excited_B * J * (J + 1) + (excited_C - excited_B) * (K ** 2) + (
    #             (-2 * excited_C * zeta)) * K + excited_C ** 2

    #     excited_Es.append(excited_E)
    # start2_3 = timeit.default_timer()
    # linelist['excited_Es'] = excited_Es

    # wavenos = []
    # for i in range(len(linelist.index)):
    #     wavenumber = origin + excited_Es[i] - ground_Es[i]
    #     wavenos.append(wavenumber)

    # linelist['wavenos'] = wavenos



#%%

   # for J, K, delta_J, delta_K in zip(ground_Js, ground_Ks, delta_J, delta_K):
   #     if delta_J == -1 and delta_K == -1:
   #         HL_factor = ((J - 1 + K) * (J + K)) / (J * ((2 * J) + 1))
   #     elif delta_J == -1 and delta_K == 1:
   #         HL_factor = ((J - 1 - K) * (J - K)) / (J * ((2 * J) + 1))
   #     elif delta_J == 0 and delta_K == -1:
   #         HL_factor = (J + 1 - K) * (J + K) / (J * (J + 1))
   #     elif delta_J == 0 and delta_K == 1:
   #         HL_factor = (J + 1 + K) * (J - K) / (J * (J + 1))
   #     elif delta_J == 1 and delta_K == -1:
   #         HL_factor = (J + 2 - K) * (J + 1 - K) / ((J + 1) * ((2 * J) + 1))
   #     elif delta_J == 1 and delta_K == 1:
   #         HL_factor = (J + 2 + K) * (J + 1 + K) / ((J + 1) * ((2 * J) + 1))

   #     HL_factors.append(HL_factor)

#%%

    # for J, K, E in zip(ground_Js, ground_Ks, ground_Es):
    #     if K == 0:
    #         boltzmann_equation = (2 * ((2 * J) + 1)) * (np.exp((-h * c * E) / (k * T)))
    #     else:
    #         boltzmann_equation = (1 * ((2 * J) + 1)) * (np.exp((-h * c * E) / (k * T)))

    #     BD_factors.append(boltzmann_equation)

    # linelist['BD_factors'] = BD_factors
