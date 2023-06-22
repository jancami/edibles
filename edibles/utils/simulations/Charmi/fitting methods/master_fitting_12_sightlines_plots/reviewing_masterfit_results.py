#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 30 09:57:07 2023

@author: charmibhatt
"""

import numpy as np
import pandas as pd
import astropy.constants as const
import matplotlib.pyplot as plt
import timeit
# import scipy
# import scipy.stats as ss
from lmfit import Model
import csv
# import lmfit
import numba as nb
from pathlib import Path
from lmfit import Parameters

flux_list = np.array([])
wave_list = np.array([])
stddev_array = np.array([])
#spec_dir = Path("/Users/charmibhatt/Library/CloudStorage/OneDrive-TheUniversityofWesternOntario/UWO_onedrive/Local_GitHub/edibles/edibles/utils/simulations/Charmi/Heather's_data/")
spec_dir = Path("/Users/charmibhatt/Library/CloudStorage/OneDrive-TheUniversityofWesternOntario/UWO_onedrive/Research/Cami_2004_data/heliocentric/6614/")




combinations = pd.read_csv(r"/Users/charmibhatt/Library/CloudStorage/OneDrive-TheUniversityofWesternOntario/UWO_onedrive/Local_GitHub/edibles/edibles/utils/simulations/Charmi/Jmax=300.txt", delim_whitespace=(True))

startl = timeit.default_timer()


def get_rotational_spectrum(B, delta_B, zeta, T, sigma, origin):
    startg = timeit.default_timer()

    # print(B)
    # print(T)
    # print('-----------')

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

    '''Calculating Linelist'''

    ground_Es = []
    for J, K in zip(ground_Js, ground_Ks):
        ground_E = ground_B * J * (J + 1) + (ground_C - ground_B) * (K ** 2)
        ground_Es.append(ground_E)

    linelist['ground_Es'] = ground_Es

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

    endl = timeit.default_timer()
    #print('>>>> linelist calculation takes   ' + str(endl - startl) + '  sec')

    '''Smoothening the linelist'''

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
    simu_intenisty = 1 - 0.1 * (smooth_data[:, 1] / max(smooth_data[:, 1]))
    
    # with sns.color_palette("flare", n_colors=2):

    # plt.plot(simu_waveno, simu_intenisty, color='red')

    # for units in wavelngth
    # simu_wavelength = (1/simu_waveno)*1e8
    # model_data = np.array([simu_wavelength, simu_intenisty]).transpose()
    # model_data = model_data[::-1]

    model_data = np.array([simu_waveno, simu_intenisty]).transpose()
    
    
    endg = timeit.default_timer()
    plt.show()
    #print('>>>> full takes   ' + str(endg -startg) + '  sec') 
    
    
    return model_data

# Create an instance of the Variables class
variables = {}

# Open the lmfit report file and read its content
with open('fitting_12_sightlines.txt', 'r') as file:
    lines = file.readlines()
    parsing_variables = False
    for line in lines:
        if line.startswith("[[Variables]]"):
            parsing_variables = True
        elif parsing_variables and not line.startswith("#") and line.strip() != "":
            key_value = line.split(":")
            key = key_value[0].strip()
            if len(key_value) == 2:
                value = key_value[1].split("+/-")[0].strip()
                variables[key] = float(value)
            else:
                variables[key] = None
                
def curve_to_fit_wavenos(sightline): 
        
        file = 'hd{}_dib6614.txt'.format(sightline)
        Obs_data = pd.read_csv(spec_dir / file,
                               delim_whitespace=(True))
        Obs_data['Wavelength'] = (1 / Obs_data['Wavelength']) * 1e8
        Obs_data = Obs_data.iloc[::-1].reset_index(
            drop=True)  # making it ascending order as we transformed wavelength into wavenumbers

        # shifting to zero and scaling flux between 0.9 and 1
        min_index = np.argmin(Obs_data['Flux'])
        Obs_data['Wavelength'] = Obs_data['Wavelength'] - Obs_data['Wavelength'][min_index] 
        Obs_data['Flux'] = (Obs_data['Flux'] - min(Obs_data['Flux'])) / (1 - min(Obs_data['Flux'])) * 0.1 + 0.9

        data_to_plot = np.array([Obs_data['Wavelength'], Obs_data['Flux'] ]).transpose()
        
        # removing red wing
        # Obs_data_trp = Obs_data [(Obs_data['Wavelength'] >= -1) & (Obs_data['Wavelength']<= 1.2)]
        Obs_data_trp = Obs_data[(Obs_data['Flux'] <= 0.95)]  # trp = triple peak structure
        
        # making data evenly spaced
        x_equal_spacing = np.linspace(min(Obs_data_trp['Wavelength']), max(Obs_data_trp['Wavelength']), 25)
        y_obs_data = np.interp(x_equal_spacing, Obs_data_trp['Wavelength'], Obs_data_trp['Flux'])
        Obs_data_continuum = Obs_data [(Obs_data['Wavelength'] >= 2) & (Obs_data['Wavelength']<= 5)]
        std_dev = np.std(Obs_data_continuum['Flux'])
        
        return x_equal_spacing, y_obs_data, std_dev, data_to_plot


#sightlines = ['23180', '24398', '144470', '147165' , '147683', '149757', '166937', '170740', '184915', '185418', '185859', '203532']
sightlines = ['144217', '144470', '145502', '147165', '149757', '179406', '184915']

red_chis = []
for i, sightline in enumerate(sightlines, start=1):
    
    #print(sightline)
    print(i)
    B = variables.get('B', None)
    delta_B = variables.get('delta_B', None)
    zeta = variables.get('zeta', None)
    T = variables.get('T{}'.format(i), None)
    sigma = variables.get('sigma{}'.format(i), None)
    origin = variables.get('origin{}'.format(i), None)
    
    x_equal_spacing, y_obs_data, std_dev, data_to_plot = curve_to_fit_wavenos(sightline)
    
    model_data = get_rotational_spectrum(B, delta_B, zeta, T, sigma, origin)
    
    one_sl_y_model_data  = np.interp(x_equal_spacing, model_data[:, 0], model_data[:, 1])
    
    
    num = (one_sl_y_model_data - y_obs_data)**2
   
    chi_squared = np.sum((num)/(std_dev)**2)
    
    #print(chi_squared)
    
    reduced_chi_squared = chi_squared/(len( one_sl_y_model_data)*12 - 39)
    #reduced_chi_squared = chi_squared/(len( one_sl_y_model_data) - 6)

    print('chi_squared:  ' + str(chi_squared))
    
    red_chis.append(reduced_chi_squared)

    # indi_fits = pd.read_excel(r"/Users/charmibhatt/Library/CloudStorage/OneDrive-TheUniversityofWesternOntario/UWO_onedrive/Local_GitHub/edibles/edibles/utils/simulations/Charmi/fitting methods/individual_best_fits.xlsx")
    # indi_B = indi_fits['B'][i-1]
    # indi_delta_B = indi_fits['delta_B'][i-1]
    # indi_zeta = indi_fits['zeta'][i-1]
    # indi_sigma = indi_fits['sigma'][i-1]
    # indi_T = indi_fits['T'][i-1]
    # indi_origin = indi_fits['origin'][i-1]
    # #print(indi_fits['result_name'])
    # indi_model_data = get_rotational_spectrum(indi_B, indi_delta_B, indi_zeta, indi_T, indi_sigma, indi_origin)
    
    
    plt.figure(figsize = (15,8))
    plt.plot(data_to_plot[:,0], data_to_plot[:,1], color = 'black' , label = str('HD') + str(sightline ))
    # #plt.plot(x_equal_spacing, one_sl_y_model_data)
    plt.plot(model_data[:,0], model_data[:,1], color = 'red', label = str(r"$\chi^2$ =  ") + str('{:.3f}'.format(reduced_chi_squared)) + str('  (Altogether)') )
    
    
    # plt.plot(indi_model_data[:,0], indi_model_data[:,1], color = 'green', label = str(r"$\chi^2$ =  ") + str('{:.3f}'.format(indi_fits['redchi'][i-1])) + str( '  (Individual)'))
    # #plt.title( '  ground_B =  ' + str(B) + ' cm-1   Delta_B = ' + str(delta_B) +   '    zeta = ' +  str(zeta) + ' Temperature = ' + str(T) + '  K   $\sigma$ = ' + str(sigma) +    '    origin= ' +  str(origin)+ '\n' + '  ground_B =  ' + str(indi_B) + ' cm-1   Delta_B = ' + str(indi_delta_B) +   '    zeta = ' +  str(indi_zeta) + ' Temperature = ' + str(indi_T) + '  K   $\sigma$ = ' + str(indi_sigma) +    '    origin= ' +  str(indi_origin) ) 
    # #plt.text(0.25, 1,  ) 
    plt.title('Altogther (Cami 2004): ground_B = {:.5f} cm-1   Delta_B = {:.5f}    zeta = {:.5f} Temperature = {:.5f} K   $\sigma$ = {:.5f}    origin= {:.5f}\n\n'.format(B, delta_B, zeta, T, sigma, origin)) 
    #       'Individual: ground_B = {:.5f} cm-1   Delta_B = {:.5f}    zeta = {:.5f} Temperature = {:.5f} K   $\sigma$ = {:.5f}    origin= {:.5f}'.format(indi_B, indi_delta_B, indi_zeta, indi_T, indi_sigma, indi_origin))

    plt.legend(loc = 'lower left')
      
    print(sightline)

    # print('--------')

# print('Average reduced chi2:  ', np.sum(red_chis))



# indi_fits = pd.read_excel(r"/Users/charmibhatt/Library/CloudStorage/OneDrive-TheUniversityofWesternOntario/UWO_onedrive/Local_GitHub/edibles/edibles/utils/simulations/Charmi/fitting methods/individual_best_fits.xlsx")

# print(indi_fits['B'])
