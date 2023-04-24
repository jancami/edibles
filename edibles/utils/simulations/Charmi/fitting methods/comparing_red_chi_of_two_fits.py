#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 21 11:12:16 2023

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
import numpy as np
import pandas as pd
import astropy.constants as const
import matplotlib.pyplot as plt
import timeit
import scipy
import scipy.stats as ss
from lmfit import Model
import csv
import lmfit


Obs_data = pd.read_csv(r"/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/Heather's_data/6614_HD23180.txt", sep = ',')
Obs_data['Wavelength'] = (1/Obs_data['Wavelength'])*1e8
Obs_data = Obs_data.iloc[::-1].reset_index(drop=True) #making it ascending order as we transformed wavelength into wavenumbers



#shifting to zero and scaling flux between 0.9 and 1
min_index = np.argmin(Obs_data['Flux'])
Obs_data['Wavelength'] = Obs_data['Wavelength'] - Obs_data['Wavelength'][min_index]
Obs_data['Flux']=  (Obs_data['Flux'] - min(Obs_data['Flux'])) / (1 - min(Obs_data['Flux'])) * 0.1 + 0.9

#removing red wing
#Obs_data_trp = Obs_data [(Obs_data['Wavelength'] >= -1) & (Obs_data['Wavelength']<= 1.2)]
Obs_data_trp = Obs_data [(Obs_data['Flux'] <= 1)] #trp = triple peak structure


#making data evenly spaced
x_equal_spacing = np.linspace(min(Obs_data_trp['Wavelength']), max(Obs_data_trp['Wavelength']), len(Obs_data_trp['Wavelength']))
y_obs_data = np.interp(x_equal_spacing, Obs_data_trp['Wavelength'], Obs_data_trp['Flux'])

plt.plot(x_equal_spacing, y_obs_data)



def get_rotational_spectrum(xx, B, T, delta_B, zeta, sigma, origin):
    
    startg = timeit.default_timer()
    
    # print(B)
    # #print(T)
    # print('-----------')
   
    #rotational constants in cm-1
    ground_B = B
    ground_C = ground_B/2
    delta_C = delta_B
    excited_B = ground_B + ((delta_B/100)*ground_B)
    excited_C = ground_C + ((delta_C/100)*ground_C)
    
    combinations = pd.read_csv(r"/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/Jmax=300.txt", delim_whitespace=(True))
    
    ground_Js = combinations['ground_J']
    excited_Js = combinations['excited_J']
    ground_Ks = combinations['ground_K']
    excited_Ks = combinations['excited_K']
    
    linelist = combinations
   
    delta_J = linelist['excited_J'] - linelist['ground_J']
    delta_K = linelist['excited_K'] - linelist['ground_K']
    
    '''Calculating Linelist'''
    
    
    ground_Es = []
    for J,K in zip(ground_Js, ground_Ks):
                ground_E = ground_B*J*(J + 1) + (ground_C - ground_B)*(K**2)
                ground_Es.append(ground_E)
                
    linelist['ground_Es'] = ground_Es 
    
    excited_Es = []
    for J,K, del_K in zip(excited_Js, excited_Ks, delta_K):
        if del_K == -1:
                excited_E = excited_B*J*(J + 1) + (excited_C - excited_B)*(K**2) - ((-2*excited_C*zeta))*K + excited_C**2
        elif del_K == 1:
                excited_E = excited_B*J*(J + 1) + (excited_C - excited_B)*(K**2) + ((-2*excited_C*zeta))*K + excited_C**2

        excited_Es.append(excited_E)
                
    linelist['excited_Es'] = excited_Es
    
    wavenos = []
    for i in range(len(linelist.index)):
        wavenumber = origin + excited_Es[i] - ground_Es[i]
        wavenos.append(wavenumber)
        
    linelist['wavenos'] = wavenos

    
    HL_factors = []
    
    for J,K, delta_J, delta_K in zip(ground_Js, ground_Ks, delta_J, delta_K):
                if delta_J == -1 and delta_K == -1:
                    HL_factor = ((J - 1 + K)*(J + K))/(J*((2*J)+1))
                elif delta_J  == -1 and delta_K == 1:
                    HL_factor = ((J - 1 - K)*(J - K))/(J*((2*J)+1))
                elif delta_J == 0 and delta_K == -1:
                    HL_factor = (J+1-K)*(J+K)/ (J*(J + 1))
                elif delta_J == 0 and delta_K == 1:
                    HL_factor = (J+1+K)*(J-K) / (J*(J + 1))
                elif delta_J == 1 and delta_K == -1:
                    HL_factor = ( J + 2 - K)*( J + 1 - K)/ ((J + 1)*((2*J) +1))
                elif delta_J == 1 and delta_K == 1:
                    HL_factor = ( J + 2 + K)*( J + 1 + K)/ ((J + 1)*((2*J) + 1))
    
                HL_factors.append(HL_factor)
    
    linelist['HL_factors'] = HL_factors
                
    BD_factors = []
            
    h = const.h.cgs.value
    c = const.c.to('cm/s').value
    k = const.k_B.cgs.value

    
    for J,K,E in zip(ground_Js, ground_Ks, ground_Es):
        if K == 0:
            boltzmann_equation = (2*((2*J) + 1))*(np.exp((-h * c * E) / (k*T)))
        else:
            boltzmann_equation = (1*((2*J) + 1))*(np.exp((-h * c * E) / (k*T)))
            
        BD_factors.append(boltzmann_equation)
        
    linelist['BD_factors'] = BD_factors
    
    intensities = [] 
    for i in range(len(linelist.index)):
                strength = (HL_factors[i] * BD_factors[i])
                intensities.append(strength)
      
    linelist['intensities'] = intensities
    
    endl = timeit.default_timer()
    #print('>>>> linelist calculation takes   ' + str(endl-startl) + '  sec')
    
   
    '''Smoothening the linelist'''
    
    

   
    smooth_wavenos = np.linspace(np.min(linelist['wavenos']) - 1 ,np.max(linelist['wavenos']) + 1, 1000) # grid_size
    
    
    Wavenos_arr= np.array(linelist['wavenos'])
    Intenisty_arr = np.array(linelist['intensities'])
    

    import numba as nb

    @nb.njit(parallel=True)
    def calculate_smooth_intensities(wavenos, intensities, smooth_wavenos, sigma):
        smooth_intensities = np.zeros(smooth_wavenos.shape)
        
        for i in nb.prange(len(smooth_wavenos)):
            wavepoint = smooth_wavenos[i]
            w_int = np.exp(-(wavenos - wavepoint)**2 / (2*sigma**2)) * intensities
            smooth_intensities[i] = np.sum(w_int)
            
        return smooth_intensities
    
    # call the numba function with input data
    smooth_intensities = calculate_smooth_intensities(Wavenos_arr, Intenisty_arr, smooth_wavenos, sigma)

    
   
    
    smooth_data = np.array([smooth_wavenos, smooth_intensities]).transpose()    
    smooth_data = np.delete(smooth_data, np.where(smooth_data[:,1] <= 0.001*(max(smooth_data[:,1]))), axis = 0)
    
    simu_waveno = smooth_data[:, 0]
    simu_intenisty = 1-0.1*(smooth_data[:, 1]/max(smooth_data[:, 1]))
    
    
    #with sns.color_palette("flare", n_colors=2):
        
        
    #plt.plot(simu_waveno, simu_intenisty, color = 'red')
    
    
    #for units in wavelngth
    #simu_wavelength = (1/simu_waveno)*1e8
    #model_data = np.array([simu_wavelength, simu_intenisty]).transpose()
    #model_data = model_data[::-1]

    model_data = np.array([simu_waveno, simu_intenisty]).transpose()
    y_model_data = np.interp(x_equal_spacing, model_data[:,0], model_data[:,1])
    
    plt.figure(figsize=(15,8))
    plt.plot(x_equal_spacing, y_obs_data, color= 'green')
    plt.plot(x_equal_spacing, y_model_data)
    plt.xlabel('Wavelength')
    plt.ylabel('Normalized Intenisty')
    plt.title('Temperature = ' + str(T) + '  K  ground_B =  ' + str(ground_B) + ' cm-1  ground_C=  ' + str(ground_C) + ' cm-1  Delta_B = ' + str(delta_B) + '    $\sigma$ = ' + str(sigma) +    '    zeta = ' +  str(zeta)) 
    
    plt.show()
    
    plt.figure(figsize=(15,8))

    #with sns.color_palette("flare", n_colors=2):

    plt.stem(x_equal_spacing, y_obs_data - y_model_data)
    plt.axhline(y=0)
    endg = timeit.default_timer()
     
    #print('>>>> full takes   ' + str(endg -startg) + '  sec') 
    
    return  y_model_data

xx= x_equal_spacing
B = 0.0023
T = 100
delta_B = -0.0353
zeta = -0.4197
sigma = 0.2247
origin =  0.0041




num = (get_rotational_spectrum(xx, B, T, delta_B, zeta, sigma, origin) - y_obs_data)**2
chi_squared = np.sum((num)/(0.004)**2)
reduced_chi_squared = chi_squared/(len(y_obs_data) - 6)

print('-------')
print(reduced_chi_squared)


xx= x_equal_spacing
B = 0.0022
T = 96.09
delta_B = -0.03
zeta = -0.415
sigma = 0.22
origin =  0.008



num = (get_rotational_spectrum(xx, B, T, delta_B, zeta, sigma, origin) - y_obs_data)**2
chi_squared = np.sum((num)/(0.004)**2)
reduced_chi_squared = chi_squared/(len(y_obs_data) - 6)

print('-------')
print(reduced_chi_squared)











