#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 23:54:55 2023

@author: charmibhatt
"""
import numpy as np
import pandas as pd
import astropy.constants as const
import matplotlib.pyplot as plt
import timeit
import scipy
import scipy.stats as ss
from scipy.signal import argrelextrema



#object_names = ['23180', '24398', '144470', '147165' , '147683', '149757', '166937', '170740', '184915', '185418', '185859', '203532']
object_names = ['185418']

for on in object_names:
    print('-----------')
    print(on)
    #Fit_data = pd.read_csv(r"/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/fitting methods/results_{}_triple_peak_below_96.csv".format(on))
    Fit_data = pd.read_csv(r"/Users/charmibhatt/Library/CloudStorage/OneDrive-TheUniversityofWesternOntario/UWO_onedrive/Local_GitHub/edibles/edibles/utils/simulations/Charmi/fitting methods/results_5780_HD185418.csv")
    print(Fit_data)
    min_value = Fit_data['redchi'].min()
    Best_fit = Fit_data[Fit_data['redchi'] == min_value]
    
    #Obs_data = pd.read_csv(r"/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/Heather's_data/6614_HD{}.txt".format(on), sep = ',')
    #Obs_data = pd.read_csv(r"/Users/charmibhatt/Documents/GitHub/DIBs/5780_fitting/DIB5780_HD185418.txt", sep = ' ')
    Obs_data = pd.read_csv(r"/Users/charmibhatt/Library/CloudStorage/OneDrive-TheUniversityofWesternOntario/UWO_onedrive/Local_GitHub/edibles/edibles/utils/simulations/Charmi/fitting methods/5780_fitting/DIB5780_HD185418.txt", delim_whitespace=(True))

    #Obs_data_trp = Obs_data [(Obs_data['Wavelength'] >= 6611) & (Obs_data['Wavelength']<= 6616)]

    #Obs_data['Wavelength'] = (1/Obs_data['Wavelength'])*1e8
    # Obs_data = Obs_data.iloc[::-1].reset_index(drop=True) #making it ascending order as we transformed wavelength into wavenumbers



    #shifting to zero and scaling flux between 0.9 and 1
    min_index = np.argmin(Obs_data['Flux'])
    Obs_data['Wavelength'] = Obs_data['Wavelength'] - Obs_data['Wavelength'][min_index]
    Obs_data['Flux']=  (Obs_data['Flux'] - min(Obs_data['Flux'])) / (1 - min(Obs_data['Flux'])) * 0.1 + 0.9

    #removing red wing
    #Obs_data_trp = Obs_data [(Obs_data['Wavelength'] >= -1) & (Obs_data['Wavelength']<= 1.2)]
    Obs_data_trp = Obs_data [(Obs_data['Flux'] <= 0.95)] #trp = triple peak structure


    #making data evenly spaced
    x_equal_spacing = np.linspace(min(Obs_data['Wavelength']), max(Obs_data['Wavelength']), len(Obs_data['Wavelength']))
    y_obs_data = np.interp(x_equal_spacing, Obs_data['Wavelength'], Obs_data['Flux'])

    #Obs_data_trp = Obs_data[(Obs_data['Flux'] <= 0.96)]
    Red_wing_cutoff = min(Obs_data_trp['Wavelength'])
    Blue_wing_cutoff = max(Obs_data_trp['Wavelength'])
    # print(Red_wing_cutoff)
    # print(Blue_wing_cutoff)
    
    

    
    
def get_rotational_spectrum(xx, B, T, delta_B, zeta, sigma, origin):
    
    
    startg = timeit.default_timer()
    
    print(B)
    print(T)
    
    #rotational constants in cm-1
    ground_B = B
    ground_C = ground_B/2
    delta_C = delta_B
    excited_B = ground_B + ((delta_B/100)*ground_B)
    excited_C = ground_C + ((delta_C/100)*ground_C)
    # min_index = np.argmin(Obs_data['Flux'])
    # origin = (1/Obs_data_trp['Wavelength'][min_index])*1e8 + origin
    
    
    combinations = pd.read_csv(r"/Users/charmibhatt/Library/CloudStorage/OneDrive-TheUniversityofWesternOntario/UWO_onedrive/Local_GitHub/edibles/edibles/utils/simulations/Charmi/Jmax=300.txt", delim_whitespace=(True))
    
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
    #smooth_data = np.delete(smooth_data, np.where(smooth_data[:,1] <= 0.001*(max(smooth_data[:,1]))), axis = 0)
    
    simu_waveno = smooth_data[:, 0]
    simu_intenisty = 1-0.1*(smooth_data[:, 1]/max(smooth_data[:, 1]))
    
    
    #with sns.color_palette("flare", n_colors=2):
    
    
    
    #for units in wavelngth
    # simu_wavelength = (1/simu_waveno)*1e8
    # model_data = np.array([simu_wavelength, simu_intenisty]).transpose()
    #model_data = model_data[::-1]
    #print(model_data)

    # model_data = np.array([simu_wavelength, simu_intenisty]).transpose()
    # x_equal_spacing = np.linspace(min(Obs_data['Wavelength']), max(Obs_data['Wavelength']), len(Obs_data['Wavelength']))
    # y_model_data = np.interp(x_equal_spacing, model_data[:,0], model_data[:,1])
    
    
    
    pd.options.display.float_format = '{:.4f}'.format


    # plt.subplot(2,1,1)
    # plt.plot(xx, y_obs_data, color= 'green')

   # plt.subplot(2,1,1)
    plt.figure(figsize = (17,8))
    #Obs_data_trp['Flux']=  (Obs_data_trp['Flux'] - min(Obs_data_trp['Flux'])) / (1 - min(Obs_data_trp['Flux'])) * 0.1 + 0.9

    #plt.plot(Obs_data['Wavelength'], Obs_data['Flux'], linewidth = 3, label = 'Data')
    #plt.show()

    #plt.plot(model_data[:,0], model_data[:,1], linewidth = 2.8, label = 'Model')
    plt.plot(simu_waveno, simu_intenisty, color = 'red')
    #plt.xlim(-7.5, 7.5)
    # plt.axvline(x = Red_wing_cutoff, color = 'black', alpha = 0.5)
    # plt.axvline(x = Blue_wing_cutoff, color = 'black', alpha = 0.5)
    plt.xlabel('Wavelength ($\AA$)' , size  = 15)
    plt.ylabel('Normalized Intenisty', size = 15, labelpad = 10)
    plt.xticks(size = 15)
    plt.yticks(size = 15)
    #plt.title('Temperature = {:.4f} K, B = {:.4f} cm$^{-1}$, $\Delta$B = {:.4f}%, $\zeta$ = {:.4f}, $\sigma$ = {:.4f} cm$^{-1}$, reduced $\chi^2$ = {:.4f}'.format(T, ground_B, -1, delta_B, zeta, sigma, -1, zeta, redchi), size='x-large')
    #plt.title('Temperature = {:.4f}'.format(T) + 'K,     B = {:.4f}'.format(B) + ' cm$^{-1}$     ' + r'$\Delta B =$ {:.4f}'.format(delta_B) + '%,     'r'$\zeta^\prime  = $ {:.4f}'.format(zeta) + ',     'r'$\sigma = $ {:.4f}'.format(sigma) + 'cm$^{-1}$' + '    reduced $\chi^2$ = {:.4f}'.format(redchi), size =15)

    #plt.xlim(-5, 3)
    plt.legend(loc = 'lower left', fontsize = 12)
    
    # plt.subplot(2,1,2)    
    # plt.stem(xx, y_obs_data - y_model_data)
    # plt.axhline(y = 0)
    # plt.axvline(x = Red_wing_cutoff, color = 'black')
    # plt.axvline(x = Blue_wing_cutoff, color = 'black')
    # plt.xlim(-7.5, 7.5)
    # plt.show()
    

    endg = timeit.default_timer()
     
    print('>>>> full takes   ' + str(endg -startg) + '  sec') 
    
    #return  y_model_data

B = Best_fit['B'].item()
T = Best_fit['T'].item()
delta_B = Best_fit['delta_B'].item()
zeta = Best_fit['zeta'].item()
sigma = Best_fit['sigma'].item()
origin = Best_fit['origin'].item()
redchi = Best_fit['redchi'].item()
xx = x_equal_spacing

get_rotational_spectrum(xx, B, T, delta_B, zeta, sigma, origin)

#get_rotational_spectrum(np.linspace(-5,5), B = 0.0099, T = 102.8485, delta_B = -0.003, zeta = -0.1705, sigma = 0.05, origin = 0.056)

    # num = (get_rotational_spectrum(xx, B, T, delta_B, zeta, sigma, origin) - y_obs_data)**2
    # chi_squared = np.sum((num)/(0.003)**2)
    # reduced_chi_squared = chi_squared/(len(y_obs_data) - 6)
    
    # print('-------')
    # print(reduced_chi_squared)
