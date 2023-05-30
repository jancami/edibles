#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  3 12:20:38 2023

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
import scipy as sp
import scipy.stats as ss
from lmfit import Model
import csv
import lmfit
from lmfit import minimize, Parameters, report_fit

#object_error = [0.00217, 0.00223, 0.00122, 0.00429]
object_error = [0.00223, 0.00122]

#creating master dataset

#object_names = ['23180', '166937', '185418', '203532' ]

object_names = ['166937', '185418' ]

#object_names = ['23180', '24398', '144470', '147165' , '147683', '149757', '166937', '170740', '184915', '185418', '185859', '203532']

mdata= np.array([])
for on in object_names:  
     print(on)
     Obs_data = pd.read_csv(r"/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/Heather's_data/6614_HD{}.txt".format(on), sep = ',')
     
     Obs_data['Wavelength'] = (1/Obs_data['Wavelength'])*1e8
     Obs_data = Obs_data.iloc[::-1].reset_index(drop=True) #making it ascending order as we transformed wavelength into wavenumbers
   
   
   
    #shifting to zero and scaling flux between 0.9 and 1
     min_index = np.argmin(Obs_data['Flux'])
     Obs_data['Wavelength'] = Obs_data['Wavelength'] - Obs_data['Wavelength'][min_index]
     Obs_data['Flux']=  (Obs_data['Flux'] - min(Obs_data['Flux'])) / (1 - min(Obs_data['Flux'])) * 0.1 + 0.9
   
     # print(min(Obs_data['Wavelength']))
     # print(max(Obs_data['Wavelength']))
     # print(len(Obs_data['Wavelength']))
     # print('------------------')

    #removing red wing
    #Obs_data_trp = Obs_data [(Obs_data['Wavelength'] >= -1) & (Obs_data['Wavelength']<= 2)]
     Obs_data_trp = Obs_data [(Obs_data['Flux'] <= 0.95)] #trp = triple peak structure
   
   
    #making data evenly spaced
     x_equal_spacing = np.linspace(min(Obs_data_trp['Wavelength']), max(Obs_data_trp['Wavelength']), 25) #len(Obs_data_trp['Wavelength']))
     y_obs_data = list(np.interp(x_equal_spacing, Obs_data_trp['Wavelength'], Obs_data_trp['Flux']))
   
     #print(type(y_obs_data))
     y_obs_data = np.array(y_obs_data)
     #print(y_obs_data.shape)
     #y_obs_data = Obs_data['Wavelength']
     #print(type(y_obs_data))
     mdata = np.append(mdata, y_obs_data)

mdata = mdata.reshape(len(object_names), -1)
#mdata = np.array(mdata)     
    
print(mdata.shape)



combinations = pd.read_csv(r"/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/Jmax=300.txt", delim_whitespace=(True))


startl = timeit.default_timer()
def get_rotational_spectrum(xx, B, T, delta_B, zeta, sigma, origin):
    
    startg = timeit.default_timer()
    print('-----------')
   
    print(B)
    #print(T)
    print('-----------')
   
    #rotational constants in cm-1
    ground_B = B
    ground_C = ground_B/2
    delta_C = delta_B
    excited_B = ground_B + ((delta_B/100)*ground_B)
    excited_C = ground_C + ((delta_C/100)*ground_C)
    
    
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
    print('>>>> linelist calculation takes   ' + str(endl-startl) + '  sec')
    
   
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
        
        
    plt.plot(simu_waveno, simu_intenisty, color = 'red')
    
    
    #for units in wavelngth
    #simu_wavelength = (1/simu_waveno)*1e8
    #model_data = np.array([simu_wavelength, simu_intenisty]).transpose()
    #model_data = model_data[::-1]

    model_data = np.array([simu_waveno, simu_intenisty]).transpose()
    y_model_data = np.interp(x_equal_spacing, model_data[:,0], model_data[:,1])
    
    # plt.figure(figsize=(15,8))
    # plt.plot(x_equal_spacing, y_obs_data, color= 'green')

    # plt.plot(x_equal_spacing, y_model_data)
    # plt.xlabel('Wavelength')
    # plt.ylabel('Normalized Intenisty')
    # plt.title('Temperature = ' + str(T) + '  K  ground_B =  ' + str(ground_B) + ' cm-1  ground_C=  ' + str(ground_C) + ' cm-1  Delta_B = ' + str(delta_B) + '    $\sigma$ = ' + str(sigma) +    '    zeta = ' +  str(zeta)) 
    
    # plt.show()
    
    # plt.figure(figsize=(15,8))

    #with sns.color_palette("flare", n_colors=2):
    plt.show()
    plt.stem(x_equal_spacing, y_obs_data - y_model_data)
    # # plt.axhline(y=0)
    endg = timeit.default_timer()
     
    print('>>>> full takes   ' + str(endg -startg) + '  sec') 
    
    # num = np.sum(y_model_data - y_obs_data)
    # chi_squared = np.sum((num)/(0.004)**2)
    # reduced_chi_squared = chi_squared/(len(y_model_data) - 6)
    # print(reduced_chi_squared)
    
    return  y_model_data

def get_rotational_spectrum_select_params(params, i, x_equal_spacing):
    x_equal_spacing = x_equal_spacing
    B = params[ 'B_%i' % (i+1)].value
    T = params[ 'T_%i' % (i+1)].value
    delta_B = params[ 'delta_B_%i' % (i+1)].value
    zeta = params[ 'zeta_%i' % (i+1)].value
    sigma = params[ 'sigma_%i' % (i+1)].value
    origin = params[ 'origin_%i' % (i+1)].value
    return get_rotational_spectrum(x_equal_spacing, B, T, delta_B, zeta, sigma, origin)

def objective(params, x_equal_spacing, mdata):
    """ calculate total residual for fits to several data sets held
    in a 2-D array, and modeled by Gaussian functions"""
    ndata, nx = mdata.shape
    # print(ndata)
    # print(nx)
    resid = 0.0*mdata[:]
    
    #print(resid)
    # make residual per data set
    weighted_set = []
    for i, oe in zip(range(ndata), object_error):
    #for i in range(ndata):
        resid[i, :] = mdata[i, :] - get_rotational_spectrum_select_params(params, i, x_equal_spacing)
        weighted = list(((resid[i,:] ** 2)) / (oe ** 2))
        weighted_set.append(weighted)
    weighted_set= np.array(weighted_set)
    #print(np.sum(weighted_set))
    weighted_set_flatten = weighted_set.flatten()
    print((weighted_set_flatten.shape))

    return weighted_set_flatten
    # now flatten this to a 1D array, as minimize() needs
    #return resid.flatten()





def fit_model(B, T, delta_B, zeta, sigma, origin):

    fit_params = Parameters()
    for iy, y in enumerate(mdata):
        
        fit_params.add( 'B_%i' % (iy+1), value= B, min=0.0005,  max=0.01)
        fit_params.add( 'T_%i' % (iy+1), value= T, min=-2.7,  max=200)
        fit_params.add( 'delta_B_%i' % (iy+1), value= delta_B, min=-1, max=0)
        fit_params.add( 'zeta_%i' % (iy+1), value= zeta, min=-1,  max=1)
        fit_params.add( 'sigma_%i' % (iy+1), value= sigma, min=-0.05,  max=0.3)
        fit_params.add( 'origin_%i' % (iy+1), value= origin, min=-1, max=1)
        
    #print(fit_params)
    
    #but now constrain all values of sigma to have the same value
    # by assigning sig_2, sig_3, .. sig_5 to be equal to sig_1
    for iy in range(2,len(object_names)+1):
        fit_params['B_%i' % iy].expr='B_1'
        fit_params['delta_B_%i' %iy].expr='delta_B_1'
        fit_params['zeta_%i' % iy].expr='zeta_1'
            
    #objective(fit_params, x_equal_spacing, mdata)
    
    
    result = minimize(objective, fit_params, args=(x_equal_spacing, mdata))
    report_fit(result)
    return result


result1= fit_model(B = 0.01, T = 2.7, delta_B = -0.1, zeta = -0.1, sigma = 0.02, origin =  0)
result2= fit_model(B = 0.005, T = 32.7, delta_B = -0.4, zeta = -0.9, sigma = 0.1, origin =  0.039)
result3= fit_model(B = 0.002, T = 53.26, delta_B = -0.132, zeta = -0.312, sigma = 0.166, origin =  0.04)
result4= fit_model(B = 0.0003, T = 32.5, delta_B = -0.45, zeta = -0.01, sigma = 0.17, origin =  0.012)
result5= fit_model(B = 0.0075, T = 92.5, delta_B = -0.23, zeta = -0.23, sigma = 0.23, origin =  0.02)

results_list = [result1, result2, result3, result4, result5]

def write_results_to_csv(results_list, filename):
    with open(filename, 'w', newline='') as csvfile:
        #fieldnames = ['result_name', 'B', 'delta_B', 'zeta', 'T_1',  'T_2', 'T_3','T_4', 'sigma_1', 'sigma_2', 'sigma_3', 'sigma_4', 'origin_1', 'origin_2', 'origin_3', 'origin_4', 'chi2', 'redchi', 'func_evals', 'B_unc', 'delta_B_unc', 'zeta_unc', 'T_1_unc', 'sigma_1_unc', 'origin_1_unc', 'T_2_unc', 'sigma_2_unc', 'origin_2_unc', 'B_init', 'T_init', 'delta_B_init', 'zeta_init' , 'sigma_init', 'origin_init' ]
        fieldnames = ['result_name', 'B', 'delta_B', 'zeta', 'T_1',  'T_2',  'sigma_1', 'sigma_2',  'origin_1', 'origin_2', 'chi2', 'redchi', 'func_evals', 'B_unc', 'delta_B_unc', 'zeta_unc', 'T_1_unc', 'sigma_1_unc', 'origin_1_unc', 'T_2_unc', 'sigma_2_unc', 'origin_2_unc', 'B_init', 'T_init', 'delta_B_init', 'zeta_init' , 'sigma_init', 'origin_init' ]

        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for i, result in enumerate(results_list):
            result_name = f'result{i+1}'
            params = result.params
            row = {
                'result_name': result_name,
                'B_init': params['B_1'].init_value,
                'T_init': params['T_1'].init_value,
                'delta_B_init': params['delta_B_1'].init_value,
                'zeta_init': params['zeta_1'].init_value,
                'sigma_init': params['sigma_1'].init_value,
                'origin_init': params['origin_1'].init_value,
                'B': params['B_1'].value,
                'T_1': params['T_1'].value,
                'T_2': params['T_2'].value,
                # 'T_3': params['T_3'].value,
                # 'T_4': params['T_4'].value,
                'delta_B': params['delta_B_1'].value,
                'zeta': params['zeta_1'].value,
                'sigma_1': params['sigma_1'].value,       
                'sigma_2': params['sigma_2'].value,
                # 'sigma_3': params['sigma_3'].value,
                # 'sigma_4': params['sigma_4'].value,
                'origin_1': params['origin_1'].value,
                'origin_2': params['origin_2'].value,
                # 'origin_3': params['origin_3'].value,
                # 'origin_4': params['origin_4'].value,
                'chi2': result.chisqr,
                'redchi': result.redchi,
                'func_evals': result.nfev,
                'B_unc': params['B_1'].stderr,
                'T_1_unc': params['T_1'].stderr,
                'T_2_unc': params['T_2'].stderr,
                # 'T_3_unc': params['T_3'].stderr,
                # 'T_4_unc': params['T_4'].stderr,
                'delta_B_unc': params['delta_B_1'].stderr,
                'zeta_unc': params['zeta_1'].stderr,
                'sigma_1_unc': params['sigma_1'].stderr,
                'sigma_2_unc': params['sigma_2'].stderr,
                # 'sigma_3_unc': params['sigma_3'].stderr,
                # 'sigma_4_unc': params['sigma_4'].stderr,
                'origin_1_unc': params['origin_1'].stderr,
                'origin_2_unc': params['origin_2'].stderr,
                # 'origin_3_unc': params['origin_3'].stderr,
                # 'origin_4_unc': params['origin_4'].stderr

            }
            writer.writerow(row)
            
            
            
write_results_to_csv(results_list, 'results_166937_185418_no_sqrt.csv')     
            


'''Checking resuts'''

# Obs_data = pd.read_csv(r"/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/Heather's_data/6614_HD23180.txt", sep = ',')
# Obs_data['Wavelength'] = (1/Obs_data['Wavelength'])*1e8
# Obs_data = Obs_data.iloc[::-1].reset_index(drop=True) #making it ascending order as we transformed wavelength into wavenumbers

# #shifting to zero and scaling flux between 0.9 and 1
# min_index = np.argmin(Obs_data['Flux'])
# Obs_data['Wavelength'] = Obs_data['Wavelength'] - Obs_data['Wavelength'][min_index]
# Obs_data['Flux']=  (Obs_data['Flux'] - min(Obs_data['Flux'])) / (1 - min(Obs_data['Flux'])) * 0.1 + 0.9

# # plt.figure(figsize = (15,8))
# plt.plot(Obs_data['Wavelength'], Obs_data['Flux'])
# # get_rotational_spectrum(x_equal_spacing, B = 0.00256, T = 79.69, delta_B = -0.0694, zeta = -0.3422, sigma= 0.17, origin = 0.099)
# # plt.show()

# Obs_data_trp = Obs_data [(Obs_data['Wavelength'] >= 2.5) & (Obs_data['Wavelength']<= 5)]
# print(np.std(Obs_data_trp['Flux']))
    

# Obs_data = pd.read_csv(r"/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/Heather's_data/6614_HD185418.txt", sep = ',')
# Obs_data['Wavelength'] = (1/Obs_data['Wavelength'])*1e8
# Obs_data = Obs_data.iloc[::-1].reset_index(drop=True) #making it ascending order as we transformed wavelength into wavenumbers

# #shifting to zero and scaling flux between 0.9 and 1
# min_index = np.argmin(Obs_data['Flux'])
# Obs_data['Wavelength'] = Obs_data['Wavelength'] - Obs_data['Wavelength'][min_index]
# Obs_data['Flux']=  (Obs_data['Flux'] - min(Obs_data['Flux'])) / (1 - min(Obs_data['Flux'])) * 0.1 + 0.9

# plt.figure(figsize = (15,8))
# plt.plot(Obs_data['Wavelength'], Obs_data['Flux'])
# get_rotational_spectrum(x_equal_spacing, B = 0.00256, T = 81.21, delta_B = -0.0694, zeta = -0.3422, sigma= 0.230, origin = 0.0725)
# # plt.xlim(2.5,5)
# # plt.ylim(0.99,1.01)
# plt.show()

