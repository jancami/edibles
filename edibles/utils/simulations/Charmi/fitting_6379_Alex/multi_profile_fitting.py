# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 11:44:35 2023

@author: alexr
"""
import functions as fn
import numpy as np
from lmfit import Model
import csv
import matplotlib.pyplot as plt
from lmfit import Parameters
import timeit

sightlines = ['24398','144470','147165','147683','149757','166937',
              '170740','184915','185418','203532','185859']

method = 'leastsq'
flux_list = np.array([])
wave_list = np.array([])
stddev_array = np.array([])
Jmax = 300
combinations = fn.allowed_perperndicular_transitions(Jmax)


#%%

for sightline in sightlines:
    
    Obs_data, x_equal_spacing, y_data_to_fit, std_dev = fn.obs_curve_to_fit(sightline)

    flux_list = np.concatenate((flux_list, y_data_to_fit))
    wave_list = np.concatenate((wave_list, x_equal_spacing))
    
    one_sl_stddev = [std_dev] * len(x_equal_spacing)
    stddev_array = np.concatenate((stddev_array, one_sl_stddev))
  
plt.plot(wave_list, flux_list)

#%% Generate multiple spectra

def get_multi_spectra( **params_list):
    """
    Calculating a model for each sight line using 'get_rotational_spectrum'.
   
    Always using the same molecular parameters, but different T.
    Args:
        xx:
        B:
        T1:
        T2:
        delta_B:
        zeta:
        sigma:
        origin:

    Returns:
    np.array
        Model flux array. Fluxes for both sight lines are appended to one 1D array.
    """
    
   
    print('---------')
    
    B = params_list['B']
    delta_B = params_list['delta_B']
    zeta = params_list['zeta']
    
   
    # first_T_index = 3
    # last_T_index = first_T_index + len(sightlines) 
 
    # first_sigma_index = last_T_index  
    # last_sigma_index = first_sigma_index + len(sightlines) 
 
    # first_origin_index = last_sigma_index  
    # last_origin_index = first_origin_index +len(sightlines) 
 
    # T_values = params_list[first_T_index:last_T_index]
    # sigma_values = params_list[first_sigma_index:last_sigma_index]
    # origin_values = params_list[first_origin_index:last_origin_index]
    
    T_values = [params_list[f'T{i+1}'] for i in range(len(sightlines))]
    sigma_values = [params_list[f'sigma{i+1}'] for i in range(len(sightlines))]
    origin_values = [params_list[f'origin{i+1}'] for i in range(len(sightlines))]

    all_y_model_data = np.array([])
    
    start = timeit.default_timer()
    
    for T, sigma, origin, sightline in zip(T_values, sigma_values, origin_values, sightlines):
        
        x_model_data, y_model_data = fn.get_rotational_spectrum(B, delta_B, zeta, T, sigma, origin, combinations, transition = 'perpendicular', bell = False)
        
        Obs_data, x_equal_spacing, y_data_fit, std_dev = fn.obs_curve_to_fit(sightline)
        
        # plt.plot(x_model_data, y_model_data, label = 'Model')
        # Obs_data = Obs_data[Obs_data['Flux']<=0.95]
        # plt.plot(Obs_data['Wavelength'], Obs_data['Flux'], label = 'Raw obs, HD{}'.format(sightline))
        # plt.plot(x_equal_spacing, y_data_fit, label = 'Interpolated obs')
        # plt.legend()
        # plt.figsize = (1,1.5)
        # plt.show()
        
        one_sl_y_model_data  = np.interp(x_equal_spacing, x_model_data, y_model_data)
        
        all_y_model_data = np.concatenate((all_y_model_data, one_sl_y_model_data))
        
    end = timeit.default_timer()
    print('Time for profile of this iteration of all sightlines ' + str(end - start))
    print('==========')
    print('==========')
    print('==========')
    # plt.plot(common_grid_for_all, all_y_model_data)
    # plt.show()
    return all_y_model_data


#%% Fit multiple spectra to a model

def fit_model_multi(B, delta_B, zeta, T, sigma, origin, combinations, transition, Jmax):
    mod = Model(get_multi_spectra, 
                independent_vars=['combinations', 'transition', 'Jmax']) 
    
    
    
    params_list = [B, delta_B, zeta]
    
    T_list = [T] * len(sightlines)
    sigma_list = [sigma] * len(sightlines)
    origin_list = [origin] * len(sightlines)

    params_list.extend(T_list)
    params_list.extend(sigma_list)
    params_list.extend(origin_list)
    
    
    print(params_list)
    
    
    
    first_T_index = 3
    last_T_index = first_T_index + len(sightlines) 
 
    first_sigma_index = last_T_index  
    last_sigma_index = first_sigma_index + len(sightlines) 
 
    first_origin_index = last_sigma_index  
    last_origin_index = first_origin_index +len(sightlines) 
 
    params = Parameters()
    params.add('B', value = B, min = 0.0005, max = 0.05)
    params.add('delta_B', value = delta_B, min = -1, max =0)
    params.add('zeta', value = zeta, min = -1, max = 1)
    
    for i, param_value in enumerate(params_list[first_T_index:last_T_index]):
        params.add(f'T{i+1}', value=param_value, min = 2.7, max = 500)
        
    for i, param_value in enumerate(params_list[first_sigma_index:last_sigma_index]):
        params.add(f'sigma{i+1}', value=param_value, min = 0.05, max = 0.3)
        
    for i, param_value in enumerate(params_list[first_origin_index:last_origin_index]):
        params.add(f'origin{i+1}', value=param_value, min = -1, max = 1)
        
   
    result = mod.fit(flux_list,
                     params, 
                     xx=wave_list, 
                     weights = 1/stddev_array, 
                     method = method,  
                     combinations = combinations, 
                     transition = transition, 
                     Jmax=Jmax) #, fit_kws={'ftol': 1e-2, 'xtol': 1e-2} )
    print(result.fit_report())
    
    # def plot_best_fit(result, x_equal_spacing, y_obs_data):
    #     plt.figure()
    #     plt.scatter(x_equal_spacing, y_obs_data, label='Observations')
    #     plt.plot(x_equal_spacing, result.best_fit, 'r-', label='Best Fit')
    #     plt.xlabel('x')
    #     plt.ylabel('y')
    #     plt.legend()
    #     plt.show()
            
    #plot_best_fit(result, x_equal_spacing, y_obs_data)
    
    return result

#%% 

def write_results_to_csv(results_list, filename):
    with open(filename, 'w', newline='') as csvfile:
        fieldnames = ['result_name', 
                      'B_init', 'delta_B_init', 'zeta_init' , 'T1_init', 'sigma1_init' ,'origin1_init' , 
                      'B',   'delta_B', 'zeta', 'T1','sigma1', 'origin1', 'chi2', 'redchi', 'func_evals', 
                      'B_unc', 'delta_B_unc', 'zeta_unc', 'T1_unc',  'sigma1_unc', 'origin1_unc']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for i, result in enumerate(results_list):
            result_name = f'result{i+1}'
            params = result.params
            row = {
                'result_name': result_name,
                'B_init': params['B'].init_value,
                'delta_B_init': params['delta_B'].init_value,
                'zeta_init': params['zeta'].init_value,
                'T1_init': params['T1'].init_value,
                'sigma1_init': params['sigma1'].init_value,
                'origin1_init': params['origin1'].init_value,
                'B': params['B'].value,
                'delta_B': params['delta_B'].value,
                'zeta': params['zeta'].value,
                'T1': params['T1'].value,
                'sigma1': params['sigma1'].value,
                'origin1': params['origin1'].value,
                'chi2': result.chisqr,
                'redchi': result.redchi,
                'func_evals': result.nfev,
                'B_unc': params['B'].stderr,
                'delta_B_unc': params['delta_B'].stderr,
                'zeta_unc': params['zeta'].stderr,
                'T1_unc': params['T1'].stderr,
                'sigma1_unc': params['sigma1'].stderr,
                'origin1_unc': params['origin1'].stderr,

            }
            writer.writerow(row)

#%% 
B = 0.0014747
delta_B = 0.1460773
zeta = 0.188564
T = 299.91933
sigma = 0.1934417
origin = -0.258784

#%%

result = fit_model_multi(B, delta_B, zeta, T, sigma, origin, combinations, transition = 'parallel', Jmax = Jmax)
