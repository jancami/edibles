#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 15 18:35:34 2023

@author: charmibhatt
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 12:28:52 2023

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

combinations = pd.read_csv(r"/Users/charmibhatt/Library/CloudStorage/OneDrive-TheUniversityofWesternOntario/UWO_onedrive/Local_GitHub/edibles/edibles/utils/simulations/Charmi/Jmax=300.txt", delim_whitespace=(True))
import numba as nb


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
    print('>>>> linelist calculation takes   ' + str(endl - startl) + '  sec')

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
    #smooth_data = np.delete(smooth_data, np.where(smooth_data[:, 1] <= 0.001 * (max(smooth_data[:, 1]))), axis=0)

    simu_waveno = smooth_data[:, 0]
    simu_intenisty = 1 - 0.1 * (smooth_data[:, 1] / max(smooth_data[:, 1]))
    
    # with sns.color_palette("flare", n_colors=2):

    # plt.plot(simu_waveno, simu_intenisty, color='red')

    # for units in wavelngth
    # simu_wavelength = (1/simu_waveno)*1e8
    # model_data = np.array([simu_wavelength, simu_intenisty]).transpose()
    # model_data = model_data[::-1]

    model_data = np.array([simu_waveno, simu_intenisty]).transpose()
    print(model_data)
    
    plt.plot(model_data[:,0], model_data[:,1])
    # plt.figure(figsize=(15,8))
    # plt.plot(x_equal_spacing, y_obs_data, color= 'green')

    # plt.plot(x_equal_spacing, y_model_data)
    # plt.xlabel('Wavelength')
    # plt.ylabel('Normalized Intenisty')
    # plt.title('Temperature = ' + str(T) + '  K  ground_B =  ' + str(ground_B) + ' cm-1  ground_C=  ' + str(ground_C) + ' cm-1  Delta_B = ' + str(delta_B) + '    $\sigma$ = ' + str(sigma) +    '    zeta = ' +  str(zeta)) 

    # plt.show()

    # plt.figure(figsize=(15,8))

    # with sns.color_palette("flare", n_colors=2):

    # plt.stem(x_equal_spacing, y_obs_data - y_model_data)
    # plt.axhline(y=0)
    endg = timeit.default_timer()

    print('>>>> full takes   ' + str(endg - startg) + '  sec')
    print('==========')
    return model_data






# Obs_data = pd.read_csv(r"/Users/charmibhatt/Library/CloudStorage/OneDrive-TheUniversityofWesternOntario/UWO_onedrive/Local_GitHub/edibles/edibles/utils/simulations/Charmi/fitting methods/5780_fitting/DIB5780_HD185418.txt", delim_whitespace=(True))
# #Obs_data['Wavelength']  = Obs_data[:,0] 
# #Obs_data['Wavelength'] = (1/Obs_data['Wavelength'])*1e8
# #Obs_data = Obs_data.iloc[::-1].reset_index(drop=True) #making it ascending order as we transformed wavelength into wavenumbers

# print(Obs_data['Flux'] )

# #shifting to zero and scaling flux between 0.9 and 1
# min_index = np.argmin(Obs_data['Flux'])

# Obs_data['Wavelength'] = Obs_data['Wavelength'] - Obs_data['Wavelength'][min_index]
# Obs_data['Flux']=  (Obs_data['Flux'] - min(Obs_data['Flux'])) / (1 - min(Obs_data['Flux'])) * 0.1 + 0.9

# plt.plot(Obs_data['Wavelength'], Obs_data['Flux'])
# #removing red wing
# #Obs_data_trp = Obs_data [(Obs_data['Wavelength'] >= -1) & (Obs_data['Wavelength']<= 2)]
# Obs_data_trp = Obs_data [(Obs_data['Flux'] <= 0.95)] #trp = triple peak structure


# #making data evenly spaced
# x_equal_spacing = np.linspace(min(Obs_data_trp['Wavelength']), max(Obs_data_trp['Wavelength']), len(Obs_data_trp['Wavelength']))
# y_obs_data = np.interp(x_equal_spacing, Obs_data_trp['Wavelength'], Obs_data_trp['Flux'])

#plt.plot(x_equal_spacing, y_obs_data)




def fit_model(B, T, delta_B, zeta, sigma, origin):
    mod = Model(get_rotational_spectrum) #, independent_vars = ['b', 'T']) #make sure independent variable of fitting function (that you made) is labelled as x
    params = mod.make_params( B =B, T=T, delta_B = delta_B, zeta = zeta, sigma = sigma, origin = origin)
    
    params['B'].min = 0.00005 
    params['B'].max = 0.05
    params['T'].min = 2.7
    params['T'].max = 1000
    # params['origin'].min = -2
    # params['origin'].max = 2
    # params['delta_B'].min = -2
    # params['delta_B'].max = 2
    params['zeta'].min = -1
    params['zeta'].max = 1
    # params['sigma'].min = 0.05
    # params['sigma'].max = 0.5

    result = mod.fit(y_obs_data, params, xx= Obs_data['Wavelength'], weights = 1/0.004) #method = 'leastsq', fit_kws={'ftol': 1e-12, 'xtol': 1e-12}
    print(result.fit_report())
    return result


get_rotational_spectrum( B = 5.1e-5, T = 479.13, delta_B = -1.83, zeta = -0.25, sigma = 0.86, origin =  -0.13)
# get_rotational_spectrum(x_equal_spacing, B = 0.00276, T = 78.23, delta_B = -0.07, zeta = -0.36, sigma = 0.1917, origin =  0.0047)

#result1= fit_model(B = 0.0023, T = 100, delta_B = -0.0353, zeta = -0.4197, sigma = 0.2247, origin =  0.0041)
  
#result1= fit_model(B = 0.01, T = 2.7, delta_B = -0.1, zeta = -0.1, sigma = 0.02, origin =  0)
# result2= fit_model(B = 0.005, T = 32.7, delta_B = -0.4, zeta = -0.9, sigma = 0.1, origin =  0.039)
# result3= fit_model(B = 0.002, T = 62.5, delta_B = -0.8, zeta = -0.5, sigma = 0.3, origin =  0.072)
# result4= fit_model(B = 0.0003, T = 32.5, delta_B = -0.45, zeta = -0.01, sigma = 0.17, origin =  0.012)
# result5= fit_model(B = 0.0075, T = 92.5, delta_B = -0.23, zeta = -0.23, sigma = 0.23, origin =  0.02)

# results_list = [result1, result2, result3, result4, result5]

def write_results_to_csv(results_list, filename):
    with open(filename, 'w', newline='') as csvfile:
        fieldnames = ['result_name', 'B_init', 'T_init', 'delta_B_init', 'zeta_init' , 'sigma_init' ,'origin_init' ,  'B', 'T',  'delta_B', 'zeta', 'sigma', 'origin', 'chi2', 'redchi', 'func_evals', 'B_unc', 'T_unc', 'delta_B_unc', 'zeta_unc', 'sigma_unc', 'origin_unc']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for i, result in enumerate(results_list):
            result_name = f'result{i+1}'
            params = result.params
            row = {
                'result_name': result_name,
                'B_init': params['B'].init_value,
                'T_init': params['T'].init_value,
                'delta_B_init': params['delta_B'].init_value,
                'zeta_init': params['zeta'].init_value,
                'sigma_init': params['sigma'].init_value,
                'origin_init': params['origin'].init_value,
                'B': params['B'].value,
                'T': params['T'].value,
                'delta_B': params['delta_B'].value,
                'zeta': params['zeta'].value,
                'sigma': params['sigma'].value,
                'origin': params['origin'].value,
                'chi2': result.chisqr,
                'redchi': result.redchi,
                'func_evals': result.nfev,
                'B_unc': params['B'].stderr,
                'T_unc': params['T'].stderr,
                'delta_B_unc': params['delta_B'].stderr,
                'zeta_unc': params['zeta'].stderr,
                'sigma_unc': params['sigma'].stderr,
                'origin_unc': params['origin'].stderr,

            }
            writer.writerow(row)
            
            
            
#write_results_to_csv(results_list, 'results_5780_HD185418.csv')     
            

# '''Lm fit'''

# mod = Model(get_rotational_spectrum) #, independent_vars = ['b', 'T']) #make sure independent variable of fitting function (that you made) is labelled as x
# #params = mod.guess(flux_data, x = np.linspace(0.005,0.01,5))
# params = mod.make_params(verbose = True, B = 0.01, T = 2.7, origin = 0, delta_B = -0.2, zeta = -0.5, sigma = 0.1953)



# params['B'].min = 0.0005 
# params['B'].max = 0.01
# params['T'].min = 2.7
# params['T'].max = 100
# params['origin'].min = -2
# params['origin'].max = 2
# params['delta_B'].min = -1
# params['delta_B'].max = 0
# params['zeta'].min = -1
# params['zeta'].max = 1
# params['sigma'].min = 0.05
# params['sigma'].max = 0.2

# result = mod.fit(y_obs_data, params, xx= x_equal_spacing, weights = 1/0.004, xtol = 1e-2, ftol = 1e-2) #, b = 0.005, T = 3, weights = 1/0.7)


# #plt.plot(x_obs_data, y_obs_data, label = 'Data')
# plt.plot(x_equal_spacing, result.best_fit, label = 'Fit')
# plt.legend()

# print(mod.param_names, mod.independent_vars)
# print(result.fit_report())

# minimizer = lmfit.Minimizer(mod, params)

# result2 = mod.fit(y_obs_data, params, xx= x_equal_spacing, weights = 1/0.004, method = 'nelder') #, b = 0.005, T = 3, weights = 1/0.7)
# print(result2.fit_report())
