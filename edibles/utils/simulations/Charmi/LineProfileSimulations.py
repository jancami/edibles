#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 16:22:21 2023

@author: charmibhatt
"""

import pandas as pd
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
import uncertainties as unc
import uncertainties.umath as umath 
import numba as nb
from pathlib import Path
from lmfit import Parameters      

Jmax = 300

def allowed_perperndicular_transitions(Jmax):
    
    'Take in Jmax and calculates the all allowed transitions. Here Jmax = Kmax'

    '''P Branch'''
    P_branch_Js = list(range(1, Jmax + 1))
    all_P_branch_Js = [j for j in P_branch_Js for _ in range(j)]
    P_branch_Jprimes = [j - 1 for j in all_P_branch_Js if j != 0]
    pP_Branch_K = [j - i for j in P_branch_Js for i in range(j)]
    pP_Branch_Kprime = [k - 1 for k in pP_Branch_K]
    rP_Branch_K = [i for j in P_branch_Js for i in sorted(range(0, j), reverse=True)]
    rP_Branch_Kprime = [k + 1 for k in rP_Branch_K]
    
    '''Q Branch'''
    Q_branch_Js = list(range(0, Jmax + 1))
    all_Q_branch_Js = [j for j in Q_branch_Js if j != 0 for _ in range(j)]
    Q_branch_Jprimes = all_Q_branch_Js[:]
    pQ_Branch_K = [j - i for j in Q_branch_Js for i in range(j)]
    pQ_Branch_Kprime = [k - 1 for k in pQ_Branch_K]
    rQ_Branch_K = [i for j in Q_branch_Js for i in sorted(range(0, j), reverse=True)]
    rQ_Branch_Kprime = [k + 1 for k in rQ_Branch_K]
    
    '''R Branch'''
    R_branch_Js = list(range(0, Jmax))
    all_R_branch_Js = [j for j in R_branch_Js if j == 0 or j != 0 for _ in range(j + 1)]
    R_branch_Jprimes = [j + 1 for j in all_R_branch_Js if j <= Jmax - 1]
    pR_Branch_K = [j - (i - 1) for j in R_branch_Js for i in range(j + 1)]
    pR_Branch_Kprime = [k - 1 for k in pR_Branch_K]
    rR_Branch_K = [i for j in R_branch_Js for i in sorted(range(0, j + 1), reverse=True)]
    rR_Branch_Kprime = [k + 1 for k in rR_Branch_K]
    
    '''Combine results'''
    Allowed_Js = (all_P_branch_Js * 2) + (all_Q_branch_Js * 2) + (all_R_branch_Js * 2)
    Allowed_Jprimes = (P_branch_Jprimes * 2) + (Q_branch_Jprimes * 2) + (R_branch_Jprimes * 2)
    Allowed_Ks = pP_Branch_K + rP_Branch_K + pQ_Branch_K + rQ_Branch_K + pR_Branch_K + rR_Branch_K
    Allowed_Kprimes = pP_Branch_Kprime + rP_Branch_Kprime + pQ_Branch_Kprime + rQ_Branch_Kprime + pR_Branch_Kprime + rR_Branch_Kprime
    
    columns = {'ground_J': Allowed_Js, 'excited_J': Allowed_Jprimes, 'ground_K': Allowed_Ks, 'excited_K': Allowed_Kprimes}
    combinations = pd.DataFrame(columns)
  
    
    combinations['delta_J'] = combinations['excited_J'] - combinations['ground_J']
    combinations['delta_K'] = combinations['excited_K'] - combinations['ground_K']
    
    
    delta_J_values = combinations['delta_J']
    delta_K_values = combinations['delta_K']
    
    label = [
        'pP' if delta_J == -1 and delta_K == -1
        else 'rP' if delta_J == -1 and delta_K == 1
        else 'pQ' if delta_J == 0 and delta_K == -1
        else 'rQ' if delta_J == 0 and delta_K == 1
        else 'pR' if delta_J == 1 and delta_K == -1
        else 'rR'
        for delta_J, delta_K in zip(delta_J_values, delta_K_values)
    ]
    
    combinations['label'] = label
        
    return combinations


def get_rotational_spectrum(B, delta_B, zeta, T, sigma, origin):
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

    # endl = timeit.default_timer()
    # print('>>>> linelist calculation takes   ' + str(endl - startl) + '  sec')

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
    
    model_data = np.array([simu_waveno, simu_intenisty]).transpose()


    # for units in wavelength
    # simu_wavelength = (1/simu_waveno)*1e8
    # model_data = np.array([simu_wavelength, simu_intenisty]).transpose()
    # model_data = model_data[::-1]

    y_model_data = model_data[:,1]
    endg = timeit.default_timer()

    print('>>>> Time taken to simulate thi profile  ' + str(endg - startg) + '  sec')
    print('==========')
    return y_model_data 


# B = 0.00445875
# delta_B = -0.11518035
# zeta = -0.14142679 
# T = 105.221344
# sigma = 0.16178154 
# origin = 0.06567328

#linelist, model_data = get_rotational_spectrum(B, delta_B, zeta, T, sigma, origin)


def curve_to_fit_wavenos(sightline): 
        spec_dir = Path("/Users/charmibhatt/Library/CloudStorage/OneDrive-TheUniversityofWesternOntario/UWO_onedrive/Local_GitHub/edibles/edibles/utils/simulations/Charmi/Heather's_data/")

        file = '6614_HD{}.txt'.format(sightline)
        Obs_data = pd.read_csv(spec_dir / file,
                                sep=',')
        Obs_data['Wavelength'] = (1 / Obs_data['Wavelength']) * 1e8
        Obs_data = Obs_data.iloc[::-1].reset_index(
            drop=True)  # making it ascending order as we transformed wavelength into wavenumbers

        # shifting to zero and scaling flux between 0.9 and 1
        min_index = np.argmin(Obs_data['Flux'])
        Obs_data['Wavelength'] = Obs_data['Wavelength'] - Obs_data['Wavelength'][min_index] 
        Obs_data['Flux'] = (Obs_data['Flux'] - min(Obs_data['Flux'])) / (1 - min(Obs_data['Flux'])) * 0.1 + 0.9
        x_obs_data = Obs_data['Wavelength']
        # removing red wing
        # Obs_data_trp = Obs_data [(Obs_data['Wavelength'] >= -1) & (Obs_data['Wavelength']<= 1.2)]
        Obs_data_trp = Obs_data[(Obs_data['Flux'] <= 0.95)]  # trp = triple peak structure

        # making data evenly spaced
        x_equal_spacing = np.linspace(min(Obs_data_trp['Wavelength']), max(Obs_data_trp['Wavelength']), 25)
        y_obs_data = np.interp(x_equal_spacing, Obs_data_trp['Wavelength'], Obs_data_trp['Flux'])

        Obs_data_continuum = Obs_data [(Obs_data['Wavelength'] >= 2) & (Obs_data['Wavelength']<= 5)]
        std_dev = np.std(Obs_data_continuum['Flux'])
        
        return x_obs_data, y_obs_data, std_dev
    
curve_to_fit_wavenos('166937')
    
sightline = '166937'
def fit_model(B, delta_B, zeta, T, sigma, origin):
    mod = Model(get_rotational_spectrum) #, independent_vars = ['b', 'T']) #make sure independent variable of fitting function (that you made) is labelled as x
    params = mod.make_params( B =B, delta_B = delta_B, zeta = zeta, T=T,sigma = sigma, origin = origin)
    

    x_obs_data, y_obs_data, std_dev = curve_to_fit_wavenos(sightline)
    #print(std_dev)
    result = mod.fit(y_obs_data, params) #, x_obs_data) #, weights = 1/std_dev)
    print(result.fit_report())
    return result

result = fit_model(B = 0.005,  delta_B = -0.4, zeta = -0.9, T = 32.7, sigma = 0.1, origin =  0.039)

