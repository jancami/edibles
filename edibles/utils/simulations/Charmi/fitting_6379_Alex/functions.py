# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 12:29:44 2023

@author: alexr

All the functions used in the analysis of the 6379 DIB using rovribrational spectroscopy. Modified from original code written by Charmi Bhatt. 

WARNING: Some functions make use of the 'beepy' module. This is a module that has a number of pre-loaded sound effects, that
I use to alert me at times when, for example, the code has finished running. It relies on having the module simpleaudio installed, 
so if they are not found many of these functions will produce a related error message! These are not essential so can be commented out 
or deleted, however if you want to use them they are easily installed through pip! 

Fun note: bp.beep(sound = 7) or bp.beep(sound = 'wilhelm') is a sound effect of man screaming. 
"""

#%% Imports

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import timeit
import astropy.constants as const
import numba as nb
from lmfit import Model
import beepy as bp

# import cProfile

#%% Making equal spaced grid of wavenumbers according to edibles resolution

def make_grid(lambda_start, lambda_end, resolution=None, oversample=None):
    ''' 
    Makes a grid of equal spaced wavenumbers to fit to the interpolation of observational data, which can then be matched to the fitting of a model. 
    '''
    # check keywords
    if oversample is None:
        oversample = 40.0
    if resolution is None:
        resolution = 1500.0

    lambda_start = np.float64(lambda_start)
    lambda_end = np.float64(lambda_end)

    # produce grid
    R = resolution * oversample
    
    # print('R = ' , R)
    n_points = (
        round(
            (np.log(lambda_end / lambda_start)) / (np.log(-(1 + 2 * R) / (1 - 2 * R)))
        )
        + 1
    )
    # print('n_points = ' , n_points)
    f = -(1 + 2 * R) / (1 - 2 * R)
    
    # print('f = ', f)
    factor = f ** np.arange(n_points)
    # print('factor = ' , factor)
    wave = np.full(int(n_points), lambda_start, dtype=np.float64)
    # print('wave = ' , wave)
    grid = wave * factor
    # print('grid = ', grid)
    return grid

#%% Extracting observational data to plot

def curve_to_fit_wavenos(sightline):
    '''
    This has been renamed to obs_curve_to_plot!!
    '''
    print('Replace curve_to_fit_wavenos function with obs_curve_to_plot!')
    bp.beep(sound = 'robot_error')

def obs_curve_to_plot(sightline, wavenos = True, scaled = True, zero = None): 
    
    '''
    RENAMED from curve_to_fit_wavenos 
    
    Takes in observations of a 6379 dib from Heather MacIsaac's data and returns data in a form that can be plotted. The path to find the file is defined in the function.
    
    Args: 
        sightline (str): star identifier to fill in file name '6379_HD{}_avg_spectra.csv'
        
        wavenos (bool): if true (default) returns the data in units of wavenumber (cm^-1) centred around zero, else in Angstroms at true value
        
        scaled (bool): if true (default) scales the flux between 0.9 and 1
        
        zero (string): Sets the point to which the plot is aligned to. Allowed values: 
            
            None (default): Returns the original unshifted data.
            
            min: if wavenos is true, data is shifted so 0 cm^-1 is at the minimum flux value of the plot.
            
            6379: If wavenos is false, minimum flux value is at 6379.
            
            Any other value will give the original data
        
    Returns:
        x_obs_data (pandas series): Wavenumber or wavelength values
        
        y_obs_data (pandas series): Flux values
        
        std_dev (float): Standard deviation of the continuum
        
        x_axis (str): appropriate x-axis label for a plot of the data (waveno or wavelength) 
    '''
    
    spec_dir = Path(r"C:\Users\alexr\edibles\edibles\utils\simulations\Charmi\fitting_6379_Alex\Heather_MacIsaac_6379_data")
    file = '6379_HD{}_avg_spectra.csv'.format(sightline)
    Obs_data = pd.read_csv(spec_dir / file,
                                sep=',')
    x_axis = 'Angstrom / $\AA$'
    
    if wavenos == True:
        Obs_data['Wavelength'] = (1 / Obs_data['Wavelength']) * 1e8
        Obs_data = Obs_data.iloc[::-1].reset_index(
            drop=True)  # making it ascending order as we transformed wavelength into wavenumbers
        x_axis = 'Wavenumber / cm$^{-1}$'
    
        
    if wavenos == True: 
        if zero == 'min':
            min_index = np.argmin(Obs_data['Flux'])
            Obs_data['Wavelength'] = Obs_data['Wavelength'] - Obs_data['Wavelength'][min_index] 
        else:
            Obs_data['Wavelength'] = Obs_data['Wavelength'] - (1/6379)*1e8
    else:
        if zero == '6379':
            min_index = np.argmin(Obs_data['Flux'])
            Obs_data['Wavelength'] = Obs_data['Wavelength'] - Obs_data['Wavelength'][min_index] +6379
        else:
            Obs_data['Wavelength'] = Obs_data['Wavelength']
   
    if scaled == True:  
        # shifting to zero and scaling flux between 0.9 and 1
        Obs_data['Flux'] = (Obs_data['Flux'] - min(Obs_data['Flux'])) / (1 - min(Obs_data['Flux'])) * 0.1 + 0.9
    y_obs_data = Obs_data['Flux']
    x_obs_data = Obs_data['Wavelength']

    Obs_data_continuum = Obs_data [(Obs_data['Wavelength'] >= 2) & (Obs_data['Wavelength']<= 5)]
    std_dev = np.std(Obs_data_continuum['Flux'])
        
    return x_obs_data, y_obs_data, std_dev, x_axis

    
# # obs_curve_to_plot test:

if __name__ == "__main__":
    sightline = '185859'

    x_obs, y_obs, std, x_label = obs_curve_to_plot(sightline, wavenos= False, zero ='6379')
    # print(len(x_obs))

    fig, ax = plt.subplots()
    ax.plot(x_obs, y_obs, label = 'HD{}'.format(sightline))
    ax.xaxis.set_major_locator(plt.MultipleLocator(1))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(0.5))
    ax.set_xlabel(x_label)
    ax.set_ylabel('Flux')
    ax.legend()
    plt.show()

    
#%% Extracting observational data to fit

def obs_curve_to_fit(sightline, fitrange = '0.95'): 
    
    '''
    Takes in observations of a 6379 dib from Heather MacIsaac's data and returns data in a form that can be passed to fit_model in order to give a best fit model. 
    Not to be confused with obs_curve_to_plot which returns data in a different form which can be plotted! 
    This function focusses on the peak structure and interpolates according to the resolution of the data to create a basis for modelling. 
    The zero point of the resulting data set is at 6379 Angstroms (15676.4 cm^-1)
    
    Args: 
        sightline (str): star identifier to fill in file name '6379_HD{}_avg_spectra.csv'
        
    Returns:
        Obs_data (pandas DataFrame): Data frame of the observational data from the file, columns 'Wavelength' and 'Flux'.
        
        x_equal_spacing (numpy array): equal spaced wavenumber values (according to the data's resolution) in the region of the spectrum with the absorption feature
        
        y_data_fit (numpy array): Corresponding flux values to x_equal_spacing, found through interpolation of the spectrum
        
        std_dev (float): standard deviation of the continuum
    '''
    
    spec_dir = Path(r"C:\Users\alexr\edibles\edibles\utils\simulations\Charmi\fitting_6379_Alex\Heather_MacIsaac_6379_data")
    file = '6379_HD{}_avg_spectra.csv'.format(sightline)
    Obs_data = pd.read_csv(spec_dir / file,
                                sep=',', )
    Obs_data.columns.values[0] = 'Initial index'

    # shifting to zero and scaling flux between 0.9 and 1
    Obs_data['Flux'] = (Obs_data['Flux'] - min(Obs_data['Flux'])) / (1 - min(Obs_data['Flux'])) * 0.1 + 0.9
    
    if fitrange == '0.95':
        # Old lambda start and end for fitting only flux <= 0.95
        Obs_data_trp = Obs_data[(Obs_data['Flux'] <= 0.95)]  # trp = triple peak structure, focussing on the peak structure
        lambda_start = min(Obs_data_trp['Wavelength']) # Start and end of the section of the spectra we are focussing on
        lambda_end = max(Obs_data_trp['Wavelength'])
    
    Obs_data['Wavelength'] = (1 / Obs_data['Wavelength']) * 1e8 # Converting to waveno.
    Obs_data = Obs_data.iloc[::-1].reset_index(
        drop=True)  # making it ascending order as we transformed wavelength into wavenumbers
    # min_index = np.argmin(Obs_data['Flux'])  ## In case the zero point needs to be changed to the minimum of the peak
    Obs_data['Wavelength'] = Obs_data['Wavelength'] - (1/6379)*1e8  #Obs_data['Wavelength'][min_index] ## Zero point at 6379A in cm^-1
    Obs_data_trp = Obs_data[(Obs_data['Flux'] <= 0.95)] # Redefine the triple peak section, now in cm^-1
    
    if fitrange == '0.95':
        Obs_data_for_fit = Obs_data_trp
    
    elif fitrange == 'bluewing':
        # Defining wavenumber range for fitting flux <= 0.95 and blue wing
        fit_index_zero = Obs_data_trp.iloc[0].name
        fit_index_end = fit_index_zero + 20
        Obs_data_for_fit = Obs_data.iloc[fit_index_zero:fit_index_end]
        # print(len(Obs_data_for_fit))
        lambda_end = 1/(Obs_data.iloc[fit_index_zero]['Wavelength'] + (1/6379)*1e8) * 1e8 # min(Obs_data_trp['Wavelength']) # Start and end of the section of the spectra we are focussing on
        lambda_start = 1/(Obs_data.iloc[fit_index_end]['Wavelength'] +(1/6379)*1e8) * 1e8 #max(Obs_data_trp['Wavelength'])
        # print(lambda_start), print(lambda_end)
        # print(Obs_data_for_fit), print(Obs_data_trp)
    
    # making data evenly spaced
    common_grid_for_all = make_grid(lambda_start, lambda_end, resolution=107000, oversample=2)
    common_grid_for_all = (1 / common_grid_for_all) * 1e8
    common_grid_for_all = common_grid_for_all - (1/6379)*1e8 # Zero point is at 6379A in cm^-1
    common_grid_for_all = common_grid_for_all[::-1] # Reverse direction
    x_equal_spacing = common_grid_for_all
    # x_equal_spacing = np.linspace(min(Obs_data_trp['Wavelength']), max(Obs_data_trp['Wavelength']), 100) ## Old interpolation was 100 points between max and min
    y_data_fit = np.interp(x_equal_spacing, Obs_data_for_fit['Wavelength'], Obs_data_for_fit['Flux']) # Interpolation of flux values

    Obs_data_continuum = Obs_data [(Obs_data['Wavelength'] >= 2) & (Obs_data['Wavelength']<= 5)]
    std_dev = np.std(Obs_data_continuum['Flux']) #Standard deviation of the continuum
        
    return Obs_data, x_equal_spacing, y_data_fit, std_dev


# # obs_curve_to_fit test:
if __name__ == "__main__":
    sightline = '185859'
    
    Obs_data, x_equal_spacing, y_data_fit, std_dev = obs_curve_to_fit(sightline, fitrange='yuu')

    fig, ax = plt.subplots()
    # ax.plot(x_equal_spacing, y_data_fit, label = 'HD{}'.format(sightline), marker = 'o')
    ax.plot(Obs_data['Wavelength'], Obs_data['Flux'], color = 'r', label = 'Full, HD{}'.format(sightline),  linewidth = '2.5')
    ax.plot(x_equal_spacing, y_data_fit, color = 'cyan', label = 'Trp', linestyle = '-.')
    # ax.xaxis.set_major_locator(plt.MultipleLocator(1))
    # ax.xaxis.set_minor_locator(plt.MultipleLocator(0.5))
    ax.set_xlabel('Waveno')
    ax.set_ylabel('Flux')
    ax.legend()
    plt.show()

#%% Allowed Perp Transitions Calculation

def allowed_perperndicular_transitions(Jmax):
    
    '''
    Take in Jmax and calculates the all allowed perpendicular transitions. Here Jmax = Kmax.
    
    Args: 
        Jmax(int).
    
    Returns: 
        combinations (pandas Dataframe): allowed transitions.
    '''

    # start = timeit.default_timer()

    # P Branch
    P_branch_Js = list(range(1, Jmax + 1))
    all_P_branch_Js = [j for j in P_branch_Js for _ in range(j)]
    P_branch_Jprimes = [j - 1 for j in all_P_branch_Js if j != 0]
    pP_Branch_K = [j - i for j in P_branch_Js for i in range(j)]
    pP_Branch_Kprime = [k - 1 for k in pP_Branch_K]
    rP_Branch_K = [i for j in P_branch_Js for i in sorted(range(0, j), reverse=True)]
    rP_Branch_Kprime = [k + 1 for k in rP_Branch_K]
    
    # Q Branch
    Q_branch_Js = list(range(0, Jmax + 1))
    all_Q_branch_Js = [j for j in Q_branch_Js if j != 0 for _ in range(j)]
    Q_branch_Jprimes = all_Q_branch_Js[:]
    pQ_Branch_K = [j - i for j in Q_branch_Js for i in range(j)]
    pQ_Branch_Kprime = [k - 1 for k in pQ_Branch_K]
    rQ_Branch_K = [i for j in Q_branch_Js for i in sorted(range(0, j), reverse=True)]
    rQ_Branch_Kprime = [k + 1 for k in rQ_Branch_K]
    
    # R Branch
    R_branch_Js = list(range(0, Jmax))
    all_R_branch_Js = [j for j in R_branch_Js if j == 0 or j != 0 for _ in range(j + 1)]
    R_branch_Jprimes = [j + 1 for j in all_R_branch_Js if j <= Jmax - 1]
    pR_Branch_K = [j - (i - 1) for j in R_branch_Js for i in range(j + 1)]
    pR_Branch_Kprime = [k - 1 for k in pR_Branch_K]
    rR_Branch_K = [i for j in R_branch_Js for i in sorted(range(0, j + 1), reverse=True)]
    rR_Branch_Kprime = [k + 1 for k in rR_Branch_K]
    
    # Combine results
    Allowed_Js = (all_P_branch_Js * 2) + (all_Q_branch_Js * 2) + (all_R_branch_Js * 2)
    Allowed_Jprimes = (P_branch_Jprimes * 2) + (Q_branch_Jprimes * 2) + (R_branch_Jprimes * 2)
    Allowed_Ks = pP_Branch_K + rP_Branch_K + pQ_Branch_K + rQ_Branch_K + pR_Branch_K + rR_Branch_K
    Allowed_Kprimes = pP_Branch_Kprime + rP_Branch_Kprime + pQ_Branch_Kprime + rQ_Branch_Kprime + pR_Branch_Kprime + rR_Branch_Kprime
    
    # Put results into a pandas dataframe
    columns = {'ground_J': Allowed_Js, 'excited_J': Allowed_Jprimes, 
               'ground_K': Allowed_Ks, 'excited_K': Allowed_Kprimes}
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
        
    # end = timeit.default_timer()
    # time = end - start
    # print('computation took ' + str(time) + 's')
    
    return combinations

# combinations  = allowed_perperndicular_transitions(300)

#%% Allowed Parallel Transitions calculation

def allowed_parallel_transitions(Jmax):
    
    '''
    Take in Jmax and calculates the all allowed Parallel transitions. Here Jmax = Kmax.
    
    Args: 
        Jmax(int): Maximum allowed angular momentum quantum number
    
    Returns: 
        combinations (pandas Dataframe): allowed transitions.
    '''
    start = timeit.default_timer()
    
    # P Branch
    P_branch_Js = list(range(1, Jmax + 1))
    all_P_branch_Js = [j for j in P_branch_Js for _ in range(j)] #Give a list of all starting Js, each j is repeated j times to allow for all the Ks
    P_branch_Jprimes = [j - 1 for j in all_P_branch_Js if j != 0] #List of all the excited Js, delta J = -1 for P branch
    qP_Branch_K = [j - i for j in P_branch_Js for i in range(j)] #List of all the available Ks, going from j down to zero, each corresponding to a j in all_P_branch_Js 
    qP_Branch_Kprime = [k for k in qP_Branch_K] #Excited Ks are equal to ground Ks by definition of parallel transition

    
    # Q Branch
    Q_branch_Js = list(range(0, Jmax + 1))
    all_Q_branch_Js = [j for j in Q_branch_Js if j != 0 for _ in range(j)]
    Q_branch_Jprimes = all_Q_branch_Js[:] #Delta J = 0 for Q branch
    qQ_Branch_K = [j - i for j in Q_branch_Js for i in range(j)]
    qQ_Branch_Kprime = [k for k in qQ_Branch_K]

    # R Branch
    R_branch_Js = list(range(0, Jmax))
    # all_R_branch_Js = [j for j in R_branch_Js if j == 0 or j != 0 for _ in range(j + 1)]
    # all_R_branch_Js = [j if j == 0 else j for j _ in range(j) for j in R_branch_Js]
    # R_branch_Js = list((range(0,Jmax)))
    all_R_branch_Js = []
    for j in R_branch_Js:
        if j ==0:
            all_R_branch_Js.append(j)   #Include each j j times unless = 0, in which case include once
        elif j!= 0:
            for i in range(j+1):
              all_R_branch_Js.append(j) #Written as a loop rather than a list comprehension for readability
    R_branch_Jprimes = [j + 1 for j in all_R_branch_Js if j <= Jmax - 1] #Delta J = +1 for R branch
    qR_Branch_K = [j - (i - 1) for j in R_branch_Js for i in range(j + 1)]
    qR_Branch_Kprime = [k for k in qR_Branch_K]

    #Combining all three branches together
    Allowed_Js = (all_P_branch_Js) + (all_Q_branch_Js) + (all_R_branch_Js)
    Allowed_Jprimes = (P_branch_Jprimes) + (Q_branch_Jprimes) + (R_branch_Jprimes)
    Allowed_Ks = qP_Branch_K +  qQ_Branch_K +   qR_Branch_K 
    Allowed_Kprimes = qP_Branch_Kprime  + qQ_Branch_Kprime  + qR_Branch_Kprime 

    #Putting results in DataFrame
    columns = {'ground_J' : Allowed_Js,'excited_J': Allowed_Jprimes, 
               'ground_K' : Allowed_Ks, 'excited_K' : Allowed_Kprimes}
    combinations = pd.DataFrame(columns)
  
    #Calculating delta values
    combinations['delta_J'] = combinations['excited_J'] - combinations['ground_J']
    combinations['delta_K'] = combinations['excited_K'] - combinations['ground_K']
    
    
    delta_J_values = combinations['delta_J']
    delta_K_values = combinations['delta_K']
    
    label = [
        'qP' if delta_J == -1 and delta_K == 0
        else 'qQ' if delta_J == 0 and delta_K == 0
        else 'qR'
        for delta_J, delta_K in zip(delta_J_values, delta_K_values)
    ]
    
    combinations['label'] = label
    end = timeit.default_timer()
    time = end - start
    print('computation took ' + str(time) + 's')
    return combinations
# combinations = allowed_parallel_transitions(3)
# print(combinations)
#%% Generate Rotational Spectrum
# 
def get_rotational_spectrum(B, delta_B, zeta, T, sigma, origin, combinations, transition = 'perpendicular', Jmax = 300, bell = True):
    
    '''
    Generates a model of the rotational spectrum of a molecule based on input parameters
    
    - Requires allowed_perperndicular_transitions to be defined.
    
    - Prints the time taken to compute the model
    
    Args (floats): 
        B (cm^-1): Ground rotational constant of the molecule. Assume 3D oblate symmetric top A = B = 2C.
        
        delta_B (%): Change in B to excited state. Assume delta_B = delta_C.
        
        zeta (cm^-1): Coriolis constant.
        
        T (K): Rotational Temperature.
        
        sigma: std of Gaussian smoothing.
        
        origin: allows for a shift along the x-axis.
        
        Jmax: Default 300, for passing to allowed_perperndicular_transitions function. (redundant now this has been removed from the fn)
        
        combinations (pandas Dataframe): Result of allowed transitions function
        
        bell (bool): default True, sounds a notification bell once the computation is complete.
        
        transition (string): Type of transition, allowed types are 'perpendicular' (default) and 'parallel'
    
    Returns:
        x_model_data, y_model_data (numpy array): calculated wavenumber and flux values for the model
        
    Notes: 
        Error message - TypeError: unsupported operand type(s) for *: 'NoneType' and 'float' in line 568, 
        
        >>>> strength = (HL_factors[i] * BD_factors[i])
        
        Ensure the transition type matches between the combinations input and the HL values calculation as if they are different, these arrays are different sizes
    '''
    
    start1 = timeit.default_timer()


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
    # 
    # end1 = timeit.default_timer()
    # print('Time to import parameters ' + str(end1-start1))
    
    # Calculating Linelist

    # start2 = timeit.default_timer()
    
    # ground_Es = []
    # excited_Es = []
    # wavenos = []
    def calculate_ground_Es(Js, Ks, B, C):
     return B * Js * (Js + 1) + (C - B) * (Ks ** 2)

    def calculate_excited_Es(Js, Ks, del_K, B, C, zeta):
        base_Es = B * Js * (Js + 1) + (C - B) * (Ks ** 2) + C ** 2 
        zeta_component = ((-2 * C * zeta)) * Ks 
        return np.where(del_K == -1, base_Es - zeta_component, base_Es + zeta_component)
     
    def calculate_wavenos(origin, excited_Es, ground_Es):
        return origin + excited_Es - ground_Es
     
     # Then we apply these functions
    linelist['ground_Es'] = calculate_ground_Es(ground_Js, ground_Ks, ground_B, ground_C)
    linelist['excited_Es'] = calculate_excited_Es(excited_Js, excited_Ks, delta_K, excited_B, excited_C, zeta)
    linelist['wavenos'] = calculate_wavenos(origin, linelist['excited_Es'], linelist['ground_Es'])

    ground_Es = linelist['ground_Es'] 

    # end2 = timeit.default_timer()
    # print('Time to make linelist '+str(end2 - start2))
    #print('Overall time to calculate linelist ' + str(end2-start2))

    # start3 = timeit.default_timer()
    
    #Calculating Honl-London factors to determine relative intensities

    HL_factors = []
    
    # Compute HL_factors directly in a list comprehension
    # different for if the transition is parallel or perpendicular
    if transition == 'perpendicular': # Added 2* on bottom of all fractions??, don't need to as relative intensities, but I have to compare to parallel
        HL_factors = [((J - 1 + K) * (J + K)) / (2*J * ((2 * J) + 1)) if (delta_J == -1 and delta_K == -1) else
                      ((J - 1 - K) * (J - K)) / (2*J * ((2 * J) + 1)) if (delta_J == -1 and delta_K == 1) else
                      (J + 1 - K) * (J + K) / (2*J * (J + 1)) if (delta_J == 0 and delta_K == -1) else
                      (J + 1 + K) * (J - K) / (2*J * (J + 1)) if (delta_J == 0 and delta_K == 1) else
                      (J + 2 - K) * (J + 1 - K) / (2*(J + 1) * ((2 * J) + 1)) if (delta_J == 1 and delta_K == -1) else
                      (J + 2 + K) * (J + 1 + K) / (2*(J + 1) * ((2 * J) + 1)) if (delta_J == 1 and delta_K == 1) else
                      None for J, K, delta_J, delta_K in zip(ground_Js, ground_Ks, delta_J, delta_K)]
    elif transition == 'parallel':
        HL_factors = [(J ** 2 - K ** 2) / (J * ((2 * J) + 1)) if (delta_J == -1 and delta_K == 0) else
                      (K ** 2) / (J * (J + 1)) if (delta_J == 0 and delta_K == 0) else
                      ((J + 1) ** 2  - K ** 2) / ((J + 1) * ((2 * J) + 1)) if (delta_J == 1 and delta_K == 0) else
                      None for J, K, delta_J, delta_K in zip(ground_Js, ground_Ks, delta_J, delta_K)]
    linelist['HL_factors'] = HL_factors

    # end3 = timeit.default_timer()
    # print('Time to calculate HL Factors ' + str(end3-start3))

    # Calculate populations of each level with Boltzmann eqn
    # start4 = timeit.default_timer()
    BD_factors = []

    h = const.h.cgs.value
    c = const.c.to('cm/s').value
    k = const.k_B.cgs.value

    # Conversion to numpy arrays
    ground_Js_np = np.array(ground_Js)
    ground_Ks_np = np.array(ground_Ks)
    ground_Es_np = np.array(ground_Es)
    
    # Calculation of static part
    static_part = (-h * c) / (k * T)
    
    # For the condition K == 0, the equation becomes twice the general equation
    # we can pre-calculate that part
    factor = (2 * ground_Js_np + 1) * np.exp(static_part * ground_Es_np)
    
    # Wherever ground_Ks == 0, double the result
    boltzmann_equation = np.where(ground_Ks_np == 0, 2 * factor, factor)
    
    # Convert the numpy array back to list if necessary
    BD_factors = boltzmann_equation.tolist()

    intensities = []
    for i in range(len(linelist.index)):
        strength = (HL_factors[i] * BD_factors[i])
        intensities.append(strength)

    linelist['intensities'] = intensities
    
    # end4 = timeit.default_timer()
    # print('Time to calculate BD Factors ' + str(end4-start4))

    # Smoothening the linelist using a Gaussian smoothing with std sigma and the numba decorator
    
    # start5 = timeit.default_timer()    
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
    smooth_intensities = calculate_smooth_intensities(Wavenos_arr, Intenisty_arr, smooth_wavenos, sigma)
    # end5 = timeit.default_timer()
    # print('Time for numba ' + str(end5-start5))
    

    # call the numba function with input data

    # start6 = timeit.default_timer()
    smooth_data = np.array([smooth_wavenos, smooth_intensities]).transpose()
    smooth_data = np.delete(smooth_data, np.where(smooth_data[:, 1] <= 0.001 * (max(smooth_data[:, 1]))), axis=0)

    simu_waveno = smooth_data[:, 0]
    simu_intenisty = 1 - 0.1 * (smooth_data[:, 1] / max(smooth_data[:, 1])) # Scale to between 0.9 and 1
    
    model_data = np.array([simu_waveno, simu_intenisty]).transpose()

    # end6 = timeit.default_timer()
    # print('Time to for convolution ' + str(end6-start6))

    # for units in wavelength
    # simu_wavelength = (1/simu_waveno)*1e8
    # model_data = np.array([simu_wavelength, simu_intenisty]).transpose()
    # model_data = model_data[::-1]

    y_model_data = model_data[:,1]
    x_model_data = model_data[:,0]
    endg = timeit.default_timer()
    
    # print('>>>> Time taken to smooth profile  ' + str(endg - linelist_time) + '  sec')
    print('>>>> Time taken to simulate this profile  ' + str(endg - start1) + '  sec')
    print('==========')
    if bell:
        print('\a')
    return x_model_data, y_model_data#, linelist

if __name__ == "__main__":
    Jmax_test = 320
    B = 0.0016
    delta_B = 0.2
    zeta = -0.35
    T = 70
    sigma = 0.1953
    origin = 0.22
    # combinations  = allowed_parallel_transitions(Jmax_test)
    combinations  = allowed_perperndicular_transitions(Jmax_test)

    xs, ys = get_rotational_spectrum(B, delta_B, zeta, T, sigma, origin, combinations, bell = 1, transition='perpendicular')
    
    # fig, ax = plt.subplots()
    # ax.plot(xs, ys, label = 'Model fit')
    # ax.xaxis.set_major_locator(plt.MultipleLocator(1))
    # ax.xaxis.set_minor_locator(plt.MultipleLocator(0.5))
    # ax.set_xlabel('Wavenumber / cm$^{-1}$')
    # ax.set_ylabel('Flux')
    # ax.legend()    
    
    
    # #Testing different Jmax values with just the linelist
    # fig, ax = plt.subplots(figsize=(32,20), facecolor = 'none')
    
    # for label in (ax.get_xticklabels() + ax.get_yticklabels()):
    # 	label.set_fontsize(32)
    
    # ax.plot(linelist['wavenos'],linelist['intensities'], color = 'forestgreen')
    # ax.set_xlabel('Wavenumber/$cm^{-1}$', fontsize = 32)
    # ax.set_ylabel('Intensity (Arbitrary Units)', fontsize = 32)
    # # plt.xlim((-0.2,4))
    # # plt.title('Jmax = ' + str(Jmax_test), fontsize = 40)
    # ax.annotate('P\n$\Delta J = -1$', (-0.6,190), xycoords='data', fontsize = 50)
    # ax.annotate('Q\n$\Delta J = 0$', (0.1,160), xycoords='data', fontsize = 50)
    # ax.annotate('R\n$\Delta J = +1$', (0.85,190), xycoords='data', fontsize = 50)
    
    # # plt.savefig('Jmax_tests_parallel/Jmax = {}_params2.png'.format(Jmax_test), bbox_inches = 'tight')
    # # plt.savefig("C:\\Users\\alexr\\OneDrive - Durham University\\GRI Mitacs Summer 23\\Project\\Presentation\\Figures\\PQR.png", bbox_inches = 'tight')
    
    # plt.show()



#%% Function to pass to fitting model

def model_curve_to_fit(x_equal_spacing, B, delta_B, zeta, T, sigma, origin, combinations, sightline, transition, Jmax):
    ''' 
    Function to pass into fit_model, generates model data with the appropriate dimensions
    '''
    x_model_data, y_model_data = get_rotational_spectrum(B, delta_B, zeta, T, sigma, origin, combinations, transition, Jmax, bell = False)   
    
    # x_obs_data, y_obs_data, std_dev, x_axis = obs_curve_to_fit(sightline)
    Obs_data, x_equal_spacing, y_data_fit, std_dev = obs_curve_to_fit(sightline,fitrange='bluewing')
    plt.plot(x_model_data, y_model_data, label = 'Model')
    # Obs_data = Obs_data[Obs_data['Flux']<=0.95]
    # plt.plot(Obs_data['Wavelength'], Obs_data['Flux'], label = 'Raw obs, HD{}'.format(sightline))
    plt.plot(x_equal_spacing, y_data_fit, label = 'Interpolated obs, HD{}'.format(sightline))
    plt.legend()
    plt.figsize = (1,1.5)
    plt.show()
    y_model_data = np.interp(x_equal_spacing, x_model_data, y_model_data)
    print(B)
    print(delta_B)
    print(zeta)
    print(T)
    print(sigma)
    print(origin)
    

    return y_model_data

#%% Model fitting
def fit_model(B, delta_B, zeta, T, sigma, origin, combinations, sightline, transition, Jmax):
    '''
    lmfit fitting of observational data through attempting to find a best fit. 
    Plots the result and makes a trumpet sound when complete

    Parameters
    ----------
    Parameters to pass to get_rotational_spectrum

    Returns
    -------
    result: ModelResult, Result of the fitting process

    '''
    start = timeit.default_timer()
    mod = Model(model_curve_to_fit, 
                independent_vars=['sightline','x_equal_spacing', 'combinations', 'transition', 'Jmax'], 
                param_names=['B','delta_B','zeta','T','sigma','origin']) 
    
    params = mod.make_params( B = B, delta_B = delta_B, zeta = zeta, T=T,sigma = sigma, origin = origin)
    params['B'].min = 0.0005 
    params['B'].max = 0.01
    params['T'].min = 2.7
    params['T'].max = 300
    params['origin'].min = -2
    params['origin'].max = 2
    params['delta_B'].min = -1
    params['delta_B'].max = 1
    params['zeta'].min = -3
    params['zeta'].max = 1
    params['sigma'].min = 0.05
    params['sigma'].max = 0.3
    
    
    print(sightline)
    print(params)
    Obs_data, x_equal_spacing, y_data_fit, std_dev = obs_curve_to_fit(sightline, fitrange='bluewing')
    print(std_dev)
    print(len(y_data_fit))
    print(len(x_equal_spacing))
    print(x_equal_spacing)
    print(y_data_fit)
    plt.plot(x_equal_spacing, y_data_fit, label = 'Observational fit data')
    plt.legend()    
    plt.show()
    print(Jmax)

    
    result = mod.fit(y_data_fit, params, weights = 1/std_dev, 
                     x_equal_spacing = x_equal_spacing, 
                     sightline=sightline, 
                     combinations = combinations, 
                     transition = transition, 
                     Jmax=Jmax
                     )
    print(result.fit_report())
    end = timeit.default_timer()
    print('Time taken to generate model ' + str(end - start))
    
    bp.beep(sound = 6) #Trumpet fanfare sound to say that it has finished!
    
    print('Sightline is HD{}'.format(sightline))
    print('Jmax = {}'.format(Jmax))
    
    # Calculate residuals between best fit model and observation
    def residual(y_obs_data, y_model_data, sigma):
        y_obs_data = np.array(y_obs_data)
        y_model_data = np.array(y_model_data)
        r = (y_obs_data - y_model_data)/sigma
        return r
    r = residual(y_data_fit, result.best_fit, std_dev)
    
    # Plot best fit model, observation and residuals on same axis
    def plot_best_fit(result, x_equal_spacing, y_obs_data):
        plt.figure(1).add_axes((0,0.2,0.6,0.5))
        plt.scatter(x_equal_spacing, y_data_fit, label='Observations')
        plt.plot(x_equal_spacing, result.best_fit, 'r-', label='Best Fit')
        plt.legend()
        plt.ylabel('Flux')

        plt.figure(1).add_axes((0,0,0.6,0.2)) #residual plot
        plt.plot(x_equal_spacing, r, linestyle = ':')
        plt.ylabel('Residual')
        plt.xlabel('Wavenumber / cm$^{-1}$')
        plt.show()
    plot_best_fit(result, x_equal_spacing, y_data_fit)
    
    return result


# if __name__ == "__main__":
#     Jmax = 300
#     sightline = '185418'

#     result = fit_model(B = 0.002, T = 22.5, delta_B = -0.45, zeta = -0.01, sigma = 0.17, origin =  0.012, sightline = sightline)
#     plt.figure(figsize = (15,8))


#%% Chi Squared function

def chi_sq(model_ys, obs_ys, std):
    '''
    Function to calcualate redued chi-squared value of model compared to observations. 
    
    model_ys and obs_ys must be of the same dimensions, this is checked and will raise an AssertionError if not.

    Parameters
    ----------
    model_ys : array-like
        Generated model fluxes, after having been passed through the np.interp process. 
    obs_ys : array-like
        Observed data.
    std : float
        Standard deviation of observations.

    Returns
    -------
    red_chi_sq : float
        Reduced chi-square value.

    '''
    assert len(model_ys) == len(obs_ys), 'Ensure the model and observed data are of the same dimension!'
    
    chi_sq = 0
    for i in range(len(obs_ys)):
        difference = (model_ys[i] - obs_ys[i])**2
        chi_sq += difference
    
    dof  = len(obs_ys) - 6
    chi_sq = chi_sq/(std**2)
    red_chi_sq = chi_sq/dof
    return red_chi_sq



#%% 
# print('\a')

