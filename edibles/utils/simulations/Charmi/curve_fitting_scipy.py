import numpy as np
import pandas as pd
import astropy.constants as const
import matplotlib.pyplot as plt
import timeit
import scipy
import scipy.stats as ss
from lmfit import Model
# from edibles.utils.edibles_oracle import EdiblesOracle
# from edibles.utils.edibles_spectrum import EdiblesSpectrum
# import warnings
# from astropy.modeling import models
# from astropy import units as u
# from specutils.spectra import Spectrum1D
# from specutils.fitting import fit_generic_continuum
# import seaborn as sns
# from scipy.signal import argrelextrema
# from lmfit import Model


full_start = timeit.default_timer()

Obs_data = pd.read_csv(r"/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/Heather's_data/6614_HD_166937.txt", sep = ',')

y_obs_data =  np.array(Obs_data['Flux'])
x_obs_data = np.array(Obs_data['Wavelength'])

# x_for_model = np.linspace(min(x_obs_data), max(x_obs_data), len(x_obs_data))
# y_obs_data = np.interp(x_for_model, x_obs_data, y_obs_data) 

plt.Figure(figsize=(15,8))

plt.plot(x_obs_data, y_obs_data)


min_index = np.argmin(y_obs_data)

central_peak_wavelength = x_obs_data[min_index]

origin = (1/central_peak_wavelength)*1e8
#origin = 15120.9
print(origin)
#origin = 15120.9
#Jmax = 300 (Kmax = Jmax (i.e all K allowed))
combinations = pd.read_csv(r"/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/Jmax=300.txt", delim_whitespace=(True))




#%%
def get_rotational_spectrum(xx, B, T):
    
    startl = timeit.default_timer()
    print(B)
    print(T)
    x_obs_data = xx
    ground_B = B
    delta_B = -0.17
    delta_C = -0.17
    zeta = -0.49
    sigma = 0.1953
    
    #rotational constants in cm-1
    ground_C = B/2
    delta_C = delta_C
    excited_B = ground_B + ((delta_B/100)*ground_B)
    excited_C = ground_C + ((delta_C/100)*ground_C)
    
    global combinations
    ground_Js = combinations['ground_J']
    excited_Js = combinations['excited_J']
    ground_Ks = combinations['ground_K']
    excited_Ks = combinations['excited_K']
    
    linelist = combinations
   
    delta_J = linelist['excited_J'] - linelist['ground_J']
    delta_K = linelist['excited_K'] - linelist['ground_K']
    
    '''Calculating Linelist'''
    #%%
    
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
    
    # normalized_intensities = (intensities / max(intensities))
    # linelist['normalized_intensities'] = normalized_intensities
    
   
    
    endl = timeit.default_timer()
    print('>>>> linelist calculation takes   ' + str(endl-startl) + '  sec')
    #%%
   
    '''Smoothening the linelist'''
    
    #%%
   
    smooth_wavenos = np.linspace(np.min(linelist['wavenos']) - 1 ,np.max(linelist['wavenos']) + 1, 1000) # grid_size)
    smooth_intensities = np.zeros(smooth_wavenos.shape)
    
    startg = timeit.default_timer()
    
    for idx,wavepoint in np.ndenumerate(smooth_wavenos):
        w_int = ss.norm.pdf(linelist['wavenos'], wavepoint, sigma) * (linelist['intensities']) 
        smooth_intensities[idx] = w_int.sum()
    

    endg = timeit.default_timer()
    
    print('>>>> gaussian takes   ' + str(endg -startg) + '  sec') 
    
    smooth_data = np.array([smooth_wavenos, smooth_intensities]).transpose()    
    smooth_data = np.delete(smooth_data, np.where(smooth_data[:,1] <= 0.001*(max(smooth_data[:,1]))), axis = 0)
    
    simu_waveno = smooth_data[:, 0]
    simu_wavelength = (1/simu_waveno)*1e8
    simu_intenisty = 1-0.1*(smooth_data[:, 1]/max(smooth_data[:, 1]))
    
    
    #%%
    
    P_Branch = linelist[(linelist['label'].str[1] == "P")]
    Q_Branch = linelist[(linelist['label'].str[1] == "Q")]
    R_Branch = linelist[(linelist['label'].str[1] == "R")]

    
    
    peak_p = linelist[linelist['intensities'] == max(P_Branch['intensities'])]
    peak_q = linelist[linelist['intensities'] == max(Q_Branch['intensities'])]
    peak_r = linelist[linelist['intensities'] == max(R_Branch['intensities'])]

    #linelist.to_excel(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Calculated_linelist_kerr_condition_c.xlsx", index=False)

    #with sns.color_palette("flare", n_colors=2):
        
    plt.plot(simu_wavelength, simu_intenisty, color = 'red')
    
    # plt.xlabel('Wavelength')
    # plt.ylabel('Normalized Intenisty')
    # plt.title('Temperature = ' + str(T) + '  K  ground_B =  ' + str(ground_B) + ' cm-1  ground_C=  ' + str(ground_C) + ' cm-1  Delta_B = ' + str(delta_B) + '    Delta_C = ' + str(delta_C) +    '    zeta = ' +  str(zeta)) 
    
    model_data = np.array([simu_wavelength, simu_intenisty]).transpose()
    model_data = model_data[::-1]
    
    x_for_model = np.linspace(min(x_obs_data), max(x_obs_data), len(x_obs_data))
    
    y_model_data = np.interp(x_for_model, model_data[:,0], model_data[:,1])
    
    
    
    return  y_model_data
    
# Bs = np.linspace(0.0030, 0.0045, 10)  
# Ts = np.linspace(45, 65, 20)  

# # get_rotational_spectrum(x_obs_data, B = 0.00411, T = 48.16)
# # get_rotational_spectrum(x_obs_data, B = 0.00305, T = 62.5)

# num = (get_rotational_spectrum(x_obs_data, B = 0.00389999, T = 52) - y_obs_data)**2
# chi_squared = np.sum((num)/(0.004)**2)
# reduced_chi_squared = chi_squared/(len(y_obs_data) - 2)
# #reduced_chis.append(reduced_chi_squared)
# print(reduced_chi_squared)

# for B in Bs:
#     for T in Ts:   
#         num = (get_rotational_spectrum(x_obs_data, B, T) - y_obs_data)**2
#         chi_squared = np.sum((num)/(0.004)**2)
#         reduced_chi_squared = chi_squared/(len(y_obs_data) - 2)
#         #reduced_chis.append(reduced_chi_squared)
#         print(reduced_chi_squared)
        
# chi_cor = np.array([B, T, reduced_chi_squared]).transpose()
# print(chi_cor)


mod = Model(get_rotational_spectrum) #, independent_vars = ['b', 'T']) #make sure independent variable of fitting function (that you made) is labelled as x
#params = mod.guess(flux_data, x = np.linspace(0.005,0.01,5))
params = mod.make_params(verbose = True, B = 0.0029, T = 59)


# params['b'].min = 0.005 
# params['b'].max = 0.01
# params['T'].min = 2.7
# params['T'].max = 100

result = mod.fit(y_obs_data, params, xx= x_obs_data, weights = 1/0.004) #, b = 0.005, T = 3, weights = 1/0.7)

#plt.plot(x_obs_data, y_obs_data, label = 'Data')
plt.plot(x_obs_data, result.best_fit, label = 'Fit')
plt.legend()

print(mod.param_names, mod.independent_vars)
print(result.fit_report())


'''scipy'''
# guess = [0.003, 19]
#popt, pcov = scipy.optimize.curve_fit(get_rotational_spectrum, xdata = Obs_data['Wavelength'], ydata = Obs_data['Flux'], p0 = guess) #, bounds = ([0.003, 19], [0.0033, 33]))
# # #popt, pcov = scipy.optimize.least_squares(get_rotational_spectrum, xdata = Obs_data['Wavelength'], ydata = Obs_data['Flux'], p0 = guess, bounds = ([2.7, 0.005], [5, 0.007]), xtol = 0.05, ftol = 0.05)

# print(popt)
# print(pcov)

# xx = x_obs_data
# plt.plot(x_obs_data, get_rotational_spectrum(xx, *popt), 'r-',
#           label='fit: B=%5.3f, T=%5.3f' % tuple(popt))
# plt.legend()

# full_end = timeit.default_timer()
# print(full_end - full_start)