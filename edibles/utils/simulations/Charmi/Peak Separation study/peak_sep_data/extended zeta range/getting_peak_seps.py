# -*- coding: utf-8 -*-
"""
Created on Sat Nov 19 18:35:09 2022

@author: Charmi Bhatt
"""


import numpy as np
import pandas as pd
import astropy.constants as const
import matplotlib.pyplot as plt
import timeit
import scipy.stats as ss
from edibles.utils.edibles_oracle import EdiblesOracle
from edibles.utils.edibles_spectrum import EdiblesSpectrum
import warnings
from astropy.modeling import models
from astropy import units as u
from specutils.spectra import Spectrum1D
from specutils.fitting import fit_generic_continuum
import seaborn as sns
from scipy.signal import argrelextrema



origin = 0 
#Jmax = 300 (Kmax = Jmax (i.e all K allowed))
combinations = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Jmax=300.txt", delim_whitespace=(True))

pr_peak_seps =[]
pq_peak_seps =[]
qr_peak_seps =[]
qr_pq_diff = []
Temp = []
newold = []


startl = timeit.default_timer()

#%%
def get_rotational_spectrum(T, ground_B, delta_B, delta_C, zeta, sigma):
    
    ground_C = ground_B/2
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
    
    normalized_intensities = (intensities / max(intensities))
    linelist['normalized_intensities'] = normalized_intensities
    
   
    
    endl = timeit.default_timer()
    print('>>>> linelist calculation takes   ' + str(endl-startl) + '  sec')
    #%%
   
    '''Smoothening the linelist'''
    
    #%%
   
    
    #given that Resolution = 100,000 at wavelength (lambda) = 6614A, 
    #delta_lambda = 0.06614 and wavelength_stepsize = 0.033 (2 peaks per FWHM)
    #similarly for waveno (i.e 15120), delta_nu = 0.15 and thus waveno_stepsize = 0.075  
    # waveno_stepsize = 0.075
    # grid_size = int(((np.max(linelist['wavenos']) + 0.5) - (np.min(linelist['wavenos']) - 0.5))/waveno_stepsize)  

    smooth_wavenos = np.linspace(np.min(linelist['wavenos']) - 1 ,np.max(linelist['wavenos']) + 1, 1000) # grid_size)
    smooth_intensities = np.zeros(smooth_wavenos.shape)
    
    startg = timeit.default_timer()
    
    for idx,wavepoint in np.ndenumerate(smooth_wavenos):
        w_int = ss.norm.pdf(linelist['wavenos'], wavepoint, sigma) * (linelist['intensities']) 
        smooth_intensities[idx] = w_int.sum()
    
    # def gaussian(x, mu, sig):
    #         return (np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.))))/np.sqrt(2*np.pi*np.power(sig, 2.))

    # smooth_gau = np.zeros(smooth_wavenos.shape)
    # for i in linelist.index:
    #     smooth_gau = smooth_gau + (linelist['normalized_intensities'][i])*gaussian(smooth_wavenos, linelist['wavenos'][i], sigma)
        
    endg = timeit.default_timer()
    
    print('>>>> gaussian takes   ' + str(endg -startg) + '  sec') 
    
    smooth_data = np.array([smooth_wavenos, smooth_intensities]).transpose()    
    smooth_data = np.delete(smooth_data, np.where(smooth_data[:,1] <= 0.0001*(max(smooth_data[:,1]))), axis = 0)
    
    
    smooth_data = np.array([smooth_wavenos, smooth_intensities]).transpose()
    #print(smooth_data.shape)
    smooth_data = np.delete(smooth_data, np.where(smooth_data[:,1] <= 0.0001*(max(smooth_data[:,1]))), axis = 0)
    #print(smooth_data.shape)
    
    np.savetxt('smooth_1.txt',  smooth_data, delimiter= " ")
    
    x = smooth_data[:,0]
    y = np.array(1-0.1*(smooth_data[:,1]/max(smooth_data[:,1])))

    
    minima = [argrelextrema(y, np.less)]
    minima_ind = minima[0][0]
    minima_ind_plus_one = np.array(minima_ind) + 1
    minima_ind_minus_one = np.array(minima_ind) - 1
    
    P_parabola_x = (x[minima_ind_minus_one[0]], x[minima_ind[0]], x[minima_ind_plus_one[0]])
    Q_parabola_x = (x[minima_ind_minus_one[1]], x[minima_ind[1]], x[minima_ind_plus_one[1]])
    R_parabola_x = (x[minima_ind_minus_one[2]], x[minima_ind[2]], x[minima_ind_plus_one[2]])
    
    
    P_parabola_y = (y[minima_ind_minus_one[0]], y[minima_ind[0]], y[minima_ind_plus_one[0]])
    Q_parabola_y = (y[minima_ind_minus_one[1]], y[minima_ind[1]], y[minima_ind_plus_one[1]])
    R_parabola_y = (y[minima_ind_minus_one[2]], y[minima_ind[2]], y[minima_ind_plus_one[2]])
    
    
       
    
    
    vertex = []    
    def vertex_of_parabola(x,y):
        x_1 = x[0]
        x_2 = x[1]
        x_3 = x[2]
        y_1 = y[0]
        y_2 = y[1]
        y_3 = y[2]
    
        a = y_1/((x_1-x_2)*(x_1-x_3)) + y_2/((x_2-x_1)*(x_2-x_3)) + y_3/((x_3-x_1)*(x_3-x_2))
    
        b = (-y_1*(x_2+x_3)/((x_1-x_2)*(x_1-x_3))
              -y_2*(x_1+x_3)/((x_2-x_1)*(x_2-x_3))
              -y_3*(x_1+x_2)/((x_3-x_1)*(x_3-x_2)))
    
        c = (y_1*x_2*x_3/((x_1-x_2)*(x_1-x_3))
            +y_2*x_1*x_3/((x_2-x_1)*(x_2-x_3))
            +y_3*x_1*x_2/((x_3-x_1)*(x_3-x_2)))
        
        vertex.append(-b/(2*a))
    
    
    
    vertex_of_parabola(P_parabola_x, P_parabola_y)
    vertex_of_parabola(Q_parabola_x, Q_parabola_y)
    vertex_of_parabola(R_parabola_x, R_parabola_y)
    
    
    peak_sep_pq = (vertex[1] - vertex[0]) 
    peak_sep_qr = (vertex[2] - vertex[1])
    peak_sep_pr = (vertex[2] - vertex[0])
    
    # print(peak_sep_pq)
    # print(peak_sep_qr)
    # print(peak_sep_pr)

    pr_peak_seps.append(peak_sep_pr)
    pq_peak_seps.append(peak_sep_pq)
    qr_peak_seps.append(peak_sep_qr)
    qr_pq_diff.append(peak_sep_qr - peak_sep_pq)
    
    # print(pr_peak_seps)
    # print(pq_peak_seps)
    # print(qr_peak_seps)
    
    new = peak_sep_pr
    old = (x[minima_ind][-1] - x[minima_ind][0])
        
    newold.append(new - old)    
    
    Temp.append(T)
    
    
    #%%
    
    
    
    # with sns.color_palette("flare", n_colors=2):        
        
    #     'Calculated'
    plt.plot(x,y) #, label = str(peak_sep_pr))
        #plt.scatter(x[minima_ind], y[minima_ind], marker = '*')
        





Ts = (7, 10) #, 20, 30, 40, 50, 60, 70)
#T = 70
ground_B = 0.05
delta_B = 0
delta_C = 0
zeta = 0
sigma = 0.1953
conditions = 'condition c' 

#get_rotational_spectrum(T, ground_B, delta_B, delta_C, zeta, sigma)


for T in Ts:
    get_rotational_spectrum(T, ground_B, delta_B, delta_C, zeta, sigma)


peak_sep_data = np.column_stack([Temp, pq_peak_seps, qr_peak_seps, pr_peak_seps, newold, qr_pq_diff])
np.savetxt('peak_seps_zeta={}_B={}_delta_B ={}.txt'.format(zeta, ground_B, delta_B), peak_sep_data)
    
print(qr_pq_diff)

# #kerr's conditions   
   
# # Ts = (8.9, 20.2, 61.2, 101.3)    
# # ground_Bs = (0.01913, 0.00947, 0.00336, 0.00286)
# # delta_Bs = (-0.85, -0.42, -0.17, -0.21)
# # zeta = (-0.46, -0.43, -0.49, -0.54)
# # sigma = (0.1358, 0.1571, 0.1953, 0.1995)
# # conditions = ('condition a', 'condition b', 'condition c', 'condition d')

# # for T, B, d, z, sig, con in zip(Ts, ground_Bs, delta_Bs, zeta, sigma, conditions):
# #     get_rotational_spectrum(T, B, d, z, sig, con) 
        





'''Previously used codes'''

#%%

#pgopher_ll = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Kerr's conditions\condition_d\kerr's_condition_d_pgopher_linelist.txt", delim_whitespace=True)



# pgopher_ll['Intensity'] = pd.to_numeric(pgopher_ll['Intensity'], errors = 'coerce')
# print(pgopher_ll['Intensity'])



# pgopher = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Kerr's conditions\condition_d\dddd.txt", delim_whitespace=(True))

# pgopher_position = pgopher['position']
# pgopher_strength = 1 - 0.1*(pgopher['strength']/max(pgopher['strength']))

# print(pgopher['position'])
# print(pgopher['strength'])
# plt.figure(figsize=(30,6))
# plt.stem(pgopher_position, pgopher_strength,  label = 'pgopher', bottom = 1)
# plt.title('Pgopher Kerr condition d')
# plt.xlim(15118, 15122)
# plt.legend()

# def area_under_curve(x,y):
          
    #        sum_of_areas = 0
    #        for i in range(1, len(x)):
    #            h = smooth_wavenos[i] - smooth_wavenos[i-1]
    #            sum_of_areas += h * (smooth_intensities[i-1] + smooth_intensities[i]) / 2
        
    #        return sum_of_areas 

    # print('area under gaussian curve is  ' + str(area_under_curve(smooth_wavenos, smooth_intensities)))
    # print('area under scipy curve is  ' + str(area_under_curve(smooth_wavenos, scipy_smooth)))
    # print('sum of normalized intenisities is  ' + str(np.sum(normalized_intensities)))
    # print('---------------')
    
# for i in linelist.index:
    #     smooth_intensities = smooth_intensities + normalized_intensities[i]*gaussian(smooth_wavenos, wavenos[i], wavenos[i]/(2.355*resolution))  
    # endg = timeit.default_timer()
    # print('>>>> gaussian takes   ' + str(endg -startg) + '  sec')    
    # smooth_norm_intensities = 1 - 0.1*(smooth_intensities/max(smooth_intensities))
    
    # startss = timeit.default_timer()    
    # for i in linelist.index:
    #     scipy_smooth = scipy_smooth + normalized_intensities[i]*ss.norm.pdf(smooth_wavenos, wavenos[i], wavenos[i]/(2.355*resolution))    
    # endss = timeit.default_timer()
    # print('>>>> scipy takes   ' + str(endss -startss) + '  sec')
    # print('-------------')
    # scipy_smooth_norm = 1 - 0.1*(scipy_smooth/max(scipy_smooth))
    


# wavelength = []
#     for i in range(len(linelist['wavenos'])):
#         wavelength.append(1/linelist['wavenos'].iloc[i]*1e8)
     
#     wavelength_spacing = 0.033
#     grid_size = int(((np.max(wavelength) + 0.05) - (np.min(wavelength) - 0.05))/wavelength_spacing)  

 # print('max wavelength is:' + str(np.max(wavelength)))
 #    print('min wavelength is:' + str(np.min(wavelength)))
 #    print('Max - min wavelength is:' + str(np.max(wavelength)-np.min(wavelength)))
 #    print('grid size is: ' + str(grid_size))
    
    
 #    print('delta lambda is :')
 #    for i in range(len(smooth_wavelength[0:3])):
 #        print(smooth_wavelength[i+1] - smooth_wavelength[i])


 # plt.figure(figsize=(30,6))  
    
    
    # smooth_norm_intensities = (smooth_intensities/max(smooth_intensities))
    
#linelist.to_excel(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\work.xlsx", index=False)
      
 #linelist = linelist[(linelist['label'] == 'rR')]
 #linelist = linelist[(linelist['ground_K'] <= 5)]
 
 # for i in range(len(smooth_wavenos[0:3])):
 #         print(smooth_wavenos[i+1] - smooth_wavenos[i])
    
    
#kerr_1996

# Ts = (8.9, 20.2, 61.2, 101.3)    
# ground_Bs = (0.01913, 0.00947, 0.00336, 0.00286)
# delta_Bs = (-0.85, -0.42, -0.17, -0.21)

# for T, B, d in zip(Ts, ground_Bs, delta_Bs):
#     get_rotational_spectrum(T, B, d)

#startplot =  timeit.default_timer()
# Ts = np.linspace(10, 100, 1)
# #ground_Bs = (0.001, 0.003, 0.007, 0.01, 0.03, 0.07, 0.1)
# ground_B = 0.1
# delta_Bs = (-1, -2, -3)


# for T in Ts:
#     plt.figure(figsize=(22,6))
#     for delta_B in delta_Bs:
#         get_rotational_spectrum(T, ground_B , delta_B)
#     plt.show()

# endplot =  timeit.default_timer()

# print(endplot - startplot)

# T = 8.9
# ground_B = 0.01913
# delta_B = -0.85
# zeta = -0.46
# sigma = 0.1358
# conditions = 'condition a'

# T = 20.2
# ground_B = 0.00947
# delta_B = -0.42
# zeta = -0.43
# sigma = 0.1571
# conditions = 'condition b'

# T = 61.2
# ground_B = 0.00336
# delta_B = -0.17
# zeta = -0.49
# sigma = 0.1953
# conditions = 'condition c'

# T = 101.3
# ground_B = 0.00286
# delta_B = -0.21
# zeta = -0.54
# sigma = 0.1995
# conditions = 'condition d'

# P_Branch = linelist[(linelist['label'].str[1] == "P")]
# Q_Branch = linelist[(linelist['label'].str[1] == "Q")]
# R_Branch = linelist[(linelist['label'].str[1] == "R")]



# peak_p = linelist[linelist['intensities'] == max(P_Branch['intensities'])]
# peak_q = linelist[linelist['intensities'] == max(Q_Branch['intensities'])]
# peak_r = linelist[linelist['intensities'] == max(R_Branch['intensities'])]

# #linelist.to_excel(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Calculated_linelist_kerr_condition_c.xlsx", index=False)
