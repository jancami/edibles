# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import astropy.constants as const
import matplotlib.pyplot as plt
import timeit
import scipy.stats as ss
from astropy.convolution import Gaussian1DKernel, convolve


plt.figure(figsize=(50,6))




def get_rotational_spectrum(T, ground_B, delta_B):
    
    ground_C = ground_B/2
    delta_C = delta_B
    
    origin = 15120
    Jmax = 400 #Kmax = Jmax (i.e all K allowed)
    resolution = 100000
    
    startc = timeit.default_timer()
    
    '''Calculating Linelist'''
    #%%
    excited_B = ground_B + ((delta_B/100)*ground_B)
    excited_C = ground_C + ((delta_C/100)*ground_C)
    P_branch_Js = list((range(1,Jmax+1)))
    all_P_branch_Js = []
    for j in P_branch_Js:
        for i in range(j):
          all_P_branch_Js.append(j)
        
    P_branch_Jprimes = []
    for j in all_P_branch_Js:
        if j != 0:
            P_branch_Jprimes.append(j-1)
            
    pP_Branch_K = []        
    for j in P_branch_Js:
        stages = list(range(0,j))
        for i in stages:
            pP_Branch_K.append(j-i)
            
    
    pP_Branch_Kprime = []
    for k in pP_Branch_K:
        pP_Branch_Kprime.append(k-1)
            
    rP_Branch_K = []        
    for j in P_branch_Js:
        stages = list(range(0,j))
        stages.sort(reverse=True)
        for i in stages:
            rP_Branch_K.append(i)
    
    
    rP_Branch_Kprime = []
    for k in rP_Branch_K:
        rP_Branch_Kprime.append(k+1)
        
    
    '''Q Branch'''
    
    Q_branch_Js = list((range(0,Jmax+1)))
    
    all_Q_branch_Js = []
    for j in Q_branch_Js:
        # if j ==0:
        #     all_Q_branch_Js.append(j)
        if j!= 0:
            for i in range(j):
              all_Q_branch_Js.append(j)
          
    Q_branch_Jprimes = []
    for j in all_Q_branch_Js:
            Q_branch_Jprimes.append(j)
            
    pQ_Branch_K = []        
    for j in Q_branch_Js:
        stages = list(range(0,j))
        for i in stages:
            pQ_Branch_K.append(j-i)
            
    
    pQ_Branch_Kprime = []
    for k in pQ_Branch_K:
        pQ_Branch_Kprime.append(k-1)
        
    rQ_Branch_K = []        
    for j in Q_branch_Js:
        stages = list(range(0,j))
        stages.sort(reverse=True)
        for i in stages:
            rQ_Branch_K.append(i)
    
    
    rQ_Branch_Kprime = []
    for k in rQ_Branch_K:
        rQ_Branch_Kprime.append(k+1)
        
            
    
            
    '''R Branch'''
            
    R_branch_Js = list((range(0,Jmax)))
    all_R_branch_Js = []
    for j in R_branch_Js:
        if j ==0:
            all_R_branch_Js.append(j)
        elif j!= 0:
            for i in range(j+1):
              all_R_branch_Js.append(j)
                  
    R_branch_Jprimes = []
    for j in all_R_branch_Js:
        if j <= Jmax-1:
            R_branch_Jprimes.append(j+1)
            
    pR_Branch_K = []        
    for j in R_branch_Js:
        stages = list(range(0,j+1))
        # if j!= 0:
        for i in stages:
            pR_Branch_K.append(j-(i-1))
        
    
    pR_Branch_Kprime = []
    for k in pR_Branch_K:
        pR_Branch_Kprime.append(k-1)
        
    rR_Branch_K = []        
    for j in R_branch_Js:
        stages = list(range(0,j+1))
        stages.sort(reverse=True)
        for i in stages:
            rR_Branch_K.append(i)
    
    
    rR_Branch_Kprime = []
    for k in rR_Branch_K:
        rR_Branch_Kprime.append(k+1)
        
            
    
    
    
    Allowed_Js = (all_P_branch_Js*2) + (all_Q_branch_Js*2) + (all_R_branch_Js*2)
    Allowed_Jprimes = (P_branch_Jprimes*2) + (Q_branch_Jprimes*2) + (R_branch_Jprimes*2)
    Allowed_Ks = pP_Branch_K + rP_Branch_K + pQ_Branch_K + rQ_Branch_K +  pR_Branch_K + rR_Branch_K
    Allowed_Kprimes = pP_Branch_Kprime + rP_Branch_Kprime + pQ_Branch_Kprime + rQ_Branch_Kprime + pR_Branch_Kprime + rR_Branch_Kprime
    
    columns = {'ground_J' : Allowed_Js,'excited_J': Allowed_Jprimes, 'ground_K' : Allowed_Ks, 'excited_K' : Allowed_Kprimes}
    linelist = pd.DataFrame(columns)
    
    linelist['delta_J'] = linelist['excited_J'] - linelist['ground_J']
    linelist['delta_K'] = linelist['excited_K'] - linelist['ground_K']
    
    label = []
    
    for i in range(len(linelist['ground_J'])):
        if linelist['excited_J'][i] - linelist['ground_J'][i] == -1 and linelist['excited_K'][i] - linelist['ground_K'][i] == -1:
            label.append('pP')
        if linelist['excited_J'][i] - linelist['ground_J'][i] == -1 and linelist['excited_K'][i] - linelist['ground_K'][i] == 1:
            label.append('rP')
        if linelist['excited_J'][i] - linelist['ground_J'][i] == 0 and linelist['excited_K'][i] - linelist['ground_K'][i] == -1:
            label.append('pQ')
        if linelist['excited_J'][i] - linelist['ground_J'][i] == 0 and linelist['excited_K'][i] - linelist['ground_K'][i] == 1:
            label.append('rQ')
        if linelist['excited_J'][i] - linelist['ground_J'][i] == 1 and linelist['excited_K'][i] - linelist['ground_K'][i] == -1:
            label.append('pR')
        if linelist['excited_J'][i] - linelist['ground_J'][i] == 1 and linelist['excited_K'][i] - linelist['ground_K'][i] == 1:
            label.append('rR')
    
    linelist['Label'] = label
    
    ground_Js = linelist['ground_J']
    excited_Js = linelist['excited_J']
    ground_Ks = linelist['ground_K']
    excited_Ks = linelist['excited_K']
    delta_J = linelist['delta_J']
    delta_K = linelist ['delta_K']
    
    
    print('----------------------')
    print('Jmax is  ' + str(Jmax))
    print('length of linelist  ' + str(len(linelist)))
    endc = timeit.default_timer()
    print('---------------')
    print('>>>> combination calclulation takes   ' + str(endc-startc) + '  sec')
    startl = timeit.default_timer()
    
    ground_Es = []
    for J,K in zip(ground_Js, ground_Ks):
                ground_E = ground_B*J*(J + 1) + (ground_C - ground_B)*(K**2)
                ground_Es.append(ground_E)
                
    linelist['ground_Es'] = ground_Es 
    
    excited_Es = []
    for J,K in zip(excited_Js, excited_Ks):
                excited_E = excited_B*J*(J + 1) + (excited_C - excited_B)*(K**2)
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
    
    normalized_intensities =  intensities / max(intensities)
    linelist['normalized_intensities'] = normalized_intensities
    
   
    
    endl = timeit.default_timer()
    print('>>>> linelist calculation takes   ' + str(endl-startl) + '  sec')
    
    #%%
    
    

    
    def gaussian(x, mu, sig):
                return (np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.))))/np.sqrt(2*np.pi*np.power(sig, 2.))
    
    

    smooth_wavenos = np.linspace(np.min(wavenos) - 0.5 ,np.max(wavenos) + 0.5, 1000)
    smooth_intensities = np.zeros(smooth_wavenos.shape)
    scipy_smooth = np.zeros(smooth_wavenos.shape)
    astropy_smooth = np.zeros(smooth_wavenos.shape)
    
    
    startg = timeit.default_timer()    
    # Go over each point in the output grid, then multiply all line intensities by the proper Gaussian
    # centered on that data point, and sum it all up. That should be the intensity at that point. 
    for idx,wavepoint in np.ndenumerate(smooth_wavenos):
        w_int = ss.norm.pdf(wavenos,wavepoint,wavepoint/(2.355*resolution)) * normalized_intensities
        smooth_intensities[idx] = w_int.sum()
    #for i in linelist.index:
    #    smooth_intensities = smooth_intensities + normalized_intensities[i]*gaussian(smooth_wavenos, wavenos[i], wavenos[i]/(2.355*resolution))  
    endg = timeit.default_timer()
    print('>>>> gaussian takes   ' + str(endg -startg) + '  sec')    
    smooth_norm_intensities = (smooth_intensities/max(smooth_intensities))
    
    
    #startss = timeit.default_timer()    
    #for i in linelist.index:
    #    scipy_smooth = scipy_smooth + normalized_intensities[i]*ss.norm.pdf(smooth_wavenos, wavenos[i], wavenos[i]/(2.355*resolution))    
    #endss = timeit.default_timer()
    #print('>>>> scipy takes   ' + str(endss -startss) + '  sec')
    #print('-------------')
    #scipy_smooth_norm = 1 - 0.1*(scipy_smooth/max(scipy_smooth))
    
    
    def area_under_curve(x,y):
          
           sum_of_areas = 0
           for i in range(1, len(x)):
               h = smooth_wavenos[i] - smooth_wavenos[i-1]
               sum_of_areas += h * (smooth_intensities[i-1] + smooth_intensities[i]) / 2
        
           return sum_of_areas 

    print('area under gaussian curve is  ' + str(area_under_curve(smooth_wavenos, smooth_intensities)))
    #print('area under scipy curve is  ' + str(area_under_curve(smooth_wavenos, scipy_smooth)))
    print('sum of normalized intenisities is  ' + str(np.sum(normalized_intensities)))
    print('---------------')
    
    
    
    
    wavelength = []
    for i in range(len(wavenos)):
        wavelength.append(1/wavenos[i]*1e8)
        
    smooth_wavelength = 1/smooth_wavenos*1e8
    
#%%
    
    #%%
    
    plt.figure(figsize=(30,6))
    plt.stem(wavelength, normalized_intensities/max(normalized_intensities),  label='calculated', linefmt='y', markerfmt='yo')
    #plt.yscale('log')
    #plt.ylim(.01,1)
    #plt.plot(smooth_wavelength, scipy_smooth_norm, color='black')
    plt.plot(smooth_wavelength, smooth_intensities/max(smooth_intensities), color='blue')
    plt.title('Calculated: T = ' + str(T) + 'K ,  ground_B =  ' + str(ground_B))
    #plt.xlim(15118, 15122)
    plt.show()
    
    
    

T = 101.3
B_ground = 0.00286
delta_B = -0.21
 
get_rotational_spectrum(T, B_ground, delta_B)

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










       
    
#kerr_1996

# Ts = (8.9, 20.2, 61.2, 101.3)    
# ground_Bs = (0.01913, 0.00947, 0.00336, 0.00286)
# delta_Bs = (-0.85, -0.42, -0.17, -0.21)

# for T, B, d in zip(Ts, ground_Bs, delta_Bs):
#     get_rotational_spectrum(T, B, d)