# -*- coding: utf-8 -*-
import sys
import numpy as np
import pandas as pd
import astropy.constants as const
import matplotlib.pyplot as plt
from edibles.utils.edibles_oracle import EdiblesOracle
from edibles.utils.edibles_spectrum import EdiblesSpectrum
import warnings
from astropy.modeling import models
from astropy import units as u
from specutils.spectra import Spectrum1D
from specutils.fitting import fit_generic_continuum
import timeit


plt.figure(figsize=(50,6))




def get_rotational_spectrum(T, ground_B, delta_B):
    
    ground_C = ground_B/2
    delta_C = delta_B
    
    origin = 15120
    Jmax = 300 #Kmax = Jmax (i.e all K allowed)
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
    
    print('Jmax is  ' + str(Jmax))
    print('length of linelist  ' + str(len(linelist)))
    endc = timeit.default_timer()
    print('---------------')
    print('>>>> combnination calclulation takes   ' + str(endc-startc) + '  sec')
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
    
    normalized_intensities =  1 - 0.1*(intensities / max(intensities))
    linelist['normalized_intensities'] = normalized_intensities
    
   
    
    endl = timeit.default_timer()
    print('>>>> linelist calculation takes   ' + str(endl-startl) + '  sec')
    startg = timeit.default_timer()
    #%%
    def gaussian(x, mu, sig):
                return (np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.))))/np.sqrt(2*np.pi*np.power(sig, 2.))
    
    
            
    smooth_wavenos = np.linspace(np.min(wavenos) - 0.5 ,np.max(wavenos) + 0.5, 10000)
    smooth_intensities = np.zeros(smooth_wavenos.shape)
    
    
    
    
    # for i in linelist.index:
    #     smooth_intensities = smooth_intensities + normalized_intensities[i]*gaussian(smooth_wavenos, wavenos[i], wavenos[i]/(2.355*resolution))
    
    
    # smooth_norm_intensities = 1 - 0.1*(smooth_intensities/max(smooth_intensities))
    
    # endg = timeit.default_timer()
    # print('>>>> gaussian takes   ' + str(endg -startg) + '  sec')
    # print('-------------')
    wavelength = []
    for i in range(len(wavenos)):
        wavelength.append(1/wavenos[i]*1e8)
        
    smooth_wavelength = 1/smooth_wavenos*1e8
    
#%%
    starName = 'HD 166937'
    #put lower range of wavelengths to extract from edibles data
    minWave = 6612
    
    #put upper range of wavelengths to extract from edibles data
    maxWave = 6616
    
    pythia = EdiblesOracle()
    rawList = pythia.getFilteredObsList(object = [starName], MergedOnly = True, WaveMin = minWave, WaveMax = maxWave)
    fnames = rawList.tolist()
    obs = len(fnames)
    
    
    sp = EdiblesSpectrum(fnames[0])
        
    sp.getSpectrum(xmin = max(minWave, np.min(sp.raw_wave)+1)
                        , xmax = min(maxWave, np.max(sp.raw_wave)-1))
    
                       
    #data = np.array([sp.bary_wave, sp.bary_flux]).transpose()
    leftEdge = 0
    rightEdge = 0
        
    if minWave <= np.min(sp.raw_wave):
        leftEdge = 1
        #print('Left edge detected')
    if maxWave >= np.max(sp.raw_wave):
        rightEdge = 1
    
    data = np.delete(np.array([sp.bary_wave, sp.bary_flux]).transpose(), 
                                np.logical_or(sp.bary_wave <= np.min(sp.bary_wave) + 40.0*leftEdge, 
                                              sp.bary_wave >= np.max(sp.bary_wave) - 40.0*rightEdge), 0)
    
    v = -6.5
    
    data[:, 0] = data[:, 0]*(1+v/299792.458)
    
    x1 = data[:,0]
    y1= data[:,1]
    
    spectrum1 = Spectrum1D(flux = y1*u.dimensionless_unscaled, spectral_axis = x1*u.angstrom)
    
    with warnings.catch_warnings():  # Ignore warnings
        warnings.simplefilter('ignore')
        g1_fit = fit_generic_continuum(spectrum1, model = models.Legendre1D(degree = 5))
    
    
    data[:,1] = y1/g1_fit(x1*u.angstrom)
    
    plt.figure(figsize=(20,6))
    #plt.plot(data[:, 0] + 0.5, data[:, 1]/max(data[:, 1]))
    
    
    #%%
    
    plt.figure(figsize=(30,6))
    plt.stem(wavenos, normalized_intensities,  label='calculated', bottom = 1, linefmt='y', markerfmt='yo')
    #plt.plot(smooth_wavelength, (smooth_norm_intensities))
    plt.title('Calculated: T = ' + str(T) + 'K ,  ground_B =  ' + str(ground_B))
    plt.xlim(15118, 15122)
    plt.show()
    
    
    


 
get_rotational_spectrum(8.9, 0.01913, -0.85)

pgopher = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Kerr's conditions\condition a\pgopher_kerr_condition_a_Jmax_300_A1g_E1u.txt", delim_whitespace=(True))

pgopher_position = pgopher['position']
pgopher_strength = 1 - 0.1*(pgopher['strength']/max(pgopher['strength']))

plt.figure(figsize=(30,6))
plt.stem(pgopher_position, pgopher_strength,  label = 'pgopher', bottom = 1)
plt.title('Pgopher Kerr condition a')
plt.xlim(15118, 15122)
plt.legend()










       
    
#kerr_1996

# Ts = (8.9, 20.2, 61.2, 101.3)    
# ground_Bs = (0.01913, 0.00947, 0.00336, 0.00286)
# delta_Bs = (-0.85, -0.42, -0.17, -0.21)

# for T, B, d in zip(Ts, ground_Bs, delta_Bs):
#     get_rotational_spectrum(T, B, d)