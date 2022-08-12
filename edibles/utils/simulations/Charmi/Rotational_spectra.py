"""
>find wa   overplot for different molecules(objects)

"""

import pandas as pd
import matplotlib.pyplot as plt
import astropy.constants as const
from astropy.convolution import convolve, Gaussian1DKernel
import numpy as np
import timeit

class Rotational_Spectra:
    
    '''
    Give a summary here  
    '''
    
    
    def __init__(self, ground_B, ground_C, delta_B, delta_C, temperature, origin, Jmax=None, Kmax=None):
        
        
        self.ground_B = ground_B
        self.ground_C = ground_C
        self.delta_B = delta_B
        self.delta_C = delta_C
        self.Jmax = Jmax
        self.Kmax = Kmax
        self.T = temperature
        self.origin = origin
        
        self.excited_B = self.ground_B - ((self.delta_B/100)*self.ground_B)
        self.excited_C = self.ground_C - ((self.delta_C/100)*self.ground_C)
        
        self.selection_rules = pd.DataFrame()
        self.linelist = pd.DataFrame()
        # self.ground_J = []
        # self.ground_K = []
        self.ground_energies = []
        self.excited_energies = []
        self.wavenos = []
        self.HL_factors = []
        self.BD_factors = []
        self.intensities = []
        self.normalized_intensities = []
        self.BD_with_partition = []
        
        
    def possible_JandKs(self):
            Jmax = self.Jmax

            #here Kmax = Jmax (i.e all K allowed)
            
            
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
            selection_rules = pd.DataFrame(columns)
            
            selection_rules['delta_J'] = selection_rules['excited_J'] - selection_rules['ground_J']
            selection_rules['delta_K'] = selection_rules['excited_K'] - selection_rules['ground_K']
            
            label = []
            
            for i in range(len(selection_rules['ground_J'])):
                if selection_rules['excited_J'][i] - selection_rules['ground_J'][i] == -1 and selection_rules['excited_K'][i] - selection_rules['ground_K'][i] == -1:
                    label.append('pP')
                if selection_rules['excited_J'][i] - selection_rules['ground_J'][i] == -1 and selection_rules['excited_K'][i] - selection_rules['ground_K'][i] == 1:
                    label.append('rP')
                if selection_rules['excited_J'][i] - selection_rules['ground_J'][i] == 0 and selection_rules['excited_K'][i] - selection_rules['ground_K'][i] == -1:
                    label.append('pQ')
                if selection_rules['excited_J'][i] - selection_rules['ground_J'][i] == 0 and selection_rules['excited_K'][i] - selection_rules['ground_K'][i] == 1:
                    label.append('rQ')
                if selection_rules['excited_J'][i] - selection_rules['ground_J'][i] == 1 and selection_rules['excited_K'][i] - selection_rules['ground_K'][i] == -1:
                    label.append('pR')
                if selection_rules['excited_J'][i] - selection_rules['ground_J'][i] == 1 and selection_rules['excited_K'][i] - selection_rules['ground_K'][i] == 1:
                    label.append('rR')
            
            selection_rules['Label'] = label
            self.selection_rules = selection_rules
            
            return self.selection_rules
            
    def ground_Js(self):
        
       self.possible_JandKs()
       self.ground_J = self.selection_rules['ground_J']
       return self.ground_J
   
    def ground_Ks(self):
        
       self.possible_JandKs()
       self.ground_K = self.selection_rules['ground_K']
       return self.ground_K
   
    def excited_Js(self):
        
       self.possible_JandKs()
       self.ground_J = self.selection_rules['excited_J']
       return self.ground_J
   
    def excited_Ks(self):
         
       self.possible_JandKs()
       self.ground_J = self.selection_rules['excited_K']
       return self.ground_J
   
    def Branch_labels(self):
        self.possible_JandKs()
        self.label = self.selection_rules['Label']
        return self.label
         
     
    def ground_energy_levels(self):
        
        
        for J,K in zip(self.ground_Js(), self.ground_Ks()):
            ground_E = self.ground_B*J*(J + 1) + (self.ground_C - self.ground_B)*(K**2)
            self.ground_energies.append(ground_E)
            
        
        return self.ground_energies
            
    def excited_energy_levels(self):
        
        #excited_energies = []
        
        for J,K in zip(self.excited_Js(), self.excited_Ks()):
            excited_E = self.excited_B*J*(J + 1) + (self.excited_C - self.excited_B)*(K**2)
            self.excited_energies.append(excited_E)
            
        # BJ(J+1) + (C-B)K**2   
        # self.excited_energies = excited_energies
        # self.linelist['excited_energies'] = self.excited_energies
            
        return self.excited_energies
            
    def wavenumber(self):

        #wavenos = []

        for i in range(len(self.ground_energies)):
            waveno = self.origin + self.excited_energy_levels()[i] - self.ground_energy_levels()[i]
            self.wavenos.append(waveno)
            
        #self.wavenos = wavenos
        #self.linelist['wavenumber'] = self.wav       
        
        return self.wavenos
    

    def honl_london_factors(self):
        
        self.HL_factors = []
                
        for J,K, Jprime, Kprime in zip(self.ground_Js(), self.ground_Ks(), self.excited_Js(), self.excited_Ks()):
            if Jprime -J == -1 and Kprime - K == -1:
                HL_factor = ((J - 1 + K)*(J + K))/(J*((2*J)+1))
            elif Jprime -J  == -1 and Kprime - K == 1:
                HL_factor = ((J - 1 - K)*(J - K))/(J*((2*J)+1))
            elif Jprime -J  == 0 and Kprime - K == -1:
                HL_factor = (J+1-K)*(J+K)/ (J*(J + 1))
            elif Jprime -J == 0 and Kprime - K == 1:
                HL_factor = (J+1+K)*(J-K) / (J*(J + 1))
            elif Jprime -J  == 1 and Kprime - K == -1:
                HL_factor = ( J + 2 - K)*( J + 1 - K)/ ((J + 1)*((2*J) +1))
            elif Jprime -J  == 1 and Kprime - K == 1:
                HL_factor = ( J + 2 + K)*( J + 1 + K)/ ((J + 1)*((2*J) + 1))
                
            self.HL_factors.append(HL_factor)
            
        # self.HL_factors = HL_factors
        # self.linelist['Hl_factors'] = self.HL_factors
        
            
        return self.HL_factors


    def Boltzmann(self):
        
        self.BD_factors = []
        
        h = const.h.cgs.value
        c = const.c.to('cm/s').value
        k = const.k_B.cgs.value
        
        #for i in range(len(linelist.index)):
        for J,K,E in zip(self.ground_Js(), self.ground_Ks(), self.ground_energy_levels()):
            if K == 0:
                boltzmann_equation = (2*((2*J) + 1))*(np.exp((-h * c * E) / (k*self.T)))
            else:
                boltzmann_equation = (1*((2*J) + 1))*(np.exp((-h * c * E) / (k*self.T)))
                
            self.BD_factors.append(boltzmann_equation)
            
        # self.BD_factors = BD_factors
        # self.linelist['BD_factors'] = self.BD_factors
                
        return self.BD_factors
    
    
    def intensity(self):
        
        for i in range(len(self.ground_Js())):
            strength = (self.honl_london_factors()[i] * self.Boltzmann()[i])
            self.intensities.append(strength)
            
        
        return self.intensities 
    
    def normalized_intensity(self):     
        
        #self.intensity()
        self.normalized_intensities = (self.intensities / max(self.intensities))
        #self.normalized_intensities = 1 - self.normalized_intensities
        
        return self.normalized_intensities
        
    def get_linelist(self):
        
        self.linelist['ground_J'] = self.ground_Js()
        self.linelist['excited_J'] = self.excited_Js()
        self.linelist['ground_K'] = self.ground_Ks()
        self.linelist['excited_K'] = self.excited_Js()  
        self.linelist['label'] = self.Branch_labels()
        self.linelist['groundE'] = self.ground_energy_levels()
        self.linelist['excitedE'] = self.excited_energy_levels()
        self.linelist['wavenos'] = self.wavenumber()
        self.linelist['HL_factors'] = self.honl_london_factors()
        self.linelist['BD_factors'] = self.Boltzmann()
        self.linelist['intensities'] = self.intensity()
        self.linelist['normalized intensities'] = self.normalized_intensity()
        
        return self.linelist
    
    def plot_rotational_spectra(self):
        
        # self.wavenumber()
        # self.normalized_intensity()
        
        plt.figure(figsize=(20,6))
        plt.stem(self.wavenos, 1 - self.normalized_intensities, linefmt = 'y', markerfmt = 'yo' , bottom=1)
    
    
    
    def smooth_spectra(self, plot=True): 
        
        waves = []
        smoothed_data = []
        
        def gaussian(x, mu, sig):
            return (np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.))))/np.sqrt(2*np.pi*np.power(sig, 2.))

        for i in range(len(self.linelist)):
            wavelength = (1/self.wavenos[i])*1e8
            waves.append(wavelength)
            
        x = np.linspace(np.min(waves) - 0.1 ,np.max(waves) + 0.1, 10000)
        y = np.zeros(x.shape)
        
        for i in range(len(waves)):
            y = y + self.normalized_intensities[i]*gaussian(x, waves[i], waves[i]/(2.355*100000))
        
        if plot:
            plt.figure(figsize=(20,6))
            plt.stem(waves, - self.normalized_intensities)
            plt.plot(x, - (y/max(y)), color = 'black', linewidth = 3)#, label= ('Jmax= ' + str(self.Jmax)))
            #plt.xlim(6613.6, 6614)
            
            plt.legend()
        
        def area_under_curve(x,y):
          
            sum_of_areas = 0
            for i in range(1, len(x)):
                h = x[i] - x[i-1]
                sum_of_areas += h * (y[i-1] + y[i]) / 2
        
            #return sum_of_areas 

        # print('area under curve is  ' + str(area_under_curve(x, y)))  
        # print('sum of normalized intenisities is  ' + str(np.sum(self.normalized_intensities)))
    
        return y



start = timeit.default_timer() 
 
a = Rotational_Spectra(0.0111, 0.010767, 3, 3, temperature=2, origin=15120, Jmax=10)
#plt.title('condition a')
a.get_linelist()
a.smooth_spectra()

stop = timeit.default_timer()
print('Jmax:  ', a.Jmax)
print('Time: ', stop - start)  


'''Previously Used codes'''
#%%
# '''>>>>>>>>> Offset <<<<<<<<<<<<<<'''
# offset = (coronene['position']-a.linelist['wavenos'])
# a.linelist['offset'] = offset
# plt.scatter(a.linelist['ground_J'], a.linelist['offset'])
# plt.title('offset in wavenumber')
# plt.xlabel('J')
# plt.ylabel('offset')
# plt.show()


# difference = a.normalized_intensities - (coronene['strength']/max(coronene['strength']))
# a.linelist['difference'] = difference
# plt.scatter(a.linelist['ground_J'], a.linelist['difference'])
# plt.title('Normalized intensity - Normlaized Pgopher strength')
# plt.xlabel('J')
# plt.ylabel('difference')
# plt.show()


# ratio = a.intensities / coronene['strength']
# a.linelist['ratio'] = ratio
# plt.scatter(a.linelist['ground_J'], a.linelist['ratio'])
# plt.title('Normalized intensity / Normlaized Pgopher strength')
# plt.xlabel('J')
# plt.ylabel('ratio')
# plt.show()

# a = Rotational_Spectra(0.0111, 0.005552, 3, 3, temperature=2, origin=15120)
# b = Rotational_Spectra(0.0111, 0.005552, 3, 3, temperature=3, origin=15120)
# c = Rotational_Spectra(0.0111, 0.005552, 3, 3, temperature=5, origin=15120)
# d = Rotational_Spectra(0.0111, 0.005552, 3, 3, temperature=7, origin=15120)

# molecules = (a, b, c, d)

# for m in molecules:
#     m.get_linelist()
#     m.smooth_spectra()
#a.linelist.to_excel(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\coronene at 15120.xlsx",  index=None)

# coronene = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\My DIBs research\Pgopher practice plots\Coronene P,Q and R Branches\coronene line lists pgopher\coronene Jmax 7.txt", delim_whitespace=(True))

# plt.figure(figsize=(20,6))
# plt.stem(coronene['position'], ( coronene['strength']/max(coronene['strength'])), linefmt = 'b', markerfmt = 'bo')

