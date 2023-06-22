"""
>find wa   overplot for different molecules(objects)

"""

import pandas as pd
import matplotlib.pyplot as plt
import astropy.constants as const
from astropy.convolution import convolve, Gaussian1DKernel
import numpy as np


#pgopher input

coronene = pd.read_csv(r'C:/Users/Charmi Bhatt/OneDrive/Desktop/My DIBs research/Pgopher practice plots/Coronene P,Q and R Branches/coronene line lists pgopher/coronene at 3kelvin.txt', delim_whitespace = True)
coronene = coronene.sort_values(by=['position'], ascending=True)

ground_J = coronene['J'] #coronene[(coronene['J'] <= 10)] 
excited_J = coronene['Jprime'] #coronene[(coronene['Jprime'] <= 10)] #
ground_K = coronene['K'] #coronene[(coronene['Jprime'] <= 10)] #
excited_K = coronene['Kprime']



JandKs = {'ground_J' : ground_J , 'excited_J' : excited_J, 'ground_K' : ground_K, 'excited_K' : excited_K }

linelist = pd.DataFrame(JandKs)

# print(len(JandKs))
# print(linelist)

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
        
        self.linelist = linelist
        self.ground_energies = []
        self.excited_energies = []
        self.wavenos = []
        self.HL_factors = []
        self.BD_factors = []
        self.intensities = []
        self.normalized_intensities = []
        self.BD_with_partition = []
        

    def ground_energy_levels(self):
        
        #ground_energies = []
        
        for J,K in zip(ground_J, ground_K):
            ground_E = self.ground_B*J*(J + 1) + (self.ground_C - self.ground_B)*(K**2)
            self.ground_energies.append(ground_E)
            
        #self.ground_energies= ground_energies
        #self.linelist['ground_energies'] = self.ground_energies
            
        return self.ground_energies
            
    def excited_energy_levels(self):
        
        #excited_energies = []
        
        for J,K in zip(excited_J, excited_K):
            excited_E = self.excited_B*J*(J + 1) + (self.excited_C - self.excited_B)*(K**2)
            self.excited_energies.append(excited_E)
            
        # BJ(J+1) + (C-B)K**2   
        # self.excited_energies = excited_energies
        # self.linelist['excited_energies'] = self.excited_energies
            
        return self.excited_energies
            
    def wavenumber(self):

        #wavenos = []

        for i in range(len(linelist.index)):
            waveno = self.origin + self.excited_energy_levels()[i] - self.ground_energy_levels()[i]
            self.wavenos.append(waveno)
            
        #self.wavenos = wavenos
        #self.linelist['wavenumber'] = self.wav       
        
        return self.wavenos
    

    def honl_london_factors(self):
        
        self.HL_factors = []
                
        for J,K, Jprime, Kprime in zip(ground_J, ground_K, excited_J, excited_K):
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
        for J,K,E in zip(ground_J, ground_K, self.ground_energy_levels()):
            if K == 0:
                boltzmann_equation = (2*((2*J) + 1))*(np.exp((-h * c * E) / (k*self.T)))
            else:
                boltzmann_equation = (1*((2*J) + 1))*(np.exp((-h * c * E) / (k*self.T)))
                
            self.BD_factors.append(boltzmann_equation)
            
        # self.BD_factors = BD_factors
        # self.linelist['BD_factors'] = self.BD_factors
                
        return self.BD_factors
    
    def with_partition_function(self):
        
        total = []
        sum = 0
        
        for i in self.Boltzmann():
            sum = sum + i
            total = sum
            
        self.BD_with_partition = total
    
        return self.BD_with_partition
    
    def intensity(self):
        
        for i in range(len(linelist.index)):
            strength = (self.honl_london_factors()[i] * self.Boltzmann()[i])
            self.intensities.append(strength)
            
        
        return self.intensities 
    
    def normalized_intensity(self):     
        
        #self.intensity()
        self.normalized_intensities = self.intensities / max(self.intensities)
        
        return self.normalized_intensities
        
    def get_linelist(self):
        
        self.linelist['groundE'] = self.ground_energy_levels()
        self.linelist['excitedE'] = self.excited_energy_levels()
        self.linelist['wavenos'] = self.wavenumber()
        self.linelist['HL_factors'] = self.honl_london_factors()
        self.linelist['BD_factors'] = self.Boltzmann()
        self.linelist['intensities'] = self.intensity()
        self.linelist['normalized intensities'] = self.normalized_intensity()
        #self.linelist['offset'] = self.find_offset()
        
        return self.linelist
    
    def plot_rotational_spectra(self):
        
        # self.wavenumber()
        # self.normalized_intensity()
        
        plt.figure(figsize=(20,6))
        plt.stem(self.wavenos, - self.normalized_intensities, linefmt = 'y', markerfmt = 'yo') #, bottom=1)
    
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
            y = y + self.normalized_intensities[i]*gaussian(x, waves[i], waves[i]/(2.355*110000))
        
        if plot:
            plt.figure(figsize=(20,6))
            plt.stem(waves, self.normalized_intensities)
            plt.plot(x, y/max(y), color = 'black', linewidth = 5, label= ('Temperature = ' + str(self.T) + 'K'))
            plt.xlim(6613.6, 6614)
            plt.legend()
        
        def area_under_curve(x,y):
          
           sum_of_areas = 0
           for i in range(1, len(x)):
               h = x[i] - x[i-1]
               sum_of_areas += h * (y[i-1] + y[i]) / 2
        
           return sum_of_areas 

        print('area under curve is  ' + str(area_under_curve(x, y)))  
        print('sum of normalized intenisities is  ' + str(np.sum(self.normalized_intensities)))
    
        return y
    
a = Rotational_Spectra(0.0111, 0.005552, 3, 3, temperature=2, origin=15120)

a.get_linelist()
a.plot_rotational_spectra()



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
