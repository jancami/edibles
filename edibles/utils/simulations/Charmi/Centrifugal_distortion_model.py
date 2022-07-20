# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 15:29:57 2022

@author: Charmi Bhatt
"""

import pandas as pd
import numpy as np
from astropy import constants as const
import matplotlib.pyplot as plt

'''July 17, conclusion: 
    
    Yay! it works... offset is found to be in the range of e-12 for all branches!! (tested at origin = 15120 and kmax = 5)
    
    >> add it to oop code
    
    '''

T = 2
origin = 15120

ground_B = 0.0111
excited_B = 0.010767
delta_B = excited_B - ground_B

ground_C = 0.005552
excited_C = 0.00538544
delta_C = excited_C - ground_C


ground_Dj = 6.6e-7
ground_Djk = 7.7e-7
ground_Dk = 8.8e-5

excited_Dj = 6.1e-7
excited_Djk = 7.1e-7
excited_Dk = 8.1e-7




ground_energy_pP = []
ground_energy_rP = []
ground_energy_pQ = []
ground_energy_rQ = []
ground_energy_pR = []
ground_energy_rR = []

excited_energy_pP = []
excited_energy_rP = []
excited_energy_pQ = []
excited_energy_rQ = []
excited_energy_pR = []
excited_energy_rR = []

deviation_pP = []
deviation_rP = []
deviation_pQ = []
deviation_rQ = []
deviation_pR = []
deviation_rR = []


ground_energies = (ground_energy_pP, ground_energy_rP, ground_energy_pQ, ground_energy_rQ, ground_energy_pR, ground_energy_rR)
excited_energies = (excited_energy_pP, excited_energy_rP, excited_energy_pQ, excited_energy_rQ, excited_energy_pR, excited_energy_rR)
deviation = (deviation_pP, deviation_rP, deviation_pQ, deviation_rQ, deviation_pR, deviation_rR)
label = ['pP', 'rP', 'pQ', 'rQ', 'pR', 'rR']




#dj + djk + dk
#coronene = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\My DIBs research\Pgopher practice plots\Coronene P,Q and R Branches\coronene line lists pgopher\coronene_centrifugal_test.txt", delim_whitespace= True)

#dj 
#coronene = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\My DIBs research\Pgopher practice plots\Coronene P,Q and R Branches\coronene line lists pgopher\coronene_centrifugal_dj.txt", delim_whitespace=(True))

#kmax = 5
#coronene = pd.read_csv(r"C:\Users\Charmi Bhatt\Desktop\Pgopher practice plots\Coronene P,Q and R Branches\coronene line lists pgopher\coronene_strengthtest2_k5.txt", delim_whitespace= True)

#15120 at kmax= 5
coronene = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\My DIBs research\Pgopher practice plots\Coronene P,Q and R Branches\coronene line lists pgopher\cfdt at 15120 and kmax=5.txt", delim_whitespace= True)





'''>>>>>>>>>>>>>>>>>> Defining P, Q and R Branches <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'''
coronene = coronene.sort_values(by=['J'], ascending=True)

#coronene = coronene[(coronene['K'] != 0)] #different scale factors found for k = 0 and K > 0
#coronene = coronene[(coronene['Jprime'] != 0)]  #to avoid 'division by zero' error

P_Branch = coronene[(coronene['label'].str[1] == "P")]
Q_Branch = coronene[(coronene['label'].str[1] == "Q")]
R_Branch = coronene[(coronene['label'].str[1] == "R")]


pP_Branch = coronene[(coronene['label'].str[1] == "P") & (coronene['label'].str[0] == "p")]
rP_Branch = coronene[(coronene['label'].str[1] == "P") & (coronene['label'].str[0] == "r")]

pQ_Branch = coronene[(coronene['label'].str[1] == "Q") & (coronene['label'].str[0] == "p")]
rQ_Branch = coronene[(coronene['label'].str[1] == "Q") & (coronene['label'].str[0] == "r")]

pR_Branch = coronene[(coronene['label'].str[1] == "R") & (coronene['label'].str[0] == "p")]
rR_Branch = coronene[(coronene['label'].str[1] == "R") & (coronene['label'].str[0] == "r")]
#rR_Branch = rR_Branch.iloc[1: , :]


#print(coronene['Kprime'])



Branches = (pP_Branch, rP_Branch, pQ_Branch, rQ_Branch, pR_Branch, rR_Branch)
number_of_branches = [0,1,2,3,4,5]

lines = [] #number of transitions
index_of_lines = []
for n in number_of_branches:
    length = len(Branches[n])
    lines.append(length)
    

'''wavenumber'''

waveno_pP = []
waveno_rP = []
waveno_pQ = []
waveno_rQ = []
waveno_pR = []
waveno_rR = []

wavenos = (waveno_pP, waveno_rP, waveno_pQ, waveno_rQ, waveno_pR, waveno_rR )

#excited level energy
for n in number_of_branches:
    for Jprime, Kprime in zip(Branches[n]['Jprime'], Branches[n]['Kprime']):
                #Eprime = excited_B*Jprime*(Jprime + 1) + (excited_C - excited_B)*(Kprime**2)
               
                
               #pgopher's
                Eprime = excited_B*Jprime*(Jprime + 1) + (excited_C - excited_B)*(Kprime**2) - (excited_Dj * (Jprime**2) * ((Jprime+1)**2)) - (excited_Djk * Jprime * (Jprime+1) * (Kprime**2)) - (excited_Dk*(Kprime**4))
                
                # print('------------------------------')
                # #print((excited_Dj * (Jprime**2) * ((Jprime+1)**2)))
                # #print((excited_Djk * Jprime * (Jprime+1) * (Kprime**2)))
                # print((excited_Dk*(Kprime**4)))
                # print('------------------------------')
                # #herzberg's                
               # Eprime = (2*excited_B*(Jprime+1)) - (2*excited_Djk*(Kprime**2)*(Jprime+1)) - (4*excited_Dj*((Jprime+1)**3))
                
                excited_energies[n].append(Eprime)
                
#ground level energy
                
for n in number_of_branches:
    for J, K in zip(Branches[n]['J'], Branches[n]['K']):
                  # E = ground_B*J*(J + 1) + (ground_C - ground_B)*(K**2)
                  
                  
                  #pgopher's
                  E = ground_B*J*(J + 1) + (ground_C - ground_B)*(K**2) - (ground_Dj * (J**2) * ((J+1)**2)) - (ground_Djk * J * (J+1) * (K**2)) - (ground_Dk*(K**4))
                  
                  
                  
                  #herzbergs
                  #E = (2*ground_B*(J+1)) - (2*ground_Djk*(K**2)*(J+1)) - (4*ground_Dj*((J+1)**3))
                  
                  ground_energies[n].append(E)
                  
#waveno = excited - ground
                  
for n,l in zip(number_of_branches, lines):   
    #if n == 0:
        index = list(range(l))
        for i in index:
            waveno = origin + excited_energies[n][i] - ground_energies[n][i]
            wavenos[n].append(waveno)
            
            
'''>>>>>>>>>>>>>>>>>> Defining deviation <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'''

            
for n,l in zip(number_of_branches, lines):   
    #if n ==0:
        index = list(range(l))
        for i in index:
            offset = (Branches[n]['position'].iloc[i] - wavenos[n][i])
            deviation[n].append(offset)

print(min(deviation[5]))
print(max(deviation[5]))

'''Plotting Deviation'''

for n in number_of_branches:
    #if n == 0:
        plt.scatter(Branches[n]['J'], deviation[n], label = label[n])
        plt.xlabel('J')
        plt.ylabel('offset')
        plt.legend()
        #plt.ylim(0,4)
        for i, l in enumerate(Branches[n]['label']):
                plt.annotate(Branches[n]['K'].iloc[i], (Branches[n]['J'].iloc[i], deviation[n][i]))
                plt.legend() 
        
        
# for n in number_of_branches:
#     #if n == 0:    
#         fig, ax = plt.subplots(figsize=(25,6))
#         ax.set(title = "Coronene  " + label[n] + "  Branch",
#                 xlabel = "Wavenumber",
#                 ylabel = "Normalized Inetensity")
            
#         ax.stem(wavenos[n], Branches[n]['strength']/(max(Branches[n]['strength'])), linefmt=('red'), markerfmt= 'ro', label = 'calculated')
#         ax.stem(Branches[n]['position'], Branches[n]['strength']/(max(Branches[n]['strength'])),linefmt=('blue'), markerfmt= 'bo', label = 'Pgopher')
#         plt.legend()
#         plt.show()
            
        
        


        
for n in number_of_branches:
    #if n == 5:
            table_contents = { 'Branch' : label[n],
                              'J' : Branches[n]['J'],                        
                              'K' : Branches[n]['K'],
                              'waveno' : wavenos[n],
                              'pgopher': Branches[n]['position'],
                              'deviation': max(deviation[n])}

   
            linelist = pd.DataFrame(data = table_contents )
            
            # print('--------------------------------')
            # print(label[n])
            # print('--------------------------------')
            # print(linelist.to_string())
            # print('--------------------------------')
            
# linelist.to_excel(r"C:\Users\Charmi Bhatt\Desktop\rR.xlsx",  index=None) #sheet_name = 'rP') #sep='\t', mode='a')
       
# for J in pP_Branch['J']:
#     if J == 2:
#         print(pP_Branch['K'])


                                                   
                       
                           
                

            