'''Code to find out allowed transitions given Jmax

To do:
>> Modify code to work for different Jmax and Kmax
>> add to rotational spectra.py
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

Jmax = 3

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
        
qP_Branch_K = []        
for j in P_branch_Js:
    stages = list(range(0,j))
    for i in stages:
        qP_Branch_K.append(j-i)
        

qP_Branch_Kprime = []
for k in qP_Branch_K:
    qP_Branch_Kprime.append(k)
        

    

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
        
qQ_Branch_K = []        
for j in Q_branch_Js:
    stages = list(range(0,j))
    for i in stages:
        qQ_Branch_K.append(j-i)
        

qQ_Branch_Kprime = []
for k in qQ_Branch_K:
    qQ_Branch_Kprime.append(k)
    

        

        
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
        
qR_Branch_K = []        
for j in R_branch_Js:
    stages = list(range(0,j+1))
    # if j!= 0:
    for i in stages:
        qR_Branch_K.append(j-(i-1))
    
qR_Branch_Kprime = []        
for j in R_branch_Js:
    stages = list(range(0,j+1))
    # if j!= 0:
    for i in stages:
        qR_Branch_Kprime.append(j-(i-1))
    

# qR_Branch_Kprime = []
# for j in qR_Branch_K:
#     qR_Branch_Kprime.append(j-(i-1))
    

    
        



Allowed_Js = (all_P_branch_Js) + (all_Q_branch_Js) + (all_R_branch_Js)
Allowed_Jprimes = (P_branch_Jprimes) + (Q_branch_Jprimes) + (R_branch_Jprimes)
Allowed_Ks = qP_Branch_K +  qQ_Branch_K +   qR_Branch_K 
Allowed_Kprimes = qP_Branch_Kprime  + qQ_Branch_Kprime  + qR_Branch_Kprime 

# print(len(Allowed_Js))
# print(len(Allowed_Jprimes))
# print(len(Allowed_Ks))
# print(len(Allowed_Kprimes))

columns = {'ground_J' : Allowed_Js,'excited_J': Allowed_Jprimes, 'ground_K' : Allowed_Ks, 'excited_K' : Allowed_Kprimes}
selection_rules = pd.DataFrame(columns)

selection_rules['delta_J'] = selection_rules['excited_J'] - selection_rules['ground_J']
selection_rules['delta_K'] = selection_rules['excited_K'] - selection_rules['ground_K']

# label = []

# for i in range(len(selection_rules['J'])):
#     if selection_rules['Jprime'][i] - selection_rules['J'][i] == -1 and selection_rules['Kprime'][i] - selection_rules['K'][i] == -1:
#         label.append('pP')
#     if selection_rules['Jprime'][i] - selection_rules['J'][i] == -1 and selection_rules['Kprime'][i] - selection_rules['K'][i] == 1:
#         label.append('rP')
#     if selection_rules['Jprime'][i] - selection_rules['J'][i] == 0 and selection_rules['Kprime'][i] - selection_rules['K'][i] == -1:
#         label.append('pQ')
#     if selection_rules['Jprime'][i] - selection_rules['J'][i] == 0 and selection_rules['Kprime'][i] - selection_rules['K'][i] == 1:
#         label.append('rQ')
#     if selection_rules['Jprime'][i] - selection_rules['J'][i] == 1 and selection_rules['Kprime'][i] - selection_rules['K'][i] == -1:
#         label.append('pR')
#     if selection_rules['Jprime'][i] - selection_rules['J'][i] == 1 and selection_rules['Kprime'][i] - selection_rules['K'][i] == 1:
#         label.append('rR')

# selection_rules['Label'] = label

print(selection_rules.to_string())

'''======================================================================='''

ground_B = 0.011
delta_B = 0.5
delta_C = 0.5
origin = 0

ground_C = ground_B/2
delta_C = delta_C
excited_B = ground_B + ((delta_B/100)*ground_B)
excited_C = ground_C + ((delta_C/100)*ground_C)

global combinations
ground_Js = selection_rules['ground_J']
excited_Js = selection_rules['excited_J']
ground_Ks = selection_rules['ground_K']
excited_Ks = selection_rules['excited_K']

linelist = selection_rules
   
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
for J,K in zip(excited_Js, excited_Ks):
            excited_E = excited_B*J*(J + 1) + (excited_C - excited_B)*(K**2) # - ((-2*excited_C*zeta))*K + excited_C**2
    
            excited_Es.append(excited_E)
            
linelist['excited_Es'] = excited_Es

wavenos = []
for i in range(len(linelist.index)):
    wavenumber = origin + excited_Es[i] - ground_Es[i]
    wavenos.append(wavenumber)
    
linelist['wavenos'] = wavenos

print(ground_Es)
print(excited_Es)
print(wavenos)

for i in range(len(wavenos)):
    plt.axhline(y=wavenos[i], xmin = 0.2, xmax=0.6, color='green')
    
#plt.xlim(0,0.6)
    
    
    
    
    
    