'''Code to find out allowed transitions given Jmax

To do:
>> Modify code to work for different Jmax and Kmax
>> add to rotational spectra.py
'''

import numpy as np
import pandas as pd

Jmax = 300

# #here Kmax = Jmax (i.e all K allowed)


# P_branch_Js = list((range(1,Jmax+1)))
# all_P_branch_Js = []
# for j in P_branch_Js:
#     for i in range(j):
#       all_P_branch_Js.append(j)
    
# P_branch_Jprimes = []
# for j in all_P_branch_Js:
#     if j != 0:
#         P_branch_Jprimes.append(j-1)
        
# pP_Branch_K = []        
# for j in P_branch_Js:
#     stages = list(range(0,j))
#     for i in stages:
#         pP_Branch_K.append(j-i)
        

# pP_Branch_Kprime = []
# for k in pP_Branch_K:
#     pP_Branch_Kprime.append(k-1)
        
# rP_Branch_K = []        
# for j in P_branch_Js:
#     stages = list(range(0,j))
#     stages.sort(reverse=True)
#     for i in stages:
#         rP_Branch_K.append(i)


# rP_Branch_Kprime = []
# for k in rP_Branch_K:
#     rP_Branch_Kprime.append(k+1)
    

# '''Q Branch'''

# Q_branch_Js = list((range(0,Jmax+1)))

# all_Q_branch_Js = []
# for j in Q_branch_Js:
#     # if j ==0:
#     #     all_Q_branch_Js.append(j)
#     if j!= 0:
#         for i in range(j):
#           all_Q_branch_Js.append(j)
      
# Q_branch_Jprimes = []
# for j in all_Q_branch_Js:
#         Q_branch_Jprimes.append(j)
        
# pQ_Branch_K = []        
# for j in Q_branch_Js:
#     stages = list(range(0,j))
#     for i in stages:
#         pQ_Branch_K.append(j-i)
        

# pQ_Branch_Kprime = []
# for k in pQ_Branch_K:
#     pQ_Branch_Kprime.append(k-1)
    
# rQ_Branch_K = []        
# for j in Q_branch_Js:
#     stages = list(range(0,j))
#     stages.sort(reverse=True)
#     for i in stages:
#         rQ_Branch_K.append(i)


# rQ_Branch_Kprime = []
# for k in rQ_Branch_K:
#     rQ_Branch_Kprime.append(k+1)
    
        

        
# '''R Branch'''
        
# R_branch_Js = list((range(0,Jmax)))
# all_R_branch_Js = []
# for j in R_branch_Js:
#     if j ==0:
#         all_R_branch_Js.append(j)
#     elif j!= 0:
#         for i in range(j+1):
#           all_R_branch_Js.append(j)
              
# R_branch_Jprimes = []
# for j in all_R_branch_Js:
#     if j <= Jmax-1:
#         R_branch_Jprimes.append(j+1)
        
# pR_Branch_K = []        
# for j in R_branch_Js:
#     stages = list(range(0,j+1))
#     # if j!= 0:
#     for i in stages:
#         pR_Branch_K.append(j-(i-1))
    

# pR_Branch_Kprime = []
# for k in pR_Branch_K:
#     pR_Branch_Kprime.append(k-1)
    
# rR_Branch_K = []        
# for j in R_branch_Js:
#     stages = list(range(0,j+1))
#     stages.sort(reverse=True)
#     for i in stages:
#         rR_Branch_K.append(i)


# rR_Branch_Kprime = []
# for k in rR_Branch_K:
#     rR_Branch_Kprime.append(k+1)
    
        



# Allowed_Js = (all_P_branch_Js*2) + (all_Q_branch_Js*2) + (all_R_branch_Js*2)
# Allowed_Jprimes = (P_branch_Jprimes*2) + (Q_branch_Jprimes*2) + (R_branch_Jprimes*2)
# Allowed_Ks = pP_Branch_K + rP_Branch_K + pQ_Branch_K + rQ_Branch_K +  pR_Branch_K + rR_Branch_K
# Allowed_Kprimes = pP_Branch_Kprime + rP_Branch_Kprime + pQ_Branch_Kprime + rQ_Branch_Kprime + pR_Branch_Kprime + rR_Branch_Kprime

# # print(len(Allowed_Js))
# # print(len(Allowed_Jprimes))
# # print(len(Allowed_Ks))
# # print(len(Allowed_Kprimes))

# columns = {'J' : Allowed_Js,'Jprime': Allowed_Jprimes, 'K' : Allowed_Ks, 'Kprime' : Allowed_Kprimes}
# selection_rules = pd.DataFrame(columns)

# selection_rules['delta_J'] = selection_rules['Jprime'] - selection_rules['J']
# selection_rules['delta_K'] = selection_rules['Kprime'] - selection_rules['K']

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

# print(selection_rules.to_string())

'''======================================================================='''

'''Combinations'''
#%%   
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
combinations = pd.DataFrame(columns)



label = []

for i in range(len(combinations['ground_J'])):
    if combinations['excited_J'][i] - combinations['ground_J'][i] == -1 and combinations['excited_K'][i] - combinations['ground_K'][i] == -1:
        label.append('pP')
    if combinations['excited_J'][i] - combinations['ground_J'][i] == -1 and combinations['excited_K'][i] - combinations['ground_K'][i] == 1:
        label.append('rP')
    if combinations['excited_J'][i] - combinations['ground_J'][i] == 0 and combinations['excited_K'][i] - combinations['ground_K'][i] == -1:
        label.append('pQ')
    if combinations['excited_J'][i] - combinations['ground_J'][i] == 0 and combinations['excited_K'][i] - combinations['ground_K'][i] == 1:
        label.append('rQ')
    if combinations['excited_J'][i] - combinations['ground_J'][i] == 1 and combinations['excited_K'][i] - combinations['ground_K'][i] == -1:
        label.append('pR')
    if combinations['excited_J'][i] - combinations['ground_J'][i] == 1 and combinations['excited_K'][i] - combinations['ground_K'][i] == 1:
        label.append('rR')

combinations['label'] = label

ground_Js = combinations['ground_J']
excited_Js = combinations['excited_J']
ground_Ks = combinations['ground_K']
excited_Ks = combinations['excited_K']
# delta_J = combinations['delta_J']
# delta_K = combinations ['delta_K']

#combinations.to_excel(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\J=300_trial.xlsx", index=False)

print(combinations.to_string())

# print('----------------------')
# print('Jmax is  ' + str(Jmax))
# print('length of combinations  ' + str(len(combinations)))
# endc = timeit.default_timer()
# print('---------------')
# print('>>>> combnination calclulation takes   ' + str(endc-startc) + '  sec')


    
