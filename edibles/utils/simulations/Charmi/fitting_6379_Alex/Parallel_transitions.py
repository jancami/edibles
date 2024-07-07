# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 16:12:34 2023

@author: alexr
"""

import pandas as pd
import timeit 


Jmax = 300

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
    columns = {'ground_J' : Allowed_Js,'excited_J': Allowed_Jprimes, 'ground_K' : Allowed_Ks, 'excited_K' : Allowed_Kprimes}
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



combinations = allowed_parallel_transitions(Jmax)
# print(combinations)
