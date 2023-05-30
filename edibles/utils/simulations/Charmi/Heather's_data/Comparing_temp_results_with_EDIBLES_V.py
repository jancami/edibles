import numpy as np
import pandas as pd
import astropy.constants as const
import matplotlib.pyplot as plt
import timeit
import scipy.stats as ss
from scipy.signal import argrelextrema
from matplotlib import cm
import numpy as np
import pandas as pd
import astropy.constants as const
import matplotlib.pyplot as plt
import timeit
import scipy as sp
import scipy.stats as ss
from lmfit import Model
import csv
import lmfit
from lmfit import minimize, Parameters, report_fit 
import uncertainties as unc
import uncertainties.umath as umath                                                                          

data = pd.read_csv('/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/fitting methods/individual_best_fits.csv') #, encoding = 'latin-1')

PR_sep_value = list(data['PR_sep'])
PR_sep_unc = list(data['PR_sep_unc'])
PR_sep = [unc.ufloat(value, uncertainty) for value, uncertainty in zip(PR_sep_value, PR_sep_unc)]
# print(PR_sep)

B = unc.ufloat(0.0222, 0.0089)
B_minus_one_sigma = 0.022 - 0.0089

B_plus_one_sigma = 0.022 + 0.0089


Temp = []
for i in range(len(PR_sep)):
    T = 0.180 * (PR_sep[i]** 2 / B)
    np.set_printoptions(precision=3)

    # print(T)
    Temp.append(T)
    
print(type(Temp))
#Temp = np.array(Temp)

#np.set_printoptions(formatter={'float': '{:0.3f}'.format})
#print(['{:.3f}'.format(x.nominal_value) for x in Temp])
for i,t in zip(range(len(PR_sep_value)), Temp):
    plus_uncertainity = (0.180*(PR_sep_value[i]**2)) /  B_minus_one_sigma
    plus_uncertainity = t.nominal_value - plus_uncertainity
    print(plus_uncertainity)
    
for i,t in zip(range(len(PR_sep_value)), Temp):
    minus_uncertainity = (0.180*(PR_sep_value[i]**2)) /  B_plus_one_sigma
    minus_uncertainity = t.nominal_value - minus_uncertainity
    print(minus_uncertainity )
    
'''PLOTTING'''''


ata = pd.read_csv('/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/fitting methods/individual_best_fits.csv') #, encoding = 'latin-1')
print(data)
data = data.drop(9)
print(data)


T_linear = data['T (linear)']
pattern = r'([\d.]+)\+/-([\d.]+)'
matches = T_linear.str.extract(pattern, expand=True)
values = matches[0].astype(float)
uncertainties = matches[1].astype(float)

# Convert values and uncertainties to lists
T_linear_list = values.tolist()
#T_linear_err = [[data['T_minus']],[data['T_plus']]]

array = pd.concat([data['T_minus'], data['T_plus']], axis=0)

# Reshape the array to (2, n)
T_linear_err = array.values.reshape(2, -1)



T_symmetric_list = data['T']
#T_symmetric_list.drop(10)

T_symmetric_err = data['T_unc']
#T_symmetric_err.drop(10)

plt.errorbar(T_symmetric_list, T_linear_list, xerr = T_symmetric_err, yerr = T_linear_err, fmt='o', color='black',
             ecolor='lightgray', elinewidth=3, capsize=0 )
plt.xlabel('T (symmetric)')
plt.ylabel('T (Linear)')
plt.title('Comparing results with EDIBLES V')

      
