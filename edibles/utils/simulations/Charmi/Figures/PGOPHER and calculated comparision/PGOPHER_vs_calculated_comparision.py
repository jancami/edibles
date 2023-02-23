# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 13:42:50 2022

@author: Charmi Bhatt
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math

# pgopher_ll = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Kerr's conditions\condition_c\kerr's_condition_c_pgopher_linelist.txt", delim_whitespace=(True))
                         
# plt.plot(pgopher_ll['Position'], pgopher_ll['Strength'])

plt.figure(figsize = (15,6))


Cal_smooth = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Figures\PGOPHER and calculated comparision\Cal_conditon_c.txt", delim_whitespace=(True))
#plt.plot(Cal_smooth.iloc[:,0], 1-0.1*(Cal_smooth.iloc[:,1]/max(Cal_smooth.iloc[:,1])), linewidth = 4, color = 'orange', alpha = 0.7, label = 'Calculated')
plt.plot(Cal_smooth.iloc[:,0], (Cal_smooth.iloc[:,1]/max(Cal_smooth.iloc[:,1])), linewidth = 4, color = 'orange', alpha = 0.7, label = 'Calculated')

#plt.xlim(-2.5,2.5)

PG_smooth = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Figures\PGOPHER and calculated comparision\PGO_conditon_c.txt", delim_whitespace=(True))
#PG_smooth = PG_smooth[(PG_smooth.iloc[:,0] <= (max(Cal_smooth.iloc[:,0])) ) & (PG_smooth.iloc[:,0] >= (min(Cal_smooth.iloc[:,0])) ) ]
plt.plot(PG_smooth.iloc[:,0], (PG_smooth.iloc[:,1]/max(PG_smooth.iloc[:,1])), color = 'black', label = 'PGOPHER')

# Calculated_normalized_intenisty = 1-0.1*(Cal_smooth.iloc[:,1]/max(Cal_smooth.iloc[:,1]))
# PGO_Calculated_intenisty = 1-0.1*(PG_smooth.iloc[:,1]/max(PG_smooth.iloc[:,1]))
# Calculated_wavenumber = Cal_smooth.iloc[:,0]
# PGO_wavenumber = PG_smooth.iloc[:,0]


# PGO_vs_cal = np.column_stack([PGO_wavenumber, PGO_Calculated_intenisty]) #, Calculated_wavenumber, Calculated_normalized_intenisty ])
# np.savetxt('PGO_conditon_c.txt', PGO_vs_cal)



plt.xlabel('Wavelength (cm$^{-1}$)', fontsize = 15)
plt.ylabel('Normalized Intenisty', fontsize = 15)
plt.legend()

plt.xticks(fontsize=15)
plt.yticks(fontsize= 15)

plt.gca().xaxis.set_major_locator(plt.MultipleLocator(1))
plt.gca().yaxis.set_minor_locator(plt.MultipleLocator(0.01))
plt.gca().xaxis.set_minor_locator(plt.MultipleLocator(0.2))


