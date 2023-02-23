# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 16:07:33 2022

@author: Charmi Bhatt
"""

import numpy as np
import pandas as pd

# values= np.array([0.9769, 1.0307, 1.0692, 1.0615, 1, 1.0461, 1.1230, 1.0230, 1.0538, 0.9769, 0.9769, 0.99230])
# uncertainties = np.array([0.0692, 0.0384, 0.0461, 0.0461, 0.0692, 0.0384, 0.0538, 0.0230, 0.0384, 0.0230, 0.0384, 0.0538])


sightline_data = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\average_B\heather_sightline_data.txt", delim_whitespace=True)

Vrq = sightline_data['Vrq']
Vqp = sightline_data['Vqp']
Vrq_err = sightline_data['Vrq_Error']
Vqp_err = sightline_data['Vqp_Error']

values = (Vrq - Vqp)/1.27
uncertainties = np.sqrt((Vrq_err**2 + Vqp_err**2))/1.27


weights = 1/(uncertainties**2)
values_times_weights = values*weights
weighted_average = np.sum(values_times_weights)/np.sum(weights)
weighted_average_err = 1/(np.sqrt(np.sum(weights)))

print(weighted_average)
print(weighted_average_err)