#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 13:54:55 2023

@author: charmibhatt
"""

import numpy as np
import pandas as pd
import astropy.constants as const
import matplotlib.pyplot as plt
import timeit
# import scipy
# import scipy.stats as ss
from lmfit import Model
import csv
# import lmfit
import numba as nb
from pathlib import Path
from lmfit import Parameters


T = [84.45, 95.03, 98.26, 115.26, 97.56, 84.71, 98.41, 88.62, 95.90, 87.57, 102.36, 86.49]
T_err = [6.82, 8.72, 9.76, 22.16, 10.30, 9.688, 9.39, 6.69, 9.31, 7.14, 10.59, 9.59]
sightlines = ['23180', '24398', '144470', '147165' , '147683', '149757', '166937', '170740', '184915', '185418', '185859', '203532']



x = np.linspace(1,12, 12)
plt.errorbar(x = x, y = T, yerr = T_err, fmt='o', color='black',
             ecolor='lightgray', elinewidth=3, capsize=0 )


plt.xticks(x, sightlines, rotation = 75)
plt.ylabel('Temperature')
plt.show()
T = [95.67, 93.59, 99.89, 105.50, 96.33, 82.07, 91.02]
T_err = [13.09, 10.72, 13.49, 15.36, 13.46, 8.135, 11.32]
sightlines = ['144217', '144470', '145502', '147165', '149757', '179406', '184915']



x = np.linspace(1,7,7)
plt.errorbar(x = x, y = T, yerr = T_err, fmt='o', color='black',
             ecolor='lightgray', elinewidth=3, capsize=0 )


plt.xticks(x, sightlines, rotation = 75)
plt.ylabel('Temperature')
plt.title('Cami et al 2004 data')


