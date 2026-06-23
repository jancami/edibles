# -*- coding: utf-8 -*-
"""
Created on Mon Jan 10 13:38:12 2022

@author: gabir
"""

import os
import matplotlib.pyplot as plt
import pandas as pd
from edibles.utils.functions import vac2air_ciddor

os.chdir(r'C:\Users\gabir\edibles\edibles\data\Labdata\CRDS')
# # #changed the current directory to the current path where PENTACENE will be found 



labspec = pd.read_csv('PENTACENE.DAT', delim_whitespace= True)
#using the vacuum to air function from util.functions to convert from vaccum to air 
#set the air values to the list labspec
labspec['wavelength'] = vac2air_ciddor(1e8/labspec['wno'])
#use matplotlib to plot labspec wavelength by labspec integer 
plt.plot(labspec.wavelength, labspec.int)

plt.xlabel('Wavelength (angstroms)') #wno would be wavenumber
plt.ylabel('intensity ')
plt.title('Pentacene')


















































