from edibles import DATADIR
from edibles import PYTHONDIR
from edibles.utils.edibles_spectrum import EdiblesSpectrum
from edibles.utils.edibles_oracle import EdiblesOracle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from PyAstronomy import pyasl
import astropy.constants as cst
from io import BytesIO


# labfilename = '/data/python/devel/141_REMPI_472-438.txt'

labfilename = "/Users/charmibhatt/Desktop/Edibles_All/Edibles_Scripts/edibles/utils/Labplots/141_REMPI_472-438.txt"
# The file contains wavenumbers (in vacuum) and intensities. Read those in as pandas dataframe. 
labspec = pd.read_csv(labfilename, delimiter=',')
# Add a column to contain the (air) wavelength in AA. 
labspec['wavelength'] = pyasl.vactoair2(1e1 * labspec['wave'],mode='ciddor')
normint = labspec.int - np.median(labspec.int)
labspec['norm'] = np.exp(-0.02 * normint / np.max(normint))


# no observattions found: 27778, 54239

objects = ["HD 22951", "HD 23016", "HD 23180", "HD 24398", "HD 41117", "HD 45314", "HD 61827", "HD 63804", "HD 73882", "HD 80558", "HD 112272", "HD 147084", "HD 147683", "HD 147888", "HD 147889", "HD 147933", "HD 148184", "HD 148937", "HD 149404", "HD 149757", "HD 152248", "HD 154043", "HD 154368", "HD 161056", "HD 165319", "HD 167971", "HD 169454", "HD 170740", "HD 170938", "HD 171957", "HD 179406", "HD 183143", "HD 185418", "HD 185859", "HD 186745", "HD 186841", "HD 203532", "HD 210121"]

V_rad_data = pd.read_csv('/Users/charmibhatt/Desktop/Edibles_All/Edibles_Scripts/edibles/utils/Labplots/V_rad_from_Haoyu_readable.csv')

# Entire wavelength range is roughly 4350-4720AA, but strongest band near 4689 DIB. 
# Match from Tim is in HD 147889, so let's try to reproduce that first.... 
image_data_list = []
plotrange = [4670,4700]
pythia = EdiblesOracle()

for object in objects: 
    List = pythia.getFilteredObsList(object=[object], OrdersOnly=True, Wave=4689)

    for i, filename in zip(range(len(List)), List): 
        sp = EdiblesSpectrum(filename)
        sp.getSpectrum(np.min(sp.raw_wave)+1,np.max(sp.raw_wave)-1)

        wave = sp.bary_wave
        #search for V_rad from the csv file provided
        row = V_rad_data.loc[V_rad_data['object'] == object]
        print(object)
        v_rad = row['V_rad'].values[0]
        print(v_rad)
        #radial velocity correction
        wave_rest = wave / (1+v_rad/cst.c.to("km/s").value)


        flux = np.clip(sp.bary_flux, 0, None) 

        bool_keep = (wave_rest > plotrange[0]) & (wave < plotrange[1])
        plotwave = wave_rest[bool_keep]
        plotflux = flux[bool_keep]

        #normalization
        normflux = plotflux / np.median(plotflux)
        plt.figure()
        plt.xlim(4680, 4695)
        ylim = [np.min(normflux), np.max(normflux)]
        plt.ylim(ylim)
        plt.plot(plotwave, normflux)
        plt.title(filename)

        plt.plot(labspec.wavelength, labspec.norm)
        #plt.show()

        
        save_plot_as = f'labplot_v_rad_{v_rad}_{sp.target}_{sp.datetime}_{i}.png'
        plt.savefig(save_plot_as, format='png')
        
        
        #Rescale lab spectrum to plot range
        #dynrange = ylim[1]-ylim[0]
    

#plt.plot(labspec.wavelength,labspec.int)
#plt.show()


