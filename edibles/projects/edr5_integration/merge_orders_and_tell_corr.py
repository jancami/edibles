"""
Merging telluric corrections and uncorrected orders to single files, using the foormat of CRIRES+.
Before merging, the orders are cropped.
"""

from importlib.resources import files
import pandas as pd
from edibles import DATADIR
import numpy as np
from astropy.io import fits
from util_functions import setting_dependent_crop, crop_spectrum
from pprint import pprint

obs_file = files('edibles') / 'data/DR5_ObsLog.csv'
obs_list = pd.read_csv(obs_file)
print(obs_list)


obs_times = obs_list['DateObs'].unique()

print(obs_times)


for obs_time in obs_times:
    sub_df = obs_list[obs_list['DateObs'] == obs_time]
    sub_df = sub_df.sort_values(by=['Order'])
    
    for i, row in sub_df.iterrows():
        file = DATADIR / row['Filename']
        with fits.open(file) as hdul:
            # pprint(hdul[0].header)
            # pprint(hdul[1].header)
            # pprint(hdul[1].data)
            data = hdul[1].data

            wave = data['WAVE']
            flux = data['FLUX']
            error = data['ERROR']
            flat = data['FLAT']

        tell_file = DATADIR / 'tell_corr' / file.name

        if tell_file.is_file():
            with fits.open(tell_file) as hdul:
                

