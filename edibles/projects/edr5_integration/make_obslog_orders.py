"""
Making an obslog for the EDR5 orders. 
"""

from astropy.io import fits
from pathlib import Path
from pprint import pprint
from edibles.projects.edr5_integration import dr5_io
import numpy as np
import pandas as pd
from importlib.resources import files
from edibles import DATADIR

edr5_dir = DATADIR
spec_dir = DATADIR / 'orders'
file_list = list(spec_dir.glob('*.fits'))
obs_log_path = files('edibles') / 'data/DR5_ObsLog.csv'

out_df = pd.DataFrame()
for file in file_list:
    with fits.open(file) as hdul:
        hdr = hdul[0].header
        data = hdul[1].data
    
    object_name = hdr['OBJECT']
    object_name = object_name.replace(' ', '').replace('HD', 'HD ')
    ra = hdr['RA']
    dec = hdr['DEC']
    date_obs = hdr['DATE-OBS']
    wave_setting, arm_setting = dr5_io.get_wave_path(hdr)
    o_index = file.name.find('O')+1
    order = file.name[o_index:].replace('.fits', '')
    spec = dr5_io.read_spec(file)
    wave_min = np.nanmin(spec[0])
    wave_max = np.nanmax(spec[0])

    file_name = file.relative_to(edr5_dir)

    out_s = pd.Series({'Object': object_name, 'RA': ra, 'DEC': dec, 'DateObs': date_obs, 'Setting': wave_setting, 'Order': order, 'WaveMin': wave_min, 'WaveMax': wave_max, 'Filename': file_name})
    out_df = pd.concat((out_df, out_s), axis=1, ignore_index=True)

out_df = out_df.T

print(out_df)
out_df.to_csv(obs_log_path, index=False)


