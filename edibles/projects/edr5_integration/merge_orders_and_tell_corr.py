"""
Merging telluric corrections and uncorrected orders to single files, using the format of CRIRES+.
Before merging, the orders are cropped.
"""

from importlib.resources import files
import pandas as pd
from edibles import DATADIR
from edibles.projects.edr5_integration import transformations
import numpy as np
from astropy.io import fits
from util_functions import setting_dependent_crop, crop_spectrum
from pprint import pprint
import matplotlib.pyplot as plt

obs_file = files('edibles') / 'data/DR5_ObsLog.csv'
obs_list = pd.read_csv(obs_file)
obs_list = obs_list.loc[obs_list['Order'] != 'ALL']
print(obs_list)

out_dir = DATADIR / 'combined'


obs_times = obs_list['DateObs'].unique()

print(obs_times)

ins_paths = ['blue', 'redl', 'redu']

for obs_time in obs_times:
    sub_df = obs_list[obs_list['DateObs'] == obs_time]
    sub_df = sub_df.sort_values(by=['Order'])

    setting = sub_df['Setting'].iloc[0]

    if setting in [346, 437]:
        out_shape = (5, 0)
    elif setting in [564, 860]:
        out_shape = (8, 0)
    

    for ins_path in ins_paths:
        out_spec = np.array([]).reshape(*out_shape)
        subsub_df = sub_df.loc[sub_df['Filename'].str.contains(ins_path), :]

        # print(subsub_df)
        if subsub_df.empty:
            continue
        for i, row in subsub_df.iterrows():
            file = DATADIR / row['Filename']
            with fits.open(file) as hdul:
                print(file)
                hdu_0 = hdul[0]
                data = hdul[1].data

                wave = data['WAVE']
                flux = data['FLUX']
                error = data['ERROR']
                flat = data['FLAT']
                order = np.full(len(wave), row['Order'])

                iter_spec = np.array([wave, flux, error, flat, order], dtype=float)

                tell_file = DATADIR / 'tell_corr' / file.name

                if tell_file.is_file():
                    try:
                        with fits.open(tell_file) as hdul_tell:
                            data_tell = hdul_tell[1].data

                    except OSError:
                        print(f'Telluric file corrupted: {tell_file}')
                        tell_spec = np.full((3, len(wave)), np.nan)

                    else:
                        m_wave = data_tell['mlambda'] * 1e4
                        m_wave = transformations.angstrom_vac_to_air(m_wave)
                        cflux = data_tell['cflux']
                        m_trans = data_tell['mtrans']
                        tell_spec = np.array([m_wave, cflux, m_trans], dtype=float)
                        hdul[1]['M_WAVE'] = m_wave
                        hdul[1]['CFLUX'] = cflux
                        hdul[1]['MTRANS'] = m_trans
                        hdul.writeto(file)

                else:
                    tell_spec = np.full((3, len(wave)), np.nan)

            if setting in [564, 860]:
                iter_spec = np.concatenate((iter_spec, tell_spec), axis=0)
            
            cl_ang = setting_dependent_crop(iter_spec, setting)
            iter_spec = crop_spectrum(iter_spec, cl_ang[0], cl_ang[1])

            out_spec = np.concatenate((out_spec, iter_spec), axis=1)

        # print(out_spec.shape)
        # plt.plot(out_spec[0], out_spec[1])

        # if setting in [564, 860]:
        #     plt.plot(out_spec[0], out_spec[6])
        #     plt.plot(out_spec[0], out_spec[7])

        # plt.show()

        # Make fits file
        col1 = fits.Column(name='WAVE', format='D', array=out_spec[0])
        col2 = fits.Column(name='FLUX', format='D', array=out_spec[1])
        col3 = fits.Column(name='FLUX_ERROR', format='D', array=out_spec[2])
        col4 = fits.Column(name='FLAT', format='D', array=out_spec[3])
        col5 = fits.Column(name='ORDER', format='J', array=out_spec[4])

        if setting in [564, 860]:
            col6 = fits.Column(name='M_WAVE', format='D', array=out_spec[5])
            col7 = fits.Column(name='CFLUX', format='D', array=out_spec[6])
            col8 = fits.Column(name='MTRANS', format='D', array=out_spec[7])

            coldefs = fits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8])

        else:
            coldefs = fits.ColDefs([col1, col2, col3, col4, col5])

        hdu = fits.BinTableHDU.from_columns(coldefs, name='SCI')

        file = DATADIR / row['Filename']

        with fits.open(file) as hdul:
            hdu_0 = hdul[0]

            out_file = file.name.replace(f'_O{row["Order"]}', '')

            out_hdul = fits.HDUList([hdu_0, hdu])
            out_hdul.writeto(out_dir / out_file, overwrite=True)






        
        





