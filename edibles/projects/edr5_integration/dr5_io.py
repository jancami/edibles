from astropy.io import fits
from pathlib import Path
import numpy as np
from edibles.projects.edr5_integration import transformations

def read_spec(spec_path:Path, bary_corr=False):
    if spec_path.parent.name == 'tell_corr':
        with fits.open(spec_path) as f:
            data = f[1].data
            hdr = f[0].header

        wave = transformations.angstrom_vac_to_air(data['lambda']*10000)
        spec = np.array([wave, data['cflux']])

    else:
        with fits.open(spec_path) as f:
            data = f[1].data
            hdr = f[0].header

        spec = np.array([data['WAVE'], data['FLUX']])

    if bary_corr:
        spec[0] = transformations.doppler_shift_wl(spec[0], hdr['BARYCORR'])
    
    return spec


def read_spec_wn(spec_path:Path):
    spec = read_spec(spec_path)
    spec[0] = transformations.angstrom_to_wavenumber(spec[0])
    
    return spec


def read_star_name_obs_time(spec_path:Path):
    with fits.open(spec_path) as f:
        hdr = f[0].header

    return hdr['OBJECT'], hdr['DATE-OBS']


def get_wave_path(hdr):
    try:
        wave = hdr['ESO INS GRAT1 WLEN']
    except KeyError:
        wave = hdr['ESO INS GRAT2 WLEN']

    setting = hdr['ESO INS PATH'].lower()
    wave = int(wave)

    return wave, setting


def get_filter_name(hdr):
    try:
        filter_name = hdr['ESO INS FILT2 NAME']
    except KeyError:
        filter_name = hdr['ESO INS FILT3 NAME']

    return filter_name
