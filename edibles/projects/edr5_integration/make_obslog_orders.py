"""
Making an obslog for the EDR5 orders. 
"""

from astropy.io import fits
from pathlib import Path
from pprint import pprint
from edibles.projects.edr5_integration import dr5_io

edr5_dir = Path('/home/Alex/spectra/EDR5')
spec_dir = Path('/home/Alex/spectra/EDR5/orders')
file_list = list(spec_dir.glob('*.fits'))

for file in file_list[:1]:
    with fits.open(file) as hdul:
        hdr = hdul[0].header
        data = hdul[1].data
        pprint(hdul[1].header)
    
    pprint(hdr)
    object_name = hdr['OBJECT']
    ra = hdr['RA']
    dec = hdr['DEC']
    date_obs = hdr['DATE-OBS']
    wave_setting, arm_setting = dr5_io.get_wave_path(hdr)


