import os
from pathlib import Path

if 'EDIBLES_DATARELEASE' in os.environ:
    DATARELEASE = os.environ['EDIBLES_DATARELEASE']
else:
    DATARELEASE = 'DR5'

if 'EDIBLES_DATADIR' in os.environ:
    DATADIR = Path(os.environ['EDIBLES_DATADIR'])
else:
    DATADIR = '/home/Alex/spectra/EDR5'

if 'EDIBLES_PYTHONDIR' in os.environ:
    EDIBLES_PYTHONDIR = Path(os.environ['EDIBLES_PYTHONDIR']) / 'edibles'
    
else:
    EDIBLES_PYTHONDIR = os.path.dirname(__file__)    

__version__ = '0.2'
