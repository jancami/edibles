import os
from pathlib import Path
from importlib.resources import files

print(os.environ)

if 'EDIBLES_DATARELEASE' in os.environ:
    DATARELEASE = os.environ['EDIBLES_DATARELEASE']
else:
    DATARELEASE = 'DR5'

if 'EDIBLES_DATADIR' in os.environ:
    DATADIR = Path(os.environ['EDIBLES_DATADIR'])
else:
    DATADIR = Path('~/EDR5')

if 'EDIBLES_PYTHONDIR' in os.environ:
    EDIBLES_PYTHONDIR = Path(os.environ['EDIBLES_PYTHONDIR']) / 'edibles'
    
else:
    EDIBLES_PYTHONDIR = files('edibles')

if 'EDIBLES_OUTPUTDIR' in os.environ:
    EDIBLES_OUTPUTDIR = Path(os.environ['EDIBLES_OUTPUTDIR'])
else:
    EDIBLES_OUTPUTDIR = Path('~/EDIBLES_output')

__version__ = '0.2'
