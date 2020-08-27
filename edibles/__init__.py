import os

if 'EDIBLES_DATARELEASE' in os.environ:
    DATARELEASE = os.environ['EDIBLES_DATARELEASE']
else:
    DATARELEASE = 'DR4'

if 'EDIBLES_DATADIR' in os.environ:
    DATADIR = os.environ['EDIBLES_DATADIR']
else:
    DATADIR = '/data/DR4'

if 'EDIBLES_PYTHONDIR' in os.environ:
    PYTHONDIR = os.environ['EDIBLES_PYTHONDIR'] + '/edibles'
else:
    PYTHONDIR = os.path.dirname(__file__)

__version__ = '0.1'
