import os

if 'EDIBLES_DATADIR' in os.environ:
    DATADIR = os.environ['EDIBLES_DATADIR']
else:
    DATADIR = '/data/DR4'


if 'EDIBLES_DATARELEASE' in os.environ:
    DATADIR = os.environ['EDIBLES_DATARELEASE']
else:
    DATADIR = 'DR4'


if 'EDIBLES_PYTHONDIR' in os.environ:
    DATADIR = os.environ['EDIBLES_PYTHONDIR']
else:
    DATADIR = '~/python/edibles/'

__version__ = ''