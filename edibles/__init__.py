import os

if os.environ['EDIBLES_DATADIR']:
    DATADIR = os.environ['EDIBLES_DATADIR']
else:
    DATADIR = '/data/DR4'


if os.environ['EDIBLES_DATARELEASE']:
    DATADIR = os.environ['EDIBLES_DATARELEASE']
else:
    DATADIR = 'DR4'


if os.environ['EDIBLES_PYTHONDIR']:
    DATADIR = os.environ['EDIBLES_PYTHONDIR']
else:
    DATADIR = '~/python/edibles/'
