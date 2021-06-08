import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path
from edibles import PYTHONDIR
from edibles.utils.edibles_spectrum import EdiblesSpectrum

# Load the benchmark data and run it through ISLineFitter. Then compare results. 

this_datadir = Path(PYTHONDIR) / "data" / "voigt_benchmarkdata"
filename = this_datadir / "omiper.m95.7698.txt"
v_resolution = 0.56 # km/s
#print(filename)
arrays = np.genfromtxt(filename,skip_header=1)

wave = arrays[:,0]
flux = arrays[:,1]
#print(wave, flux)

