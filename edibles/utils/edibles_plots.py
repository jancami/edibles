import numpy as np
import matplotlib.pyplot as plt
import astropy.constants as cst
from edibles import PYTHONDIR
from edibles.utils.edibles_oracle import EdiblesOracle
from edibles.utils.edibles_spectrum import EdiblesSpectrum
from pathlib import Path
import pandas as pd

# Case: a certain wavelength range shows some interesting features in the lab. 
# Prepare a few plots that show the "best" EDIBLES data in that wavelength range
# in a way that allows to see clearly what are DIBs and what are stellar lines. 

wrange = [4675,4695]

pythia = EdiblesOracle(verbose=1)
Targetlist = pythia.getFilteredObjects(Wave=np.median(wrange), EBV_min=.8)
for target in Targetlist:
	print(target)
	ObsList = pythia.getFilteredObsList(object=target, Wave=np.median(wrange), OrdersOnly=True)
	if len(ObsList) == 0:
		print('No Observations found in this range for ', target)
	else: 
		for filename in ObsList:
			print(filename)
			sp = EdiblesSpectrum(filename)
			wave = sp.wave
			flux = sp.flux

			idx = np.where((wave > wrange[0]) & (wave < wrange[1]))
			wave = wave[idx]
			flux = flux[idx]
			norm = flux / np.median(flux)

			plt.plot(wave,norm)
			plt.title(target)
			#plt.xlim(6707.0,6709.5)
			plt.ylim(.7,1.3)
			#plt.plot(wave, fitresult.best_fit, color='b')
		plt.show()