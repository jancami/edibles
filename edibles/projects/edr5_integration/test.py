from edibles.utils.edibles_oracle import EdiblesOracle
from edibles.utils.edibles_spectrum import EdiblesSpectrum
import matplotlib.pyplot as plt
import numpy as np

pythia = EdiblesOracle()
filelist = pythia.getFilteredObsList(object=['HD 147889'], MergedOnly=False, Wave=6707)

print(filelist)

for file in filelist:
    sp = EdiblesSpectrum(file)
    wrange = [6700, 6720]
    print('min', np.nanmin(sp.bary_wave))
    print('mmaxin', np.nanmax(sp.bary_wave))
    sp.getSpectrum(xmin=wrange[0], xmax=wrange[1])
    wave=sp.bary_wave
    flux=sp.bary_flux
    plt.plot(wave, flux)
    plt.show()

