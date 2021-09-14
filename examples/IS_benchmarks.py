import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path
from edibles import PYTHONDIR

from edibles.utils.ISLineFitter import ISLineFitter, ISLineModel

# Load the benchmark data and run it through ISLineFitter. Then compare results. 

this_datadir = Path(PYTHONDIR) / "data" / "voigt_benchmarkdata"
filename = this_datadir / "omiper.m95.7698.txt"
v_resolution = 0.56 # km/s
#print(filename)
arrays = np.genfromtxt(filename,skip_header=1)

wave = arrays[:,0]
flux = arrays[:,1]
#print(wave, flux)


# fitting
fitter = ISLineFitter(wave, flux, v_resolution=0.56)
known_n_components = 5
fit_result = fitter.fit(species='KI', windowsize=1.5, n_anchors=4,
                        known_n_components=known_n_components, WaveMax=7700, WaveMin=7695)

# results
print(fit_result.fit_report())

# plot best fit model
plt.plot(fitter.wave2fit, fitter.flux2fit)
plt.plot(fitter.wave2fit, fit_result.best_fit)
plt.xlabel("Overall")
plt.show()

# plot best fit model of each components
comps = fit_result.eval_components(x=wave)
continuum = comps["cont"]
plt.plot(fitter.wave2fit, fitter.flux2fit)
singe_component = ISLineModel(1, lam_0=fitter.air_wavelength,
                              fjj=fitter.oscillator_strength,
                              gamma=fitter.gamma,
                              v_res=fitter.v_res)
for i in range(known_n_components):
    pars = singe_component.make_params(b_Cloud0=fit_result.params["b_Cloud%i" % (i)].value,
                                       N_Cloud0=fit_result.params["N_Cloud%i" % (i)].value,
                                       V_off_Cloud0=fit_result.params["V_off_Cloud%i" % (i)].value)
    plt.plot(fitter.wave2fit, singe_component.eval(params=pars, x=fitter.wave2fit), label="Comp %i" % (i))
plt.xlabel("Individual Components")
plt.legend()
plt.show()
