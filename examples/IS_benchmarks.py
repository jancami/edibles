import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path
from edibles import PYTHONDIR

from edibles.utils.ISLineFitter import ISLineFitter, ISLineModel, measure_snr

# Load the benchmark data and run it through ISLineFitter. Then compare results. 

this_datadir = Path(PYTHONDIR) / "data" / "voigt_benchmarkdata"
filename = this_datadir / "omiper.m95.7698.txt"
v_resolution = 0.56 # km/s
#print(filename)
arrays = np.genfromtxt(filename,skip_header=1)

wave = arrays[:,0]
flux = arrays[:,1]
#print(wave, flux)

# First, create a model corresponding to Dan's literature values. 
# Parameters from Dan Welty for all components: 
v = np.array([10.50, 11.52, 13.45, 14.74, 15.72]) + 0.10
N = np.array([12.50, 10.0, 44.3, 22.5, 3.9])*1e10
b = np.array([0.60, 0.44, .72, .62, .60])
#Dans_model = ISLineModel(v_res=0.56, verbose=0, wave=wave)
#Dans_model.AddLine(lineID='KI_7698', v=v, N=N, b=b, verbose=0)

K_wave = [7698.974]
K_fjj = [3.393e-1]
K_Gamma = [3.820e7]
n_components=5
model = ISLineModel(n_components, v_res=0.56, lam_0=K_wave, fjj=K_fjj, gamma=K_Gamma, verbose=0)
pars_guess = model.guess(V_off=v[0:n_components])
idx_cloud = 0
while idx_cloud < n_components:
    pars_guess["V_off_Cloud%i" % idx_cloud].value = v[idx_cloud]
    pars_guess["b_Cloud%i" % idx_cloud].value = b[idx_cloud]
    pars_guess["N_Cloud%i" % idx_cloud].value = N[idx_cloud]
    idx_cloud = idx_cloud + 1
ymodel=model.eval(x=wave,params=pars_guess)
#ymodel=ISLineModel(v_res=0.56, s)
plt.plot(wave, flux, color='black')
plt.plot(wave, ymodel, color='red', ls='--')
plt.plot(wave, ymodel+.1, color='blue')
plt.xlim([7699.1,7699.5])
plt.show()

#model=ISLineModel(LineID='KI_7698', v=v, N=N, b=b, v_res=0.56, wave=wave)



# fitting
fitter = ISLineFitter(wave, flux, v_resolution=0.56, normalized=True)
known_n_components = 1
fit_result = fitter.fit(species='KI', windowsize=1, n_anchors=4,
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
