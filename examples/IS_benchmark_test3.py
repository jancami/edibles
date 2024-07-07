# Purpose:
# To test if Chi-Sqr reported lmfit is the same as we understand, and
# to determine if we should use SNR as weights in the fitting.

# Method:
# Use Dan Welty's data and try fit with 5 components.
# Spoon-feeding: use published results as initial guess of parameters
# Calculate Chi-Sqr as Sigma (Res * SNR) **2, and compare to result.chisqr

# Results
# ChiSqr reported in the result class of lmfit is the same as my calculation.
# The lmfit package is using weights in a wired way, and we should use SNR as weights


##########################
# Code starts here
##########################
# imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from pathlib import Path
from edibles import PYTHONDIR
from edibles.utils.ISLineFitter import ISLineFitter, ISLineModel, measure_snr

# Load the benchmark data
this_datadir = Path(PYTHONDIR) / "data" / "voigt_benchmarkdata"
filename = this_datadir / "omiper.m95.7698.txt"
v_resolution = 0.56 # km/s
arrays = np.genfromtxt(filename, skip_header=1)
wave = arrays[:, 0]
flux = arrays[:, 1]

# Measuring result from DW+2001 for omiPer
known_Voffs = [10.50, 11.52, 13.45, 14.74, 15.72]
known_Ns = [12.5, 10.0, 44.3, 22.5, 3.9] # * 10**10
known_bs = [0.6, 0.44, 0.72, 0.62, 0.3]
K_wave = [7698.974]
K_fjj = [3.393e-1]
K_Gamma = [3.820e7]

# Check data
# data is ~ 2.0 AA wide.
# peak at 7699.35 AA, shoulder at 7699.15 and 7699.5 AA
# SNR ~ 120
SNR = np.max(measure_snr(wave, flux, block_size=0.2))
plt.plot(wave, flux)
plt.plot(wave, np.ones_like(wave)*(1 + 1/SNR))
plt.plot(wave, np.ones_like(wave)*(1 - 1/SNR))
plt.ylim([(1 - 5/SNR), (1 + 5/SNR)])
plt.show(block=False)
plt.pause(3)
plt.close()
print("SNR", SNR)
print("N_points", len(wave))

# use a window of ~1.05 AA. This should not change the result
# Data and model to fit
n_components = 5
wave2fit = wave[(wave >= 7698.8) & (wave <= 7699.85)]
flux2fit = flux[(wave >= 7698.8) & (wave <= 7699.85)]
model2fit = ISLineModel(n_components, v_res=0.56,
                        lam_0=K_wave, fjj=K_fjj, gamma=K_Gamma,
                        verbose=0)

# Fit!
pars_guess = model2fit.guess(V_off=known_Voffs[0:n_components])
idx_cloud = 0
while idx_cloud < n_components:
    pars_guess["V_off_Cloud%i" % idx_cloud].value = known_Voffs[idx_cloud]
    pars_guess["b_Cloud%i" % idx_cloud].value = known_bs[idx_cloud]
    pars_guess["N_Cloud%i" % idx_cloud].value = known_Ns[idx_cloud] * 10**10
    idx_cloud = idx_cloud + 1

result = model2fit.fit(data=flux2fit,
                       params=pars_guess,
                       x=wave2fit,
                       weights=np.ones_like(wave2fit) * SNR)

# Calculate ChiSqr
res = flux2fit - result.best_fit
print("ChiSqr by lmfit: %.2f" % result.chisqr)
print("ChiSqr by calculation: %.2f" % np.sum((res * SNR)**2))
print("Reduced ChiSqr (by calculation): %.2f" % (np.sum((res * SNR)**2)/(len(wave2fit)-3*5)))

# Plot result
fig = plt.figure(figsize=(10, 6.5))
plt.gcf().subplots_adjust(hspace=0)
spec = gridspec.GridSpec(ncols=1, nrows=2,
                         height_ratios=[4, 1])
# top panel for data and model
ax0 = fig.add_subplot(spec[0])
plt.gca().xaxis.set_visible(False)
plt.scatter(wave2fit, flux2fit)
plt.plot(wave2fit, result.best_fit, color="orange")
plt.ylabel("Fitted Model")

# bottom panel for residual and SNR
ax2 = fig.add_subplot(spec[1])
plt.plot(wave2fit, res)
plt.plot(wave2fit, np.ones_like(wave2fit)/SNR, color="r", linestyle="--")
plt.plot(wave2fit, -1 * np.ones_like(wave2fit)/SNR, color="r", linestyle="--")
plt.ylabel("Residual")
plt.xlabel("Wavelength")
plt.show()

