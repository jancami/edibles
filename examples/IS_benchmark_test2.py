# Purpose:
# To test if changing window size would affect the final result, especially regarding
# the number of components to be kept.

# Method:
# Use Dan Welty's data and try fit with 3, 4, 5, and 6 components.
# Spoon-feeding: use published results as initial guess of parameters
# DW reported 5 components and the 6th components is to test how AIC/BIC might work
# Change window size in each test

# Results
#  Test  WindowSize   n_AIC   n_BIC
#    0      0.35        5       5
#    1      0.75        5       4
#    2      1.05        5       4
#    3      1.45        5       4

# Comments
# 1. Reduced Chi-Sqr is much smaller in the first test, something in the "continuum region"
#    may not have been addressed.
# 2. Fitting results depends on window size and method used (AIC or BIC) to a certain degree



##########################
# Code starts here
##########################
# imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
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
known_Voffs = [10.50, 11.52, 13.45, 14.74, 15.72, 13.00]
known_Ns = [12.5, 10.0, 44.3, 22.5, 3.9, 5.0] # * 10**10
known_bs = [0.6, 0.44, 0.72, 0.62, 0.3, 0.3]
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
plt.pause(5)
plt.close()
print("SNR", SNR)
print("N_points", len(wave))

# See how window size will affect fitting results
for test in [0,1,2,3]:
    wave2fit = wave[(wave >= 7699.15 - test*0.2) & (wave <= 7699.5 + test*0.2)]
    flux2fit = flux[(wave >= 7699.15 - test*0.2) & (wave <= 7699.5 + test*0.2)]

    AIC, BIC, ChiSqr = [], [], []
    V_all, b_all, N_all = [[]]*4, [[]]*4, [[]]*4
    for n_components in [3, 4, 5, 6]:
        model2fit = ISLineModel(n_components, v_res=0.56,
                                lam_0=K_wave, fjj=K_fjj, gamma=K_Gamma,
                                verbose=0)
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

        AIC.append(result.aic)
        BIC.append(result.bic)
        ChiSqr.append(result.chisqr)

        V_all[n_components - 3] = [result.params["V_off_Cloud%i" % idx_cloud].value
                                 for idx_cloud in range(n_components)]
        b_all[n_components - 3] = [result.params["b_Cloud%i" % idx_cloud].value
                                   for idx_cloud in range(n_components)]
        N_all[n_components - 3] = [result.params["N_Cloud%i" % idx_cloud].value / 10**10
                                   for idx_cloud in range(n_components)]

        plt.plot(wave2fit, flux2fit)
        plt.plot(wave2fit, result.best_fit)
        plt.xlabel("Test %i, %i components" % (test, n_components))
        plt.show(block=False)
        plt.pause(5)
        plt.close()

    # Print results from a single test (which weight to use)
    print("="*15 + " Statistic Result " + "="*15)
    print("test: %i" % test)
    print("Wavelength window between: %.2f and %.2f" % (np.min(wave2fit), np.max(wave2fit)))
    print("AIC: ", AIC)
    print("BIC: ", BIC)
    print("ChiSqr: ", ChiSqr)
    print("Reduced ChiSqr: ", [item/(len(wave2fit) - n_components*3) for item in ChiSqr])
    print("\n")
    print("="*15 + " Fitting Result " + "="*15)

    print("V_known: ", known_Voffs)
    print(", ".join(["%.2f" % item for item in V_all[0]]))
    print(", ".join(["%.2f" % item for item in V_all[1]]))
    print(", ".join(["%.2f" % item for item in V_all[2]]))
    print(", ".join(["%.2f" % item for item in V_all[3]]))

    print("N_known: ", known_Ns)
    print(", ".join(["%.2f" % item for item in N_all[0]]))
    print(", ".join(["%.2f" % item for item in N_all[1]]))
    print(", ".join(["%.2f" % item for item in N_all[2]]))
    print(", ".join(["%.2f" % item for item in N_all[3]]))

    print("b_known: ", known_bs)
    print(", ".join(["%.2f" % item for item in b_all[0]]))
    print(", ".join(["%.2f" % item for item in b_all[1]]))
    print(", ".join(["%.2f" % item for item in b_all[2]]))
    print(", ".join(["%.2f" % item for item in b_all[3]]))

    print("\n")

