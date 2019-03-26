from __future__ import print_function
import numpy as np 
from edibles.decompose import decompose
import edibles.fit.avoigt as fit
import edibles.fit.avoigt_fit as ft
import edibles.fit.make_grid as mg

# ============================================================================
# code from voigt_test 3
# ============================================================================
grid = mg.make_grid(5886, 5894, resolution=80000)
yy = 1 + fit.voigt(grid, lambda_peak=5890, b_eff=3.47, log_N=12.843, gamma=6.064e+07, osc_freq=0.631, resolving_power=80000)
flux = yy
noise = np.random.normal(0, 0.02, len(flux))
flux = flux + noise
# interpolate into fixed grid
wave = np.linspace(5887, 5893, 1000)
flux = np.interp(wave, grid, flux)

obj = ft.fitSpectrum()
obj.afit(wave, flux, [5890], lines=['NaI_5891'], cheb_order=1, resolving_power=80000)

# ============================================================================
# decompose code
# ============================================================================
params, error = decompose(wave, obj.fit_parm, obj.fit_err, obj.fit_norm, print_data=True, plot=False)
print(params)

# ============================================================================
# code from voigt_test 7
# ============================================================================
grid = mg.make_grid(5886, 5894, resolution=80000)
cont_const = 0.05 * (grid-5890) + 12.3
yy1 = 1 + fit.voigt(grid, lambda_peak=5890, b_eff=3.47, log_N=12.34, gamma=6.064e+07, osc_freq=0.631, resolving_power=80000)
yy2 = 1 + fit.voigt(grid, lambda_peak=5890.12, b_eff=4.13, log_N=11.17, gamma=6.064e+07, osc_freq=0.631, resolving_power=80000)
flux = cont_const + yy1 + yy2
noise = np.random.normal(0, 0.02, len(flux))
flux = flux + noise
# interpolate into fixed grid
wave = np.linspace(5887, 5893, 1000)
flux = np.interp(wave, grid, flux)

obj = ft.fitSpectrum()
obj.afit(wave, flux, [[5890], [5890.1]], lines=['NaI_5891'], cheb_order=1, resolving_power=80000)


# ============================================================================
# decompose code
# ============================================================================
params, error = decompose(wave, obj.fit_parm, obj.fit_err, obj.fit_norm, print_data=True)
print(params)
