from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
import edibles.fit.mpfit_3 as mpfit
from edibles.functions.continuum_guess import generate_continuum
from scipy.optimize import curve_fit
from astropy import constants as cst
from edibles.fit.make_grid import make_grid
from edibles.functions.voigtMathematical import voigt_math






alpha = 0.0576265588185308
gamma = 0.00048255778745462673
delta_v = 1000
x_min = 5977
x_max = 5983
cent = 5980
n_piece = 4

b_eff=3.47
Gamma=6.064e7


# generate wavelength grid with resolving power delta_v (R = c/delta_v)
R = cst.c.value / delta_v
x_nonbroad = make_grid(x_min, x_max, resolution=R)
wave = np.array(x_nonbroad)


# generate voigt data with specified parameters
flux_norm = voigt_math(wave, cent, alpha, gamma)

# Generate the continuum data
y_spline, y_points = generate_continuum((wave, flux_norm), delta_v=delta_v, n_piece=n_piece)

# plot
plt.plot(wave, flux_norm, markersize='1', label='Data')
plt.plot(wave, y_spline)
print(len(wave), len(flux_norm))
plt.show()




# add some noise to the data and try to fit the data generated beforehand
initial_guess = (delta_v)

data_noisy = flux_norm + 0.02*np.random.normal(size=len(wave))
plt.plot(wave, data_noisy)
plt.show()

# popt, pcov = curve_fit(generate_continuum, (wave,flux_norm), data_noisy, p0 = initial_guess)


# print(popt)