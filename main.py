import numpy as np
import matplotlib.pyplot as plt
from astropy import constants as cst

import sys
import os
path = os.getcwd()
os.chdir('..')
sys.path.append(os.getcwd())
os.chdir(path)

from edibles.continuum_guess import generate_continuum
from edibles.voigtMathematical import voigt_math
from edibles.fit.make_grid import make_grid
from edibles.astro_wrapper import voigt_astro
# from edibles.test_voigt_integral import normalize


# ===========
# Main script
# ===========

# set params
alpha = 0.1
gamma = 0.1
delta_v = 1000
x_min = 5977
x_max = 5983
cent = 5980
n_piece = 3

b_eff=3.5
Gamma=1.9


# generate wavelength grid with resolving power delta_v (R = c/delta_v)
R = cst.c.value / delta_v
x_nonbroad = make_grid(x_min, x_max, resolution=R)
wave = np.array(x_nonbroad)


# generate voigt data with specified parameters
flux_norm = voigt_math(wave, cent, alpha, gamma)


plt.plot(wave, flux_norm, 'grey', markersize='1', label='Data')

# Generate the continuum data
x_spline, y_spline = generate_continuum(wave, flux_norm, delta_v, n_piece)

# check x arrays
if np.array_equal(wave, x_spline) is not True:
	print('Bad x resolution!')


plt.plot(x_spline, y_spline, label='Spline fit')

# ==========
# C*V
# # ==========
# cxv = y_spline * flux_norm

# plt.plot(wave, cxv, label='C*V')



# ==========
# Astro
# ==========



x2, y2 = voigt_astro(wave, cent, b_eff, Gamma)

plt.plot(x2, y2, 'red', label='astro')



# plot

plt.xlabel('Frequency')
plt.legend()
plt.show()
