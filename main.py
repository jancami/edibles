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


# ===========
# Main script
# ===========

# set params
alpha = 0.1
gamma = 0.1
delta_v = 100.0
x_min = 7665
x_max = 7669
cent = 7667
n_piece = 3


# generate wavelength grid with resolving power delta_v (R = c/delta_v)
R = cst.c.value / delta_v
x_nonbroad = make_grid(x_min, x_max, resolution=R)
wave = np.array(x_nonbroad)


# generate voigt data with specified parameters
flux = -voigt_math(wave, cent, alpha, gamma) + 1


# Generate the continuum data
x_spline, y_spline = generate_continuum(wave, flux, delta_v, n_piece)


# check x arrays
if np.array_equal(wave, x_spline) is not True:
	print('Bad x resolution!')



# ==========
# C*V
# ==========
cxv = y_spline * flux


# ==========
# Residuals
# ==========

resid = flux - cxv


# plot
plt.plot(wave, flux, 'k.', markersize='1', label='Data')
plt.plot(x_spline, y_spline, label='Spline fit')

plt.plot(wave, cxv, label='C*V')
plt.plot(wave, resid, label='Residuals')

plt.legend()
plt.show()
