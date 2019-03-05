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
from edibles.frequency_convert import convert_to_nu
from edibles.fit.make_grid import make_grid
# from edibles.test_voigt_integral import normalize


# ===========
# Main script
# ===========

# set params
alpha = 0.1
gamma = 0.1
delta_v = 1000
x_min = 7663
x_max = 7671
cent = 7667
n_piece = 3


# generate wavelength grid with resolving power delta_v (R = c/delta_v)
R = cst.c.value / delta_v
x_nonbroad = make_grid(x_min, x_max, resolution=R)
wave = np.array(x_nonbroad)




# generate voigt data with specified parameters
flux_norm = voigt_math(wave, cent, alpha, gamma)

nu, nu0 = convert_to_nu(wave, cent)
print(nu0)

# normalize the flux values
# flux_norm = normalize(wave, flux, show=True)

# flux_norm = -flux_norm + 1

# check x arrays
if np.array_equal(nu, wave) is not True:
	print(len(wave))
	print(len(nu))
	print('Bad x resolution!')


# Generate the continuum data
x_spline, y_spline = generate_continuum(nu, flux_norm, delta_v, n_piece)


# check x arrays
if np.array_equal(nu, x_spline) is not True:
	print('Bad x resolution!')



# ==========
# C*V
# ==========
cxv = y_spline * flux_norm


# ==========
# Residuals
# ==========

# resid = flux_norm - cxv


# plot
plt.plot(nu, flux_norm, 'k.', markersize='1', label='Data')
plt.xlabel('Frequency')
plt.plot(x_spline, y_spline, label='Spline fit')
plt.plot(nu, cxv, label='C*V')
# plt.plot(nu, resid, label='Residuals')

plt.legend()
plt.show()
