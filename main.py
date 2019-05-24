import numpy as np
import matplotlib.pyplot as plt
from astropy import constants as cst
from edibles.functions.continuum_guess import generate_continuum
from edibles.functions.voigtMathematical import voigt_math
from edibles.fit.make_grid import make_grid
from edibles.functions.astro_wrapper import voigt_astro
# from edibles.tests.test_voigt_integral import normalize


# ===========
# Main script
# ===========

# set params
alpha = 0.0576265588185308
gamma = 0.00048255778745462673
delta_v = 1000
x_min = 5977
x_max = 5983
cent = 5980
scaling = 5.0
n_piece = 3

b_eff=3.47
Gamma=6.064e7




# some real data:
# NaI_5891       NaI     5891.5835  0.63100000  6.064e+07  22.990


# generate wavelength grid with resolving power delta_v (R = c/delta_v)
R = cst.c.value / delta_v
x_nonbroad = make_grid(x_min, x_max, resolution=R)
wave = np.array(x_nonbroad)


# generate voigt data with specified parameters
flux_norm = voigt_math(wave, cent, alpha, gamma, scaling)


plt.plot(wave, flux_norm, 'grey', markersize='1', label='Data')

# Generate the continuum data
y_spline, y_points = generate_continuum((wave, flux_norm), delta_v=delta_v, n_piece=n_piece)



plt.plot(wave, y_spline, label='Spline fit')

# ==========
# C*V
# # ==========
# cxv = y_spline * flux_norm

# plt.plot(wave, cxv, label='C*V')



# ==========
# Astro
# ==========



y2 = voigt_astro(wave, cent, b_eff, Gamma, scaling)

plt.plot(wave, y2, 'red', linestyle=':', label='astro')
plt.legend()
print(np.max(flux_norm))
print(np.max(y2))


difference = flux_norm - y2
diff_add = np.sum(difference)
print('Total difference between math and voigt:', diff_add)


plt.figure()
plt.plot(wave, difference, label='residual')
plt.xlabel('Frequency')
plt.legend()
plt.show()



