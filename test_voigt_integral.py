import numpy as np
import matplotlib.pyplot as plt
from astropy import constants as cst

import sys
import os
path = os.getcwd()
os.chdir('..')
sys.path.append(os.getcwd())
os.chdir(path)

from edibles.voigtMathematical import voigt_math
from edibles.fit.make_grid import make_grid


# set params
alpha = 0.1
gamma = 0.1
delta_v = 10
x_min = 7665
x_max = 7669
cent = 7667
n_piece = 3


# generate wavelength grid with resolving power delta_v (R = c/delta_v)
R = cst.c.value / delta_v
x_nonbroad = make_grid(x_min, x_max, resolution=R)
wave = np.array(x_nonbroad)


# generate voigt data with specified parameters
flux = voigt_math(wave, cent, alpha, gamma)

# compute area under curve
area = np.trapz(flux, wave)
print('Area under curve: ', area)


# loop until normalized (maybe)
i = 0
while area != 1.0:
	print('Normalizing')
	flux = flux / area
	area = np.trapz(flux, wave)

	print('Area under curve: ', area)
	print(area.is_integer())
	i = i + 1
	if i >1000:
		print('Area under curve could not be normalized to one!')
		continue


# plot
plt.plot(wave, flux, 'k.', markersize='1')

plt.legend()
plt.show()
