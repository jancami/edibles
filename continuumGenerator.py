import numpy as np
from scipy.interpolate import CubicSpline
from scipy.stats import chisquare
import matplotlib.pyplot as plt

import sys
import os
path = os.getcwd()
os.chdir('..')
sys.path.append(os.getcwd())
os.chdir(path)

from edibles.voigtMathematical import voigt_math
from edibles.fit.make_grid import make_grid


def generateContinuum(x, y, n_sections=4, move_x=False, move_y=True):
	'''
	This function fits a continuum to data separated into n sections
	where the x and y-values are the median of each section	using a cubic spline

	INPUT:
	x: [ndarray] wavelength grid
	y: [ndarray] flux values
	n_sections: [int, default=4] evenly split data into n sections
	move_x: [bool, default=False] change x values of points to fit
	move_y: [bool, default=True] change y values of points to fit

	OUTPUT:
	x_cont: [ndarray] wavelength grid points for fit continuum
	y_cont: [ndarray] flux value points for fit continuum
	'''

	x_sections = np.array_split(x, n_sections)
	y_sections = np.array_split(y, n_sections)

	x_value_array = []
	y_value_array = []
	for i in range(len(x_sections)):
		x_value = np.median(x_sections[i])
		y_value = np.median(y_sections[i])
		x_value_array.append(x_value)
		y_value_array.append(y_value)


	spline = CubicSpline(x_value_array, y_value_array)


	x_spline = np.linspace(np.min(x), np.max(x), (len(x)))
	y_spline = spline(x_spline)

	# # plot median points of each section
	# plt.plot(x_value_array, y_value_array, 'ko')


	return x_spline, y_spline


# Example

alpha, gamma = 0.1, 0.1

wave = np.linspace(-3, 3, 1000)
flux = -voigt_math(wave, alpha, gamma)

# Generate the spline curve
x_spline, y_spline = generateContinuum(wave, flux)

plt.plot(wave, flux, marker='o', markersize='3', label='Data')
plt.plot(x_spline, y_spline, label='Spline fit')
plt.legend()
plt.show()


plt.show()
