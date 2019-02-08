import numpy as np
import scipy.interpolate
import matplotlib.pyplot as plt

import sys
import os
path = os.getcwd()
os.chdir('..')
sys.path.append(os.getcwd())
os.chdir(path)

from edibles.voigtMathematical import voigt_math


def generateContinuum(x, y, n_sections=4, move_x=True, move_y=False):
	'''
	This function fits a continuum to data separated into n sections
	using a cubic spline.

	INPUT:
	x: [ndarray] wavelength grid
	y: [ndarray] flux values
	n_sections: [int, default=4] evenly split data into n sections
		(n+1 points) one point will be on either end of the data
	move_x: [bool, default=True] change x values of points to fit
	move_y: [bool, default=False] change y values of points to fit

	OUTPUT:
	x_cont: [ndarray] wavelength grid points for fit continuum
	y_cont: [ndarray] flux value points for fit continuum
	'''

	x_sections = np.split(x, n_sections)
	y_sections = np.split(y, n_sections)

	# print(x_sections)



	x_value_array = [x[0]]
	y_value_array = [y[0]]
	for i in range(len(x_sections)):

		x_value = x_sections[i][-1]
		y_value = y_sections[i][-1]
		x_value_array.append(x_value)
		y_value_array.append(y_value)


	# cont_values = np.column_stack((x_value_array, y_value_array))

	return x_value_array, y_value_array



# Example

x = np.linspace(0, 1000, 1000)


y = voigt_math(x, 0.5, 0.5)


generateContinuum(x, y)













# scipy.interpolate.CubicSpline()

# plt.plot()