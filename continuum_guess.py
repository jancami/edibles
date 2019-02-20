import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
from astropy.modeling.models import Voigt1D
from astropy import constants as cst

import sys
import os
path = os.getcwd()
os.chdir('..')
sys.path.append(os.getcwd())
os.chdir(path)

from edibles.fit.make_grid import make_grid



def generate_continuum(x, y, n_piece=None):
	'''
	This function fits a continuum to data separated into n sections
	where the x and y-values are the median of each section	using a cubic spline

	INPUT:
	x: [ndarray] wavelength grid
	y: [ndarray] flux values
	n_piece: [int, default=4] evenly split data into n sections
	move_x: [bool, default=False] change x values of points to fit
	move_y: [bool, default=True] change y values of points to fit

	OUTPUT:
	x_cont: [ndarray] wavelength grid points for fit continuum
	y_cont: [ndarray] flux value points for fit continuum
	'''
	
	
        # check n_piece param
	if n_piece is None: n_piece = 2

	
	# make resolving power 1 m/s (R = c/delta_v)
	R = cst.c.value / 1.0
	x_nonbroad = make_grid(np.min(x), np.max(x), resolution=R)
	x_spline = np.array(x_nonbroad)


	x_sections = np.array_split(x, n_piece*2)
	y_sections = np.array_split(y, n_piece*2)

	# initialize list of points to spline fit
	x_points = [np.min(x)]
	y_points = [np.median(y_sections[0])]


	for i in range(1, len(x_sections), 2):

		x_point = np.max(x_sections[i])

		# create span of points for median
		# check if last section
		if x_point == np.max(x):
			span = y_sections[i]
			y_point = np.median(span)

		else:
			span = np.append(y_sections[i], y_sections[i+1])
			y_point = np.median(span)

		x_points.append(x_point)
		y_points.append(y_point)

	spline = CubicSpline(x_points, y_points)
	y_spline = spline(x_spline)
	# plt.plot(x_points, y_points, 'kx', markersize='8', label='Points')

	return x_spline, y_spline


# # Example
# alpha, gamma = 0.1, 0.1

# wave = np.linspace(7666, 7668, 1000)
# # flux = -voigt_math(wave, alpha, gamma)
# voigt_func = Voigt1D(x_0=7667, amplitude_L=10, fwhm_L=0.1, fwhm_G=0.2)
# flux = -1 * voigt_func(wave)

# # Generate the spline curve
# x_spline, y_spline = generate_continuum(wave, flux)

# plt.plot(wave, flux, marker='o', markersize='3', label='Data')
# plt.plot(x_spline, y_spline, label='Spline fit')
# plt.legend()
# plt.show()
