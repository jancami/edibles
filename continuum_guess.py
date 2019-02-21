import numpy as np
from scipy.interpolate import CubicSpline
from astropy import constants as cst

import sys
import os
path = os.getcwd()
os.chdir('..')
sys.path.append(os.getcwd())
os.chdir(path)

from edibles.fit.make_grid import make_grid



def generate_continuum(x, y, delta_v=1.0, n_piece=None):
	'''
	This function fits a continuum to data separated into n sections
	where the x and y-values are the median of each section	using a cubic spline

	INPUT:
	x:       [ndarray] wavelength grid
	y:       [ndarray] flux values
	delta_v: [float] desired resolution of continuum (in m/s)
	n_piece: [int, default=4] evenly split data into n sections

	OUTPUT:
	x_cont:  [ndarray] wavelength grid points for fit continuum
	y_cont:  [ndarray] flux value points for fit continuum
	
	Graphic:

	X------|-:    |    :-|-----X
	|      |  :   |   :  |     |
	|      |   :  X  :   |     |
	|      |    :-|-:    |     |
	|______|______|______|_____|

	<----piece---><---piece---->
	<-sec-><-sec-><-sec-><-sec->

	'''
	
	
       # check n_piece param
	if n_piece is None: n_piece = 2

	# split x & y arrays into n_piece*2 sections
	x_sections = np.array_split(x, n_piece*2)
	y_sections = np.array_split(y, n_piece*2)

	# initialize list of points to spline fit
	x_points = [np.min(x)]
	y_points = [np.median(y_sections[0])]


	# loop through every other section (1, 3, 5...)
	# make n_piece+1 points to fit a spline through
	# create a spline point on edge of each piece
	for i in range(1, len(x_sections), 2):


		# set x_point 
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


	# make resolving power delta_v (R = c/delta_v)
	R = cst.c.value / delta_v
	x_nonbroad = make_grid(np.min(x), np.max(x), resolution=R)
	x_spline = np.array(x_nonbroad)

	spline = CubicSpline(x_points, y_points)
	y_spline = spline(x_spline)
	# plt.plot(x_points, y_points, 'kx', markersize='8', label='Points')

	return x_spline, y_spline



