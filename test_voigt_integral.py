import numpy as np
import sys


def normalize(x, y, show=False):
	'''
	This function checks if the area under the curve produced using the x
	and y arrays is equal to one. If it is not equal to one, the y axis is
	recursively normalized until its integral is 1.0.

	INPUT:
	x: [ndarray] input x axis
	y: [ndarray] input y axis
	
	OUTPUT:
	y_norm: [ndarray] output normalized y axis

	'''

	# Check to see if area under curve equals 1.0
	area = np.trapz(y, x)

	if show is True:
		print('Initial area under curve: {:.16f}'.format(area))

	# loop until normalized (maybe)
	i = 0
	while area != 1.0:

		# stop after too long
		if i >1000:
			print('Area under curve could not be normalized to one!')
			sys.exit()

		i = i + 1
		if show is True:
			print('Normalizing attempt number: {}'.format(i))

		# normalize / renormalize
		y = y / area
		area = np.trapz(y, x)

		if show is True:
			print('Area under curve: {:.16f}'.format(area))
			print('area.is_integer() results: {}'.format(area.is_integer()))
			print('')

	y_norm = y

	return y_norm
