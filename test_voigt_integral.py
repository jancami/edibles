import numpy as np
import sys


def normalize(x, y):
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

	print('Initial area under curve: {:.16f}'.format(area))

	# loop until normalized (maybe)
	i = 0
	while area != 1.0:
		i = i + 1
		print('Normalizing attempt number: {}'.format(i))
		y = y / area
		area = np.trapz(y, x)

		print('Area under curve: {:.16f}'.format(area))
		print(area.is_integer())
		if i >1000:
			print('Area under curve could not be normalized to one!')
			sys.exit()

		y_norm = y

	return y_norm
