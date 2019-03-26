from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import sys


def read2dfile(path):
	'''
	This function will read a n*2 text file and show a plot

	input:		path/to/file

	Output:		plot

	example file:
		7.662 1.453
		7.662 1.436
		7.662 1.414
			.
		 	.
		 	.

	'''

	array = np.loadtxt(path)

	wave = array[:, 0]
	flux = array[:, 1]

	plt.plot(wave, flux)
	plt.show()


fullCmdArguments = sys.argv
args = fullCmdArguments[1:]
arN = len(sys.argv)

if len(args) != 1:
	print('\nSyntax:	python read2dtxt.py path/to/file\n')

else:
	path = args[0]
	read2dfile(path)
