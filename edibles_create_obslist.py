# Set up of the list of FITS files in the data directory with the necessary information that we 
# need to supply to the oracle to work. Essentially, we will recursively list all of the FITS files
# first, then open each of them to read the header and extract the information. Finally, store everything 
# in a text file. 
import os
from astropy.io import fits
import numpy as np
import csv
from edibles_settings import *
from edibles_spectrum import *


print(datadir)
allfitsfiles = []
for path, dirs, files in os.walk(datadir):
	for file in files:
		if file.endswith(".fits"):
			fullfilename = os.path.join(path, file)
			relative_filename = fullfilename[len(datadir):]
			print(relative_filename)
			allfitsfiles.append(relative_filename)

print(len(allfitsfiles))

class obsfile: 
	def __init__ (self):
		self.filename=''
		self.object=''
		self.date_obs=''
		self.setting=''
		self.wave_min=0.
		self.wave_max=0.
		#self.order=0
		#self.merged=False

n_files = len(allfitsfiles)
full_list = [obsfile() for i in range(n_files)]

for count in range(len(allfitsfiles[0:10])): 
	full_list[count].filename = allfitsfiles[count]
	#print(datadir + full_list[count].filename)
	spec = edibles_spectrum(datadir + full_list[count].filename)
	#print(spec.header)
	full_list[count].object = spec.header["OBJECT"]
	full_list[count].date_obs = spec.header["DATE-OBS"]
	full_list[count].setting = spec.header["HIERARCH ESO INS GRAT2 WLEN"]
	wave,flux = spec.GetSpectrum()
	full_list[count].wave_min = np.min(wave)
	full_list[count].wave_max = np.max(wave)
	del spec

# Create arrays of formatted strings to print to a csv file now. 
pstrings = [['Filename', 'Object', 'DateObs', 'setting', 'WaveMin', 'WaveMax']]
for count in range(n_files):
	pstrings.append([full_list[count].filename, full_list[count].object, full_list[count].date_obs, \
			         full_list[count].setting, full_list[count].wave_min, full_list[count].wave_max])

# Time to print things out! Let's use csv format to do that. 
outfile = edibles_pythondir + "/data/" + datarelease + "_ObsLog.csv"
length_checker = np.vectorize(len)
all_lengths = length_checker(allfitsfiles)
print(np.max(all_lengths))
with open(outfile, 'w') as csvFile:
	writer=csv.writer(csvFile)
	writer.writerows(pstrings)
