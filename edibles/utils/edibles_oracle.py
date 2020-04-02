from edibles.edibles import DATADIR
from edibles.edibles import PYTHONDIR
from edibles.edibles.utils.edibles_spectrum import EdiblesSpectrum
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

class EdiblesOracle:
#	""
	#This class will pocess the EDIBLES obs log and target info files. 
	#Users can then query the oracle for observations matching specific criteria. 

	#""

	def __init__(self):
		print(DATADIR)
		filename = PYTHONDIR + '/edibles/data/DR4_ObsLog.csv'
		self.obslog = pd.read_csv(filename)
		#print(self.obslog.dtypes)
		total_rows = len(self.obslog.index)
		#print(total_rows)


	def GetObsListByWavelength(self, wave=None):
		# This function filters the list of Observations to return only those 
		# that include the requested wavelength. 

		if (wave is None):
			wave = 5000
		wave_matches = (self.obslog.WaveMin < wave) & (self.obslog.WaveMax > wave)
		ind = np.where(wave_matches)
		#print(ind)
		return self.obslog.iloc[ind].Filename


if __name__ == "__main__":
	#print("Main")
	pythia = EdiblesOracle()
	List = pythia.GetObsListByWavelength(5000)
	#print(List)
	for filename in List: 
		sp = EdiblesSpectrum(filename)
		plt.figure()
		plt.xlim(5000,5100)
		plt.plot(sp.wave, sp.flux)
		plt.show()




