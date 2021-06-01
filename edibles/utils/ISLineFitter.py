import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from edibles import DATADIR
from edibles import PYTHONDIR
class ISLineFitter():
    def __init__(self, wave, flux):
        # wave and flux from edibles spectrum or elsewhere
        self.wave = wave
        self.flux = flux

    def fit(self):
        # this will do the fitting
        pass
        
        
    def load_species_info(self,species=None):
        
        folder = Path(PYTHONDIR+"/data")
        filename = folder / "auxiliary_data/line_catalogs/edibles_linelist_atoms.csv"
        species_df=pd.read_csv(filename)
        #search for matches by species name
#        bool_object_matches = np.zeros(len(species.index),dtype=bool)
#        if species is None:
#            bool_object_matches = np.ones(len(self.ebvlog.index),dtype=bool)
        
        
        

        

        search_name=species_df[species_df['Species'].str.contains(species)]
        if len(search_name)==0:
            print('Do not have information about species '+species)
            return
        else:
            print(search_name)
            self.air_wavelength=search_name['WavelengthAir'].to_numpy()
            self.oscillator_strength=search_name['OscillatorStrength'].to_numpy()


if __name__ == "__main__":
    print("Hello Word!")
    ##Random values for flux and wave to init. Remove before push
    import random
    wave=np.linspace(0,10,11)
    flux=np.asarray((random.sample(range(100), k=len(wave))))/100
    ####################################################
    fit=ISLineFitter(wave,flux)
    a_test=fit.load_species_info('NaP')
    b_test=fit.load_species_info('Na')
    print(fit.air_wavelength)
