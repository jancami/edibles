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
        
        
    def load_species_info(self,species=None,Wave=None, WaveMin=None, WaveMax=None, OscillatorStrength=None, OscillatorStrengthMin=None, OscillatorStrengthMax=None):
        '''This method will provide a filtered list of species information that matches
        the specified criteria on sightline/target parameters as well as
        on observational criteria (e.g. wavelength range).
    
         '''
        folder = Path(PYTHONDIR+"/data")
        filename = folder / "auxiliary_data/line_catalogs/edibles_linelist_atoms.csv"
        self.species_df=pd.read_csv(filename)
        #print(self.species_df)
        #search for matches by species name
        #print('Inside the function: object is', species)

        bool_species_matches = np.zeros(len(self.species_df.index),dtype=bool)
        
        if species is None:
            bool_species_matches = np.ones(len(self.species_df.index),dtype=bool)
        elif (isinstance(species, np.ndarray) | isinstance(species, list)):
            
            for thisobject in species:

                bool_species_matches = (self.species_df.Species.str.contains(thisobject)) | (bool_species_matches)
                
        else:
            
            bool_species_matches = self.ebvlog.object == object

        bool_wave_matches = np.ones(len(self.species_df.index),dtype=bool)
        if Wave:
            bool_wave_matches = (self.species_df.WavelengthAir == Wave)
        if WaveMin:
            bool_wave_matches = (self.species_df.WavelengthAir > WaveMin) & (bool_wave_matches)
        if WaveMax:
            bool_wave_matches = (self.species_df.WavelengthAir < WaveMax) & (bool_wave_matches)
            
        
        bool_osc_matches = np.ones(len(self.species_df.index),dtype=bool)
        if OscillatorStrength:
            bool_osc_matches = (self.species_df.OscillatorStrength == OscillatorStrength)
        if OscillatorStrengthMin:
            bool_osc_matches = (self.species_df.OscillatorStrength > OscillatorStrengthMin) & (bool_osc_matches)
        if OscillatorStrengthMax:
            bool_osc_matches = (self.species_df.OscillatorStrength < OscillatorStrengthMax) & (bool_osc_matches)
            
            
            
        ind = np.where(bool_species_matches & bool_wave_matches & bool_osc_matches)[0]
        self.species_list=self.species_df['Species'].iloc[ind].to_numpy()
        self.air_wavelength=self.species_df['WavelengthAir'].iloc[ind].to_numpy()
        self.oscillator_strength=self.species_df['OscillatorStrength'].iloc[ind].to_numpy()
        return(self.species_list,self.air_wavelength,self.oscillator_strength)

if __name__ == "__main__":
    print("Hello Word!")
    ##Random values for flux and wave to init. Remove before public push, or can be changed by others.
    import random
    wave=np.linspace(0,10,11)
    flux=np.asarray((random.sample(range(100), k=len(wave))))/100
    ####################################################
    fit=ISLineFitter(wave,flux)
    test_species_info=fit.load_species_info(species=['Na'],OscillatorStrengthMin=0.1)
    print(test_species_info)
