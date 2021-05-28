import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from edibles.utils import edibles_oracle

def specie_charac_return(specie):
    '''
    Returns the characteristics of identified spiecies
    Specie must be : AlI,CaI,CaII,CrI,FeI,HeI*,KI,7LiI,6LiI,NaI,87RbI,TiII
    Return the wavelengths and f of the specie
    '''
    obslist=pd.read_csv(r'/Users/Samuel/Desktop/Programs/ImportGithub/edibles/data/auxiliary_data/line_catalogs/edibles_linelist_atoms.csv')
    
    wlength=[]
    flist=[]
    
    for i in range(len(obslist['Specie'])):
        s=obslist['Specie'][i]
        f=obslist['OscillatorStrength'][i]
        w=obslist['WavelengthAir'][i]
        if s==specie:
            wlength.append(w)
            flist.append(f)
            
    # print('Wavelenght for ',specie,':',wlength)
    # print('f for', specie,':',flist)
    
    return(wlength,flist)


class ISLineFitter():
    # def __init__(self, wave, flux):
    #     # wave and flux from edibles spectrum or elsewhere
    #     self.wave = wave
    #     self.flux = flux
    
    def __init__(self,target='', specie=''):
        'files from edibles or elsewhere'
        # self.target
        # self.specie
                

    def fit(self,target,specie):
        # this will do the fitting
        
        
        #We select all the wavelength of our target 
        WaveList=specie_charac_return(specie)[0]
        fList=specie_charac_return(specie)[1]
        
        #We make this loop for all the specie wavelenght, if we want a 
        #specified wavelenght we have to add an input
        for w in WaveList :
            
            pythia = edibles_oracle.EdiblesOracle()
            Slist = pythia.getObsListByWavelength(w, OrdersOnly=True)
        
            files = []
            for filename in Slist:
                if target in filename:
                    files.append(filename)
                    
        
        
        
        
        
        
        pass



    
    
    

        
        
        
        
        


if __name__ == "__main__":
    # print("Hello Word!")
    target='HD147889'
    specie='KI'
    
    
       
    isF=ISLineFitter()
    isF.fit(target, specie)
    
    print('Wavelenght for ',specie,':',specie_charac_return(specie)[0])
    print('f for', specie,':',specie_charac_return(specie)[1])
    
    
    
       