import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from edibles import DATADIR
from edibles import PYTHONDIR
from edibles.utils.edibles_spectrum import EdiblesSpectrum


class EdiblesOracle:
    """
    This class will process the EDIBLES obs log and target info files.
    Users can then query the oracle for observations matching specific criteria.
    """

    def __init__(self):
        print(DATADIR)
        filename = PYTHONDIR + "/data/DR4_ObsLog.csv"
        self.obslog = pd.read_csv(filename)
        filename = PYTHONDIR + "/data/sightline_data/Targets_EBV.csv"
        self.ebvlog = pd.read_csv(filename)
        
        
        print(self.ebvlog.dtypes)
        # total_rows = len(self.ebvlog.index)
        # print(total_rows)

    def getObsList(self, WaveMin=None, WaveMax=None, MergedOnly=False, OrdersOnly=False):
        '''Combine wavelength and object search
        first select files based on the object lists produced by
        getFilteredObsList and then select the correct wavelength regime'''
        
        # Boolean matches for wavelength.
        
        bool_object_matches = np.zeros(len(self.obslog.index),dtype=bool)
        for obj in self.matching_objects:
            bool_object_matches = (self.obslog.Object == obj) or bool_object_matches
         
        # Do we have to filter out merged or single-order spectra? Note that if both
        # MergedOnly and OrdersOnly are True, only the Merged spectra will be returned.

        if MergedOnly and OrdersOnly:
            print("ONLY RETURNING MERGED SPECTRA")

        bool_order = self.obslog.Order != "Z"
        if OrdersOnly is True:
            bool_order = self.obslog.Order != "ALL"
        if MergedOnly is True:
            bool_order = self.obslog.Order == "ALL"

        print(bool_order)
        ind = np.where(bool_object_matches & bool_order)

        print(ind)
        print(' result', self.obslog.iloc[ind].Filename)
        return self.obslog.iloc[ind].Filename        



    def getFilteredObsList(self,object=None, MergedOnly=False, OrdersOnly=False,EBV=None,EBV_min=None,EBV_max=None, EBV_reference=None,WaveMin=None, WaveMax=None):
     	
    	#This method will provide a filtered list of observations that match 
    	#the specified criteria on sightline/target parameters as well as
    	#on observation criteria (e.g. wavelength range). 

    	#Each sightline parameter is processed in the same way to produce boolean 
    	#arrays to filter the sightlines. In a second step, we then search the 
    	#obs log for matching objects within the additional criteria. 


        # Note that the below statement will only work if there is a single object specified. 
        # We probably need to consider an array of objects as well. 
        if object:
            bool_object_matches = (self.ebvlog.object == object)
        else:
            bool_object_matches = np.ones(len(self.ebvlog.index),dtype=bool)
            
            
        # Initialize a boolean array to match all entries in the sightline file. 
        # Work through each of the criteria and add the corresponding filter criterion. 
        bool_ebv_matches = np.ones(len(self.ebvlog.index),dtype=bool)
        if EBV:
            # Only keep sightline if the value is an exact match. 
            bool_ebv_matches = self.ebvlog.value == EBV
        if EBV_min:
            bool_ebv_matches = (self.ebvlog.value > EBV_min) & bool_ebv_matches
        if EBV_max:
            bool_ebv_matches = (self.ebvlog.value < EBV_max) & bool_ebv_matches
        if EBV_reference:
            # If reference is "All", we should not apply an additional filter. 
            # If reference is specified, filter on that reference. 
            # If no reference is specified, use the preferred value. 
            if EBV_reference=='All':
                pass
            else:
                #check if proper ref. is given [1,2] for EBV, [3,4] fpr SpT.
                bool_ebv_matches = (self.ebvlog.reference_id == EBV_reference) & bool_ebv_matches
        else:
            bool_ebv_matches = (self.ebvlog.preferred_flag == 1) & bool_ebv_matches
        


        '''
        bool_order = self.obslog.Order != "Z"
        if OrdersOnly is True:
            bool_order = self.obslog.Order != "ALL"
        if MergedOnly is True:
            bool_order = self.obslog.Order == "ALL"
        '''
      
        #print(bool_object_matches)  
        #print(bool_ebv_matches)  
        ind = np.where(bool_object_matches & bool_ebv_matches)
        self.matching_objects = self.ebvlog.object.values[ind]
      
        # Now push this list through for further filtering based on obs log

        '''
        bool_obslog_match=np.zeros(len(self.obslog.index),dtype=bool)
        obslog_objects=self.obslog.Object.values
        print(type(obslog_objects))
        print(type(self.ebvlog.iloc[ind].object))
        for i in ind:
            bool_obslog_match =  (self.obslog.Object.values == self.ebvlog.iloc[ind].object) or bool_obslog_match
        print(bool_obslog_match)
        '''
        self.getObsList(WaveMin=None, WaveMax=None, MergedOnly=False, OrdersOnly=False)
        return (self.matching_objects)
        
        


    def getObsListByWavelength(self, wave=None, MergedOnly=False, OrdersOnly=False):
        """
        This function filters the list of Observations to return only those
        that include the requested wavelength.
        We will create a set of boolean arrays that we will then combined
        as the filter.

        :param wave: Wavelength that the returned files will include
        :type wave: float
        :param MergedOnly: Only include spectra from merged orders
        :type MergedOnly: bool
        :param OrdersOnly: Only include individual spectrum orders
        :type OrdersOnly: bool

        """

        # Boolean matches for wavelength.
        if wave is None:
            wave = 5000
        bool_wave_matches = (self.obslog.WaveMin < wave) & (self.obslog.WaveMax > wave)

        # Do we have to filter out merged or single-order spectra? Note that if both
        # MergedOnly and OrdersOnly are True, only the Merged spectra will be returned.

        if MergedOnly and OrdersOnly:
            print("ONLY RETURNING MERGED SPECTRA")

        bool_order = self.obslog.Order != "Z"
        if OrdersOnly is True:
            bool_order = self.obslog.Order != "ALL"
        if MergedOnly is True:
            bool_order = self.obslog.Order == "ALL"

        ind = np.where(bool_wave_matches & bool_order)
        # print(ind)
        return self.obslog.iloc[ind].Filename
        

    def getObsListByTarget(self, target=None, MergedOnly=False, OrdersOnly=False):
    
        """
        This function filters the list of Observations to return only those
        of the requested target.
        We will create a set of boolean arrays that we will then combined
        as the filter.

        :param target: Target name that the returned files will include
        :type target: object
        :param MergedOnly: Only include spectra from merged orders
        :type MergedOnly: bool
        :param OrdersOnly: Only include individual spectrum orders
        :type OrdersOnly: bool

        """
        
        # Boolean matches for wavelength.
        if target is None:
            target = 'HD164073'
        bool_target_matches = (self.obslog.Object == target)
        
        # Do we have to filter out merged or single-order spectra? Note that if both
        # MergedOnly and OrdersOnly are True, only the Merged spectra will be returned.

        if MergedOnly and OrdersOnly:
            print("ONLY RETURNING MERGED SPECTRA")

        bool_order = self.obslog.Order != "Z"
        if OrdersOnly is True:
            bool_order = self.obslog.Order != "ALL"
        if MergedOnly is True:
            bool_order = self.obslog.Order == "ALL"

        ind = np.where(bool_target_matches & bool_order)
        # print(ind)
        return self.obslog.iloc[ind].Filename


if __name__ == "__main__":
    # print("Main")
    pythia = EdiblesOracle()
    List=pythia.getFilteredObsList(MergedOnly=True,EBV_min=0.2,EBV_max=0.8,EBV_reference=1)
    print("Results from getFilteredObsList: ")
    print(List)
    
    

    List = pythia.getObsListByWavelength(5000, MergedOnly=True)
    print(List)
    '''
    for filename in List:
        sp = EdiblesSpectrum(filename)
        plt.figure()
        plt.title(filename)
        plt.xlabel("Wavelength (" + r"$\AA$" + ")")
        plt.xlim(5000, 5100)
        plt.plot(sp.wave, sp.flux)
        plt.show()
    '''
