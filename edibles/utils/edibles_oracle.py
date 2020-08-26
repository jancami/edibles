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
        
    def getFilteredObsList(self,object=None, MergedOnly=False, OrdersOnly=False,EBV=None,EBV_min=None,EBV_max=None, EBV_reference=None):
        #load in proper files.
       
        #params: EBV, EBV_min,EBV_max, EBV_reference
        #if reference not specified, take preferrred value
        if object:
            bool_object_matches = (self.ebvlog.object == object)
        
        if EBV:
            bool_ebv_matches = (self.ebvlog.value < EBV_max) & (self.ebvlog.value > EBV_min)
    #reference flag not working yet
        '''
        #bool_prefer = self.ebvlog.preferred_flag != 100
        if EBV_reference:
            bool_prefer = self.ebvlog.preferred_flag == EBV_reference
        else:
            EBV_reference=1
            bool_prefer = self.ebvlog.preferred_flag == EBV_reference
        print(bool_prefer)
        '''
        '''
        bool_order = self.obslog.Order != "Z"
        if OrdersOnly is True:
            bool_order = self.obslog.Order != "ALL"
        if MergedOnly is True:
            bool_order = self.obslog.Order == "ALL"
        '''
      
        
        ind = np.where(bool_object_matches & bool_ebv_matches)
        # print(ind)
        return (self.ebvlog.iloc[ind].object,self.ebvlog.iloc[ind].value)
        
        
        
        #returns desired values.
        
        
        
        #creation of small data frame to be returned.
        
        
       #Returns list of filtered observations
        
        
        
        
        
        
        
        
        
        

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
    List=pythia.getFilteredObsList(object="b'HD 101065'",MergedOnly=True,EBV=0.8,EBV_min=0.7,EBV_max=1,EBV_reference=0)
    print(List)
    
    '''
    List = pythia.getObsListByWavelength(5000, MergedOnly=True)
    # print(List)
    for filename in List:
        sp = EdiblesSpectrum(filename)
        plt.figure()
        plt.title(filename)
        plt.xlabel("Wavelength (" + r"$\AA$" + ")")
        plt.xlim(5000, 5100)
        plt.plot(sp.wave, sp.flux)
        plt.show()
    '''
