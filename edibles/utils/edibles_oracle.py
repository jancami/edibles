import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path
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
        
        folder = Path(PYTHONDIR+"/data")
        filename=folder /"DR4_ObsLog.csv"
        self.obslog = pd.read_csv(filename)
        filename=folder /"sightline_data"/"Targets_EBV.csv"
        self.ebvlog = pd.read_csv(filename)
        
        
        print(self.ebvlog.dtypes)
        # total_rows = len(self.ebvlog.index)
        # print(total_rows)

    def _getObsListFilteredByObsLogParameters(self, object=None, Wave=None, WaveMin=None, WaveMax=None, MergedOnly=False, OrdersOnly=False):
        '''Filter all the observations in the ObsLog by the parameters
        contained in the obslog, i.e. by object (if specified), wavelength
        range or merged versus specific orders. '''
        
        # We will use Boolean matches for all filter criteria. 
        
        #print('Inside the function: object is', object)

        bool_object_matches = np.zeros(len(self.obslog.index),dtype=bool)
        #print(object.dtype)
        if object is None:
            pass
        elif (isinstance(object, np.ndarray) | isinstance(object, list)):
                for thisobject in object:
                    #print("Object in loop:", thisobject)
                    #print(self.obslog.Object == thisobject)
                    bool_object_matches = (self.obslog.Object == thisobject) | (bool_object_matches)
                    #print(bool_object_matches.sum())
        else: 
            print("EDIBLES Oracle is Panicking in getObsListFilteredByObsLogParameters: don't know what I'm dealing with!")

        #print('Inside the function: number of matches is ', bool_object_matches.sum())


        # Do we have to filter out merged or single-order spectra? Note that if both
        # MergedOnly and OrdersOnly are True, only the Merged spectra will be returned.
        if MergedOnly and OrdersOnly:
            print("WARNING: ONLY RETURNING MERGED SPECTRA")

        bool_order_matches = self.obslog.Order != "Z"
        if OrdersOnly is True:
            bool_order_matches = self.obslog.Order != "ALL"
        if MergedOnly is True:
<<<<<<< HEAD
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
            if isinstance(object, list):
                bool_object_matches = (np.arry(self.ebvlog.object) == object)
            else:
                bool_object_matches = (self.ebvlog.object == object)
        else:
            bool_object_matches = np.ones(len(self.ebvlog.index),dtype=bool)
=======
            bool_order_matches = self.obslog.Order == "ALL"

        # Boolean matches for wavelength.
        # If Wave is specified, we will interpret the query as requiring that particular 
        # wavelength to be included in the range covered by the observation. 
        # If WaveMin is specified, we include observations that have wavelengths
        # Above WaveMin. Similar for WaveMax
        bool_wave_matches = np.ones(len(self.obslog.index),dtype=bool)
        if Wave: 
            bool_wave_matches = (self.obslog.WaveMin < Wave) & (self.obslog.WaveMax > Wave)
        if WaveMin: 
            bool_wave_matches = (self.obslog.WaveMax > WaveMin) & (bool_wave_matches)
        if WaveMax: 
            bool_wave_matches = (self.obslog.WaveMin < WaveMax) & (bool_wave_matches)


        # Now combine all filters
        bool_filter = bool_object_matches & bool_order_matches & bool_wave_matches
        print('Total filter results: ', bool_filter.sum())
        ind = np.where(bool_filter)

        #print(ind)
        #print(' result', self.obslog.iloc[ind].Filename)
        return self.obslog.Filename.values[ind]        


    def getFilteredObsList(self,object=None, Wave=None, MergedOnly=False, OrdersOnly=False,EBV=None,EBV_min=None,EBV_max=None, EBV_reference=None,WaveMin=None, WaveMax=None):
        
        #This method will provide a filtered list of observations that match 
        #the specified criteria on sightline/target parameters as well as
        #on observation criteria (e.g. wavelength range). 

        #Each sightline parameter is processed in the same way to produce boolean 
        #arrays to filter the sightlines. In a second step, we then search the 
        #obs log for matching objects within the additional criteria. 

        bool_object_matches = np.zeros(len(self.ebvlog.index),dtype=bool)
        if object is None:
             pass
        elif (isinstance(object, np.ndarray) | isinstance(object, list)):
                for thisobject in object:
                    bool_object_matches = (self.ebvlog.object == thisobject) | (bool_object_matches)
                    #print(bool_object_matches.sum())
        else: 
            print("EDIBLES Oracle is Panicking in getFilteredObsList: don't know what I'm dealing with!")
>>>>>>> 2d49afc7e47c345f54e571f3cbbf99e94df682e7
            
    
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
        # Now process the references or "preferred" values. 
        # If reference is "All", we should not apply an additional filter. 
        # If reference is specified, filter on that reference. 
        # If no reference is specified, use the preferred value. 
        if EBV_reference is None:
            bool_ebv_matches = (self.ebvlog.preferred_flag == 1) & bool_ebv_matches
        elif EBV_reference=='All':
                pass
        else:
                #check if proper ref. is given [1,2] for EBV, [3,4] fpr SpT.
                bool_ebv_matches = (self.ebvlog.reference_id == EBV_reference) & bool_ebv_matches
        
        bool_combined_matches = bool_object_matches & bool_ebv_matches
        ind = np.where(bool_combined_matches)
        matching_objects = self.ebvlog.object.values[ind]

        print('getFilteredObslist: Found a total of ', bool_object_matches.sum(), ' object matches.')  
        print('getFilteredObslist: Found a total of ', bool_ebv_matches.sum(), ' E(B-V) matches.')  
        print('getFilteredObslist: Found a total of ', bool_combined_matches.sum(), ' combined matches.')  
        
        print(matching_objects)
        
        # Now push this list of objects through for further filtering based on obs log
        FilteredObsList = self._getObsListFilteredByObsLogParameters(object=matching_objects, Wave=Wave, WaveMin=WaveMin, WaveMax=WaveMax, MergedOnly=MergedOnly, OrdersOnly=OrdersOnly)

        #print(FilteredObsList)

<<<<<<< HEAD
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
        matching_objects = self.ebvlog.object.values[ind]
        matching_ebv = self.ebvlog.value.values[ind]
        self.matching_objects = self.ebvlog.object.values[ind]
      
        # Now push this list through for further filtering based on obs log
=======
>>>>>>> 2d49afc7e47c345f54e571f3cbbf99e94df682e7
        '''
        bool_obslog_match=np.zeros(len(self.obslog.index),dtype=bool)
        obslog_objects=self.obslog.Object.values
        #print(type(obslog_objects))
        #print(type(self.ebvlog.iloc[ind].object))
        for i in ind:
<<<<<<< HEAD
            bool_obslog_match =  (self.obslog.Object.values == self.ebvlog.iloc[ind].object.values) or bool_obslog_match
        print(bool_obslog_match)
=======

            bool_obslog_match =  (self.obslog.Object.values == self.ebvlog.iloc[ind].object.values) or bool_obslog_match
        print(bool_obslog_match)

>>>>>>> 2d49afc7e47c345f54e571f3cbbf99e94df682e7
            bool_obslog_match =  (self.obslog.Object.values == self.ebvlog.iloc[ind].object.values) | bool_obslog_match
        inds = np.where(bool_obslog_match)
        
        matching_objects_obs = self.obslog.Object.values[inds]
        print(matching_objects_obs)
<<<<<<< HEAD
        '''

        '''
=======

    
>>>>>>> 2d49afc7e47c345f54e571f3cbbf99e94df682e7
        matched_files=[]
        for i in range(len(matching_objects)):
            matched_files.append(self.getObsList(target=matching_objects[i],MergedOnly=True))
        print(matched_files)
        '''
        
<<<<<<< HEAD
        return (matching_objects,matching_ebv)
        
        
  
        
        
        self.getObsList(WaveMin=None, WaveMax=None, MergedOnly=False, OrdersOnly=False)
        return (self.matching_objects)
        
        
=======
        return (FilteredObsList)
>>>>>>> 2d49afc7e47c345f54e571f3cbbf99e94df682e7


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
<<<<<<< HEAD
    List=pythia.getFilteredObsList(object="HD 103779", MergedOnly=True,EBV_min=0.2,EBV_max=0.8,EBV_reference=1)
    List=pythia.getFilteredObsList(EBV_min=0.2,EBV_max=0.8,EBV_reference=1)
    List=pythia.getFilteredObsList(MergedOnly=True,EBV_min=0.2,EBV_max=0.8,EBV_reference=1)
=======

    List=pythia.getFilteredObsList(MergedOnly=True,EBV_min=0.2,EBV_max=0.8, object='HD 145502')
>>>>>>> 2d49afc7e47c345f54e571f3cbbf99e94df682e7
    print("Results from getFilteredObsList: ")
    print(List)

    List=pythia.getFilteredObsList(object=['HD 145502', 'HD 149757'], MergedOnly=True, Wave=6614)
    print("Results from getFilteredObsList: ")
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
