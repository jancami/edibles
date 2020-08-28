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
        filename=folder /"sightline_data"/"Targets_SpType.csv"
        self.sptypelog = pd.read_csv(filename)
        
        
        print(self.sptypelog.dtypes)
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


    def FilterEngine(self, object, log, value, unc_lower, unc_upper, reference_id):
        
        bool_object_matches = np.zeros(len(log.index),dtype=bool)
        if object is None:
             bool_object_matches = np.ones(len(log.index),dtype=bool)
        elif (isinstance(object, np.ndarray) | isinstance(object, list)):
                for thisobject in object:
                    bool_object_matches = (log.object == thisobject) | (bool_object_matches)
                    #print(bool_object_matches.sum())
        else: 
            print("EDIBLES Oracle is Panicking in getFilteredObsList: don't know what I'm dealing with!")
            
            
        # Initialize a boolean array to match all entries in the sightline file. 
        # Work through each of the criteria and add the corresponding filter criterion. 
        bool_value_matches = np.ones(len(log.index),dtype=bool)
        if value:
            # Only keep sightline if the value is an exact match. 
            bool_value_matches = log.value == value
        if unc_lower:
            bool_value_matches = (log.value > unc_lower) & bool_value_matches
        if unc_upper:
            bool_value_matches = (log.value < unc_upper) & bool_value_matches
        # Now process the references or "preferred" values. 
        # If reference is "All", we should not apply an additional filter. 
        # If reference is specified, filter on that reference. 
        # If no reference is specified, use the preferred value. 
        if reference_id is None:
            bool_value_matches = (log.preferred_flag == 1) & bool_value_matches
        elif reference_id=='All':
                pass
        else:
                #check if proper ref. is given [1,2] for EBV, [3,4] fpr SpT.
                bool_value_matches = (log.reference_id == reference_id) & bool_value_matches
        
        bool_combined_matches = bool_object_matches & bool_value_matches
        ind = np.where(bool_combined_matches)
        matching_objects = log.object.values[ind]

        print('getFilteredObslist: Found a total of ', bool_object_matches.sum(), ' object matches.')  
        print('getFilteredObslist: Found a total of ', bool_value_matches.sum(), ' E(B-V) matches.')  
        print('getFilteredObslist: Found a total of ', bool_combined_matches.sum(), ' combined matches.')  
        
        return matching_objects



    def getFilteredObsList(self,object=None, Wave=None, MergedOnly=False, OrdersOnly=False,EBV=None,EBV_min=None,EBV_max=None, EBV_reference=None, SpType=None, SpType_min=None, SpType_max=None, SpType_reference=None, WaveMin=None, WaveMax=None):
        
        #This method will provide a filtered list of observations that match 
        #the specified criteria on sightline/target parameters as well as
        #on observation criteria (e.g. wavelength range). 

        #Each sightline parameter is processed in the same way to produce boolean 
        #arrays to filter the sightlines. In a second step, we then search the 
        #obs log for matching objects within the additional criteria. 

        matching_objects_ebv = self.FilterEngine(object, self.ebvlog, EBV, EBV_min, EBV_max, EBV_reference)
        matching_objects_sptype = self.FilterEngine(object, self.sptypelog, SpType, SpType_min, SpType_max, SpType_reference)
        
        # set(lst1).intersection(lst2) to check 
        print(matching_objects_ebv)
        print(matching_objects_sptype)
        
        ####### in FilterEngine#############
        ####################################
        ####################################
        bool_object_matches = np.zeros(len(self.ebvlog.index),dtype=bool)
        if object is None:
             bool_object_matches = np.ones(len(self.ebvlog.index),dtype=bool)
        elif (isinstance(object, np.ndarray) | isinstance(object, list)):
                for thisobject in object:
                    bool_object_matches = (self.ebvlog.object == thisobject) | (bool_object_matches)
                    #print(bool_object_matches.sum())
        else: 
            print("EDIBLES Oracle is Panicking in getFilteredObsList: don't know what I'm dealing with!")
            
            
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
        
        ####################################
        ####### in FilterEngine#############
        ####################################
        
        
        print(matching_objects)
        
        # Now push this list of objects through for further filtering based on obs log
        FilteredObsList = self._getObsListFilteredByObsLogParameters(object=matching_objects, Wave=Wave, WaveMin=WaveMin, WaveMax=WaveMax, MergedOnly=MergedOnly, OrdersOnly=OrdersOnly)

        #print(FilteredObsList)

        '''
        bool_obslog_match=np.zeros(len(self.obslog.index),dtype=bool)
        obslog_objects=self.obslog.Object.values
        #print(type(obslog_objects))
        #print(type(self.ebvlog.iloc[ind].object))
        for i in ind:

            bool_obslog_match =  (self.obslog.Object.values == self.ebvlog.iloc[ind].object.values) or bool_obslog_match
        print(bool_obslog_match)

            bool_obslog_match =  (self.obslog.Object.values == self.ebvlog.iloc[ind].object.values) | bool_obslog_match
        inds = np.where(bool_obslog_match)
        
        matching_objects_obs = self.obslog.Object.values[inds]
        print(matching_objects_obs)

    
        matched_files=[]
        for i in range(len(matching_objects)):
            matched_files.append(self.getObsList(target=matching_objects[i],MergedOnly=True))
        print(matched_files)
        '''
        
        return (FilteredObsList)


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

    print("1. Results from getFilteredObsList: ")
    List=pythia.getFilteredObsList(MergedOnly=True,EBV_min=0.7,EBV_max=0.8, SpType='B0.5 III')    
    print(List)
    
#    print("2. Results from getFilteredObsList: ")
#    List=pythia.getFilteredObsList(MergedOnly=True,EBV=0.6,EBV_max=0.9)    
#    print(List)
#
#    print("3. Results from getFilteredObsList: ")
#    List=pythia.getFilteredObsList(object=['HD 145502', 'HD 149757'], MergedOnly=True, Wave=6614)
#    print(List)
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
