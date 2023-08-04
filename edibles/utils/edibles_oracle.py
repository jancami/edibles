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
    Args:
        verbose (int): what level of output do we want during runtime?
    """

    def __init__(self, verbose=False):
        #print(DATADIR)
        
        if verbose is None: 
            self.verbose=0 
        else: 
            self.verbose=verbose
        folder = Path(PYTHONDIR+"/data")
        filename=folder /"DR4_ObsLog.csv"
        self.obslog = pd.read_csv(filename)
        filename=folder /"sightline_data"/"Formatted_EBV.csv"
        self.ebvlog = pd.read_csv(filename)
        filename=folder /"sightline_data"/"Formatted_SpType.csv"
        self.sptypelog = pd.read_csv(filename)
        filename = folder /"sightline_data"/"Formatted_LogN(HI).csv"
        self.NHIlog = pd.read_csv(filename)
        filename = folder /"sightline_data"/"Formatted_LogN(H2).csv"
        self.NH2log = pd.read_csv(filename)
        filename = folder /"sightline_data"/"Formatted_f(H2).csv"
        self.fH2log = pd.read_csv(filename)
        filename = folder /"sightline_data"/"Formatted_RV.csv"
        self.rvlog = pd.read_csv(filename)
        filename = folder /"sightline_data"/"Formatted_AV.csv"
        self.avlog = pd.read_csv(filename)
        
        filename = folder /"sightline_data"/"ObservedObjects.csv"
        self.object_log = pd.read_csv(filename,names=["object"],header=0)
        
        
        #print(self.sptypelog.dtypes)
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
            bool_object_matches = np.ones(len(self.ebvlog.index),dtype=bool)
        elif (isinstance(object, np.ndarray) | isinstance(object, list)):
                for thisobject in object:
                    #print("Object in loop:", thisobject)
                    #print(self.obslog.Object == thisobject)
                    bool_object_matches = (self.obslog.Object == thisobject) | (bool_object_matches)
                    #print(bool_object_matches.sum())
        else: 
             bool_object_matches = self.ebvlog.object == object

        #print('Inside the function: number of matches is ', bool_object_matches.sum())


        # Do we have to filter out merged or single-order spectra? Note that if both
        # MergedOnly and OrdersOnly are True, only the Merged spectra will be returned.
        if MergedOnly and OrdersOnly:
            print("EDIBLES Oracle WARNING: ONLY RETURNING MERGED SPECTRA")

        bool_order_matches = self.obslog.Order != "Z"
        if OrdersOnly is True:
            bool_order_matches = self.obslog.Order != "ALL"
        if MergedOnly is True:
            bool_order_matches = self.obslog.Order == "ALL"

        #print(bool_order_matches)
        
        bool_wave_matches = np.ones(len(self.obslog.index),dtype=bool)
        if Wave: 
            bool_wave_matches = (self.obslog.WaveMin < Wave) & (self.obslog.WaveMax > Wave)
        if WaveMin: 
            bool_wave_matches = (self.obslog.WaveMax > WaveMin) & (bool_wave_matches)
        if WaveMax: 
            bool_wave_matches = (self.obslog.WaveMin < WaveMax) & (bool_wave_matches)

        ind = np.where(bool_object_matches & bool_order_matches & bool_wave_matches)
        
        if self.verbose != 0: 
            print(ind)
            print("**Filtered File List**")
            print(self.obslog.iloc[ind].Filename)
        return self.obslog.iloc[ind].Filename            


    def FilterEngine(self, object, log, value, unc_lower, unc_upper, reference_id):
        # Generic function to filter through the list of objects. 
        # Note: object should be a list or a numpy array type!

        # First, find all the objects in our log that match the specified objects.
        
        bool_object_matches = np.zeros(len(log.index),dtype=bool)
        if object is None:
             bool_object_matches = np.ones(len(log.index),dtype=bool)
        elif (isinstance(object, np.ndarray) | isinstance(object, list)):
                for thisobject in object:
                    bool_object_matches = (log.object == thisobject) | (bool_object_matches)
                    #print(bool_object_matches.sum())
        else: 
            print("EDIBLES Oracle is Panicking in FilterEngine: don't know what I'm dealing with!")
            
        # Next, find all the matches with the parameters -- but only if they are specified! 
        # Initialize a boolean array to match all entries in the sightline file. 
        # Then work through each of the criteria and add the corresponding filter criterion. 
        bool_value_matches = np.ones(len(log.index),dtype=bool)
        
        #print(bool_value_matches)
        if value is not None:
            # Only keep sightline if the value is an exact match. 
            bool_value_matches = (log.value == value)
        if unc_lower is not None:
            bool_value_matches = (log.value > unc_lower) & bool_value_matches
        if unc_upper is not None:
            print(value)
            print(unc_upper)
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
        #ind = np.where(bool_combined_matches)
        #matching_objects = log.object.values[ind]
        matching_objects_df = log.loc[bool_combined_matches, ['object','value']]

        if self.verbose: 
            print('EdiblesOracle:FilterEngine: Found a total of ', bool_object_matches.sum(), ' object match(es).')  
            print('EdiblesOracle:FilterEngine: Found a total of ', bool_value_matches.sum(), ' parameter match(es).')  
            print('EdiblesOracle:FilterEngine: Found a total of ', bool_combined_matches.sum(), ' combined match(es).')  
    
        
        return matching_objects_df

    def getFilteredObjects(self,object=None, Wave=None, \
                           EBV=None, EBV_min=None, EBV_max=None, EBV_reference=None, \
                           SpType=None, SpType_min=None, SpType_max=None, SpType_reference=None, \
                           WaveMin=None, WaveMax=None, LogNHI=None,LogNHI_min=None,LogNHI_max=None,\
                           LogNHI_reference=None,LogNH2=None,LogNH2_min=None,LogNH2_max=None, \
                           LogNH2_reference=None, fH2=None,fH2_min=None,fH2_max=None, \
                           fH2_reference=None, RV=None,RV_min=None,RV_max=None, \
                           RV_reference=None, AV=None,AV_min=None,AV_max=None, \
                           AV_reference=None):
        
        '''This method will provide a filtered list of objects that match 
        the specified criteria on sightline/target parameters as well as
        on observational criteria (e.g. wavelength range). This function consists
        of two steps: 

        | 1. Find all targets that match specified target parameters. This is done
           for each parameter using the FilterEngine function. 
        | 2. Find the objects that match all target specifications. '''
        # STEP 1: Filter objects for each of the parameters -- but only if parameters are specified!
        if (EBV or EBV_min or EBV_max or EBV_reference) is not None:
            print("EBV")
            matching_objects_ebv = self.FilterEngine(object, self.ebvlog, EBV, EBV_min, EBV_max, EBV_reference)
        else:
            matching_objects_ebv = self.object_log
        
        if (SpType or SpType_min or SpType_max or SpType_reference) is not None:
            print("SP_TYPE")
            matching_objects_sptype = self.FilterEngine(object, self.sptypelog, SpType, SpType_min, SpType_max, SpType_reference)
        else:
            matching_objects_sptype = self.object_log
            
        if (LogNHI or LogNHI_min or LogNHI_max or LogNHI_reference) is not None:
            print("LogN(HI)")
            matching_objects_logNHI = self.FilterEngine(object, self.NHIlog, LogNHI, LogNHI_min, LogNHI_max, LogNHI_reference)
        else:
            matching_objects_logNHI = self.object_log
        
        if (LogNH2 or LogNH2_min or LogNH2_max or LogNH2_reference) is not None:
            print("LogN(H2)")
            matching_objects_logNH2 = self.FilterEngine(object, self.NH2log, LogNH2, LogNH2_min, LogNH2_max, LogNH2_reference)
        else:
            matching_objects_logNH2 = self.object_log
        
        if (fH2 or fH2_min or fH2_max or fH2_reference) is not None:
            print("fH2")
            matching_objects_fH2 = self.FilterEngine(object, self.fH2log, fH2, fH2_min, fH2_max, fH2_reference)
        else:
            matching_objects_fH2 = self.object_log
            
        if (RV or RV_min or RV_max or RV_reference) is not None:
            print("RV")
            matching_objects_rv = self.FilterEngine(object, self.rvlog, RV, RV_min, RV_max, RV_reference)
        else:
            matching_objects_rv = self.object_log
        
        if (AV or AV_min or AV_max or AV_reference) is not None:
            print("AV")
            matching_objects_av = self.FilterEngine(object, self.avlog, AV, AV_min, AV_max, AV_reference)
        else:
            matching_objects_av = self.object_log
            
        
     
        # STEP 2: Find the common objects
        ebv_objects = matching_objects_ebv['object']
        sptype_objects = matching_objects_sptype['object']
        logNHI_objects = matching_objects_logNHI['object']
        logNH2_objects = matching_objects_logNH2['object']
        fH2_objects = matching_objects_fH2['object']
        rv_objects = matching_objects_rv['object']
        av_objects = matching_objects_av['object']
        
        #print(logNHI_objects.tolist())
        #print(ebv_objects.tolist())
        #print(sptype_objects.tolist())
        ##################
        if object is None:
            search_list = self.object_log["object"].to_list()
        else:
            search_list = object
        
        common_objects_set = set(search_list).intersection(set(ebv_objects.to_list()), set(sptype_objects.to_list()),
            set(logNHI_objects.to_list()), set(logNH2_objects.to_list()), set(fH2_objects.to_list()), set(rv_objects.to_list()), 
            set(av_objects.to_list()))
        
        ###################
        common_objects_list= list(common_objects_set)
        if self.verbose != 0: 
            print("EdiblesOracle:GetFilteredObjects:***Common Objects***")
            if len(common_objects_list) == 0:
                print("None")
            else:
                print(common_objects_list)

        return (common_objects_list)


    def getFilteredObsList(self,object=None, Wave=None, MergedOnly=False, OrdersOnly=False,\
                           EBV=None, EBV_min=None, EBV_max=None, EBV_reference=None, \
                           SpType=None, SpType_min=None, SpType_max=None, SpType_reference=None, \
                           WaveMin=None, WaveMax=None, LogNHI=None,LogNHI_min=None,LogNHI_max=None,\
                           LogNHI_reference=None,LogNH2=None,LogNH2_min=None,LogNH2_max=None, \
                           LogNH2_reference=None, fH2=None,fH2_min=None,fH2_max=None, \
                           fH2_reference=None, RV=None,RV_min=None,RV_max=None, \
                           RV_reference=None, AV=None,AV_min=None,AV_max=None, \
                           AV_reference=None):
        
        '''This method will provide a filtered list of observations that match 
        the specified criteria on sightline/target parameters as well as
        on observational criteria (e.g. wavelength range). This function consists
        of three steps: 

        | 1. Find all targets that match specified target parameters. This is done
           for each parameter using the FilterEngine function. 
        | 2. Find the objects that match all target specifications. 
        | 3. Find the observations that match specified parameters for only these targets. '''

        #print(getFilteredObslist.__dict__)

        # STEP 1: Filter objects for each of the parameters -- but only if parameters are specified!

        if (EBV or EBV_min or EBV_max or EBV_reference) is not None:
            if self.verbose != 0: 
                print("Filtering on E(B-V)...")
            matching_objects_ebv = self.FilterEngine(object, self.ebvlog, EBV, EBV_min, EBV_max, EBV_reference)
        else:
            matching_objects_ebv = self.object_log
        
        if (SpType or SpType_min or SpType_max or SpType_reference) is not None:
            if self.verbose != 0: 
                print("Filtering on Spectral Type...")
            matching_objects_sptype = self.FilterEngine(object, self.sptypelog, SpType, SpType_min, SpType_max, SpType_reference)
        else:
            matching_objects_sptype = self.object_log
            
        if (LogNHI or LogNHI_min or LogNHI_max or LogNHI_reference) is not None:
            if self.verbose != 0: 
                print("Filtering on N(HI)...")
            matching_objects_logNHI = self.FilterEngine(object, self.NHIlog, LogNHI, LogNHI_min, LogNHI_max, LogNHI_reference)
        else:
            matching_objects_logNHI = self.object_log
        
        if (LogNH2 or LogNH2_min or LogNH2_max or LogNH2_reference) is not None:
            if self.verbose != 0: 
                print("Filtering on N(H2)...")
            matching_objects_logNH2 = self.FilterEngine(object, self.NH2log, LogNH2, LogNH2_min, LogNH2_max, LogNH2_reference)
        else:
            matching_objects_logNH2 = self.object_log
        
        if (fH2 or fH2_min or fH2_max or fH2_reference) is not None:
            if self.verbose != 0: 
                print("Filtering on f(H2)...")
            matching_objects_fH2 = self.FilterEngine(object, self.fH2log, fH2, fH2_min, fH2_max, fH2_reference)
        else:
            matching_objects_fH2 = self.object_log
            
        if (RV or RV_min or RV_max or RV_reference) is not None:
            if self.verbose != 0: 
                print("Filtering on R(V)...")
            matching_objects_rv = self.FilterEngine(object, self.rvlog, RV, RV_min, RV_max, RV_reference)
        else:
            matching_objects_rv = self.object_log
        
        if (AV or AV_min or AV_max or AV_reference) is not None:
            if self.verbose != 0: 
                print("Filtering on A_V...")
            matching_objects_av = self.FilterEngine(object, self.avlog, AV, AV_min, AV_max, AV_reference)
        else:
            matching_objects_av = self.object_log
            
     
        # STEP 2: Find the common objects
        ebv_objects = matching_objects_ebv['object']
        sptype_objects = matching_objects_sptype['object']
        logNHI_objects = matching_objects_logNHI['object']
        logNH2_objects = matching_objects_logNH2['object']
        fH2_objects = matching_objects_fH2['object']
        rv_objects = matching_objects_rv['object']
        av_objects = matching_objects_av['object']
        
        #print(logNHI_objects.tolist())
        #print(ebv_objects.tolist())
        #print(sptype_objects.tolist())
        ##################
        if object is None:
            search_list = self.object_log["object"].to_list()
        else:
            search_list = object
        
        print("EdiblesOracle:GetFilteredObsList:Search List: ", search_list)
        print("EBV: ", ebv_objects.to_list())
        common_objects_set = set([search_list]).intersection(set(ebv_objects.to_list()))#,set(sptype_objects.to_list()))
            #set(logNHI_objects.to_list()), set(logNH2_objects.to_list()), set(fH2_objects.to_list()), set(rv_objects.to_list()),
            #set(av_objects.to_list()))

        print("Common: ", list(common_objects_set))
     
        ###################
        common_objects_list= list(common_objects_set)
        if self.verbose !=0: 
            print("EdiblesOracle:GetFilteredObsList:Common Objects")
            if len(common_objects_list) == 0:
                print("None")
            else:
                print(common_objects_list)
        
        # STEP 3
        # Now push this list of objects through for further filtering based on obs log
        FilteredObsList = self._getObsListFilteredByObsLogParameters(object=common_objects_list, Wave=Wave, WaveMin=WaveMin, WaveMax=WaveMax, MergedOnly=MergedOnly, OrdersOnly=OrdersOnly)

        #print(len(FilteredObsList))

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
            print("EDIBLES Oracle: ONLY RETURNING MERGED SPECTRA")

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
            print("EDIBLES Oracle: ONLY RETURNING MERGED SPECTRA")

        bool_order = self.obslog.Order != "Z"
        if OrdersOnly is True:
            bool_order = self.obslog.Order != "ALL"
        if MergedOnly is True:
            bool_order = self.obslog.Order == "ALL"

        ind = np.where(bool_target_matches & bool_order)
        return self.obslog.iloc[ind].Filename


if __name__ == "__main__":
    # print("Main")
    pythia = EdiblesOracle()

    # EXAMPLE 1: Get all observations for a single object. 
    List=pythia.getFilteredObsList(object=["HD 183143"], MergedOnly=True, Wave=3302.0)
    print("1. Results from getFilteredObsList: ")
    print(List)


    # EXAMPLE 2: Find all objects that match certain criteria. 
    List=pythia.getFilteredObjects(object=["HD 145502"], EBV_min=0.5, fH2_max=.3)
    print("2. Results from getFilteredObjects: ")
    print(List)

    #List=pythia.getFilteredObsList(object=["HD 103779"],MergedOnly=True,EBV_min=0.2,EBV_max=0.8,EBV_reference=3)
    #List=pythia.getFilteredObsList(EBV_min=0.2,EBV_max=0.8,EBV_reference=1)
    #List=pythia.getFilteredObsList(MergedOnly=True,EBV_min=0.2,EBV_max=0.8,EBV_reference=1)

    #print("1. Results from getFilteredObsList: ")
    #List=pythia.getFilteredObsList(MergedOnly=True,EBV_min=0.7,EBV_max=0.8, SpType='B0.5 III')    
    #List=pythia.getFilteredObsList(MergedOnly=True,EBV_min=0.2,EBV_max=0.8, object=['HD 145502'])
    ##List = pd.DataFrame(List).T
    ##List.columns = ['Object', 'EBV']
    ##print("Results from getFilteredObsList: ")
    ##print(List)

    #List=pythia.getFilteredObsList(object=['HD 145502', 'HD 149757'], MergedOnly=True, Wave=6614)
    ##List = pd.DataFrame(List).T
    ##List.columns = ['Object', 'EBV']
    ##print("Results from getFilteredObsList: ")
    ##print(List)
    
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
