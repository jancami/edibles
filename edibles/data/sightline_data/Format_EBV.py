import math
import numpy as np
import pandas as pd
from edibles.data.sightline_data.ReadValencic_data import ReadValencic

"""
This routine will parse all the sources of E(B-V) data and process them to 
produce the standard formatted output tables required for the oracle. 
The E(B-V) values are listed in 2 different files, and represent 3 different
sources (in order of preference):
1. Valencic04: Ref 3 in ParameterReferences.csv
2. EBV_SIMBAD: Ref 4 in ParameterReferences.csv
3. Tycho:      Ref 2 in ParameterReferences.csv

For Valencic, the data is contained in "Valencic2004_datafile4.txt" and can be read in with
the ReadValencic fundtion. For 2. and 3, the data is contained in "InputEBV.csv" and can be 
turned into a pandas dataframe. Note that no uncertainties are listed for SIMBAD and Tycho data
so these are set to NaN. 

We will keep track of whether the preferred flag is set is a separate boolean array that we will
set to True as we process the different sources in order of preference. 

"""
ref_id_Valencic = 3
ref_id_SIMBAD = 4
ref_id_Tycho = 2


# Read in the list of observed objects. 
filename = str("ObservedObjects.csv")
df = pd.read_csv(filename)
observed_objects = df['Object'].tolist()
n_obj = len(observed_objects)
# We just need to keep track of the preferred flag per object, so initialize boolean array. 
preferred_flag_assigned = [False for i in range(n_obj)]

# creata a pandas dataframe to hold the standard formatted output data and a counter to keep track of
# how many total entries we have. 
df_out=pd.DataFrame(columns=["object","value","unc_lower","unc_upper","reference_id","preferred_flag"])
counter = 0

# Read in and process the E(B-V) data 
valencic_data = ReadValencic()
ebv_data = pd.read_csv('InputEBV.csv', names=['Object', 'EBV_Tycho', 'EBV_Simbad'])

# Now loop over each object, and go through all data sources. 
for objloop in range(n_obj):
    this_object = observed_objects[objloop]
    # Check Valencic first
    match = valencic_data.loc[valencic_data['Object'] == this_object]
    if len(match) != 0: 
        this_EBV = match['E(B-V)'].values[0]
        this_EBV_sigma = match['Sigma_E(B-V)'].values[0]
        # Append a row to the output dataframe and set the proper Reference ID. This should always
        # be the preferred value, so set the flag to 1 and mark that this flag is set for this object. 
        df_out.loc[counter] = (this_object, this_EBV, this_EBV_sigma, this_EBV_sigma, ref_id_Valencic, 1)
        counter += 1
        preferred_flag_assigned[objloop] = True
        
    # Then check for a match in the SIMBAD/Tycho data
    match = ebv_data.loc[ebv_data['Object'] == this_object]
    if len(match) != 0: 
        # Next down our list of preferred values is EBV_SIMBAD. However, some objects do not have this
        # value present, so check first and only proceed if we do have a value. 
        # Check whether we have a real value!! 
        this_EBV_SIMBAD = float(match['EBV_Simbad'].values[0])
        if not math.isnan(this_EBV_SIMBAD): 
            if (preferred_flag_assigned[objloop] == False):
                is_preferred=1
            else:
                is_preferred=0
            #print(observed_objects[objloop], preferred_flag_assigned[objloop], preferred_flag_assigned[objloop] == False, is_preferred)
            df_out.loc[counter] = (this_object, this_EBV_SIMBAD, np.NaN, np.NaN, ref_id_SIMBAD, is_preferred)
            counter += 1
            preferred_flag_assigned[objloop] = True
        # Then the Tycho values
        this_EBV_Tycho = float(match['EBV_Tycho'].values[0])
        if not math.isnan(this_EBV_Tycho): 
            if preferred_flag_assigned[objloop] == False:
                is_preferred=1 
            else:
                is_preferred=0
            df_out.loc[counter] = (this_object, this_EBV_Tycho, np.NaN, np.NaN, ref_id_Tycho, is_preferred)
            counter += 1
            preferred_flag_assigned[objloop] = True


# Sort properly on the HD number(i.e. strictly increasing) before printing to CSV. 
# We'll create the sort key (the number only) and then add it as an extra column
# to the data frame, sort it, and remove the column again. 
sortkey = df_out['object'].tolist()
sortkey = [int(x[3:]) for x in sortkey]
df_out_print = df_out.assign(S=sortkey).sort_values('S').drop('S', 1)
# Then finally print to the formatted csv file
df_out_print.to_csv('Formatted_EBV.csv', index=False, na_rep='NaN')
