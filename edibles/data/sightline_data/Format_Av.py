import math
import numpy as np
import pandas as pd
from edibles.data.sightline_data.ReadValencic_data import ReadValencic

"""
This routine will parse all the sources of A(V) data and process them to produce the standard formatted output tables required for the oracle.
The A(V) values are listed in 1 different file(s), and represent 1 different source(s) (in order of preference):
1. Valencic04: Ref 3 in ParameterReferences.csv


For Valencic, the data is contained in "Valencic2004_datafile4.txt" and can be read in with the ReadValencic function.

We will keep track of whether the preferred flag is set is a separate boolean array that we will
set to True as we process the different sources in order of preference.

"""
ref_id_Valencic = 3


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

# Now loop over each object, and go through all data sources.
for objloop in range(n_obj):
    this_object = observed_objects[objloop]
    # Check Valencic first
    match = valencic_data.loc[valencic_data['Object'] == this_object]
    if len(match) != 0:
        this_EBV = match['A(V)'].values[0]
        this_EBV_sigma = match['Sigma_A(V)'].values[0]
        # Append a row to the output dataframe and set the proper Reference ID. This should always
        # be the preferred value, so set the flag to 1 and mark that this flag is set for this object.
        df_out.loc[counter] = (this_object, this_EBV, this_EBV_sigma, this_EBV_sigma, ref_id_Valencic, 1)
        counter += 1
        preferred_flag_assigned[objloop] = True

# Sort properly on the HD number(i.e. strictly increasing) before printing to CSV.
# We'll create the sort key (the number only) and then add it as an extra column
# to the data frame, sort it, and remove the column again.
sortkey = df_out['object'].tolist()
sortkey = [int(x[3:]) for x in sortkey]
df_out_print = df_out.assign(S=sortkey).sort_values('S').drop('S', 1)
# Then finally print to the formatted csv file
df_out_print.to_csv('Formatted_AV.csv', index=False, na_rep='NaN')

