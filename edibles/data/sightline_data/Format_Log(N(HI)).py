import math
import numpy as np
import pandas as pd
from edibles.data.sightline_data.ReadJenkins_data import ReadJenkins

"""
This routine will parse all the sources of Log(N(HI)) data and process them to produce the standard formatted output tables required for the oracle.
The Log(N(HI)) values are listed in 1 different file(s), and represent 1 different source(s) (in order of preference):
1. Jenkins09: Ref 7 in ParameterReferences.csv


For Jenkins, the data is contained in "Jenkins2009_data.txt" and can be read in with the ReadJenkins function.

We will keep track of whether the preferred flag is set is a separate boolean array that we will
set to True as we process the different sources in order of preference.

"""
ref_id_Jenkins = 7


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

# Read in and process the Log(N(HI)) data
jenkins_data = ReadJenkins()

# Now loop over each object, and go through all data sources.
for objloop in range(n_obj):
    this_object = observed_objects[objloop]
    # Check Jenkins first
    match = jenkins_data.loc[jenkins_data['Object'] == this_object]
    if len(match) != 0:
        this_value = match['log(NHI)'].values[0]
        this_value_upper = match['Upper_log(NHI)'].values[0]
        unc_up= np.around(this_value_upper - this_value,decimals = 2)
        this_value_lower = match['Lower_log(NHI)'].values[0]
        unc_low= np.around(this_value - this_value_lower, decimals = 2)
        # Append a row to the output dataframe and set the proper Reference ID. This should always
        # be the preferred value, so set the flag to 1 and mark that this flag is set for this object.
        df_out.loc[counter] = (this_object, this_value, unc_low, unc_up, ref_id_Jenkins, 1)
        counter += 1
        preferred_flag_assigned[objloop] = True

# Sort properly on the HD number(i.e. strictly increasing) before printing to CSV.
# We'll create the sort key (the number only) and then add it as an extra column
# to the data frame, sort it, and remove the column again.
sortkey = df_out['object'].tolist()
sortkey = [int(x[3:]) for x in sortkey]
df_out_print = df_out.assign(S=sortkey).sort_values('S').drop('S', 1)
# Then finally print to the formatted csv file
df_out_print.to_csv('Formatted_LogN(HI).csv', index=False, na_rep='NaN')


