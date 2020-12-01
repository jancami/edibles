import math
import numpy as np
import pandas as pd
from edibles.data.sightline_data.ReadJenkins_data import ReadJenkins

"""
This routine will parse all the sources of Log(N(H2)) and Log(N(HI)) data and process them to produce the standard formatted output tables required for the oracle.
The Log(N(H2)) and Log(N(HI)) values are listed in 1 different file(s), and represent 1 different source(s) (in order of preference):
1. Jenkins09: Ref 7 in ParameterReferences.csv


For Jenkins, the data is contained in "Jenkins2009_data.txt" and can be read in with the ReadJenkins function.

We will keep track of whether the preferred flag is set is a separate boolean array that we will
set to True as we process the different sources in order of preference.

"""
def ErrorCalc(a,b,a_err,b_err):
    error=((10**a)*np.sqrt((5*(a**2)*(100**a)*(a_err**2)+4*(10**(a+b))*(a**2)*(a_err**2)+(100**b)*(a**2)*(a_err**2)+(b**2)*(100**b)*(b_err**2))/(2*(10**a)+(10**b)**2)))/(10*(10**a)+5*(10**b))
    return(error)
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

# Read in and process the Log(N(H2)) data
jenkins_data = ReadJenkins()

# Now loop over each object, and go through all data sources.
for objloop in range(n_obj):
    this_object = observed_objects[objloop]
    # Check Jenkins first
    match = jenkins_data.loc[jenkins_data['Object'] == this_object]
    if len(match) != 0:
        nh2_value = match['log(NH2)'].values[0]
        nhI_value = match['log(NHI)'].values[0]
        
        this_value=np.around((2*10**nh2_value)/((2*10**nh2_value)+(10**nhI_value)),decimals = 2)
        if math.isnan(this_value) or math.isnan(nh2_value) or math.isnan(nhI_value) :
            pass
        else:
            NH2_upper = match['Upper_log(NH2)'].values[0]
            NH2_unc_up= NH2_upper - nh2_value
            NHI_upper = match['Upper_log(NHI)'].values[0]
            NHI_unc_up= NHI_upper - nhI_value
            unc_up=np.around(ErrorCalc(nh2_value,nhI_value,NH2_unc_up,NHI_unc_up),decimals = 2)
            
            NH2_lower = match['Lower_log(NH2)'].values[0]
            NH2_unc_low= NH2_lower - nh2_value
            NHI_lower = match['Lower_log(NHI)'].values[0]
            NHI_unc_low= NHI_lower - nhI_value
            unc_low=np.around(ErrorCalc(nh2_value,nhI_value,NH2_unc_low,NHI_unc_low),decimals = 2)
            
        
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
df_out_print.to_csv('Formatted_f(H2).csv', index=False, na_rep='NaN')



