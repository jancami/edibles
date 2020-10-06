import numpy as np
import pandas as pd
from edibles.data.sightline_data.ReadValencic_data import ReadValencic

"""
This routine will parse all the sources of E(B-V) data and process them to 
produce the standard formatted output tables required for the oracle. 
The different sources are (in terms of preferences)
1. Valencic 2004  (Ref 3 in ParameterReferences.csv)
2. SIMBAD         (Ref 1 in ParameterReferences.csv)
3. Tycho          (Ref 2 in ParameterReferences.csv)

We will initialize the preferred flag to zero, and deal with the preferred flag at the end. 
"""

# Read in the list of observed objects. 
filename = str("ObservedObjects.csv")
df = pd.read_csv(filename)
observed_objects = df['Object'].tolist()
#nobj = len(observed_objects)

# creata a pandas dataframe to hold the formatted data. 
df_out=pd.DataFrame(columns=["object","value","unc_lower","unc_upper","reference_id","preferred_flag"])
counter = 0

# Get the Valencic data 
valdata = ReadValencic()
print(valdata.columns.values)

# Do the work. Go over each object and see if it's in the Valencic file. 
for this_object in observed_objects: 
	match = valdata.loc[valdata['Object'] == this_object]
	if len(match) != 0: 
	    this_EBV = match['E(B-V)'].values[0]
	    this_EBV_sigma = match['Sigma_E(B-V)'].values[0]
	    df_out.loc[counter] = (this_object, this_EBV, this_EBV_sigma, this_EBV_sigma, 1, 0)
	    counter += 1

df_out.to_csv('Formatted_EBV.csv', index=False, na_rep='NaN')
	
