"""
small routine to reformat the multi-column sightline info to multi-row sightline info
"""

"""
remark: can we rename the target id to 'target' instead of 'object'. object is a thing in python...
But object is also a thing in SIMBAD. 

We will also set prefereed values here. 
For E(B-V): preferred = E(B-V)_Simbad

"""

import numpy as np
import glob
import pandas as pd

EBV_preferred = 'E(B-V)_Simbad'

# grab all the Input CSV files to convert
inputs=glob.glob('Input*.csv')

# initialise the reference list
ref=0
references=[]
references_counter=[]
ref_counter=1

for i in inputs:

    print(i) # prints the input file name
    
    data = np.genfromtxt(i, dtype=None, delimiter=",", names=True) # read the input CSV file

    colnames=data.dtype.names

    nr_rows = data.shape[0]

    # initialize the output target data
    object_id=[]
    value=[]
    unc_lower=[]
    unc_upper=[]
    reference_id=[]
    preferred_flag=[]

    # extract data from the input and put in the lists
    for row in np.arange(nr_rows):

        nr_cols = len(data[row])

        for col in np.arange(1,nr_cols):

           if col != 'False':

               object_id.append(data[row][0].astype(str))
               value.append(data[row][col].astype(str))
               unc_lower.append(np.nan)
               unc_upper.append(np.nan)

               reference_id.append(col+ref)
               preferred_flag.append(1)# not setting the PREFERED VALUE flag yet.

    parameter = i[5:-4]

    # creata a pandas dataframe representing the output data
    df=pd.DataFrame(list(zip(object_id,value,unc_lower,unc_upper,reference_id,preferred_flag)), columns=["object","value","unc_lower","unc_upper","reference_id","preferred_flag"])

    # write the pandas dataframe to CSV
    df.to_csv('Targets_'+parameter+'.csv', index=False, na_rep='NaN')

    # create the reference list
    for c in np.arange(1,nr_cols):
        references_counter.append(ref_counter)
        references.append(colnames[c])
        ref_counter=ref_counter+1
        print(ref_counter)

    ref=ref+nr_cols-1

# craete pandas dataframe for the reference and write to CSV
dfr=pd.DataFrame(list(zip(references_counter, references)), columns=["reference_id","source"])
dfr.to_csv('References.csv', index=False, na_rep='NaN')
