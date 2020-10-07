import pandas as pd

def ReadSIMBAD():
    # File headers are: 
    # |typed ident| identifier |typ| coord1 (ICRS,J2000/2000) | plx | radvel |Mag U |Mag B |Mag V |Mag R |Mag I |Mag J |Mag H |Mag K |spec. type    
    Headers = ['Count', 'Object', 'Name', 'Object Type', 'Coordinates', 'Parallax', 'Radial Velocity', 'U', 'B', 'V', 'R', 'I', 'J', 'H', 'K', 'Spectral Type']
    file_name=str("SIMBAD_query_observed_Oct06_2020.ascii")
    val_data=pd.read_csv(file_name,sep='|',skiprows=7, skipfooter=3, skipinitialspace=True, na_values='~', header=None,names=Headers,engine='python')
    #print(val_data)
    return(val_data)

if __name__ == "__main__":
    # print("Main")
    data = ReadSIMBAD()
    print(data)