import pandas as pd

def ReadJenkins():
    Headers=["Num","VComp","Object","Longitude","Latitude","Vmag", "SpType","SpType_ref","E(B-V)","Distance","Z","Lower_log(NHI)","log(NHI)","Flag_log(NHI)","Upper_log(NHI)","ref_log(NHI)","Lower_log(NH2)","log(NH2)","Flag_log(NH2)","Upper_log(NHI)","ref_log(NH2)"]
    
    file_name=str("Jenkins2009_data.txt")
    jen_data=pd.read_csv(file_name,sep='        |    '   ,skiprows=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18],header=None,names=Headers,engine='python')
    print(jen_data)
    return()
ReadJenkins()

