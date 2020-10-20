import pandas as pd
import numpy as np
def ReadJenkins():
    Headers=["Num","VComp","Object","Longitude","Latitude","Vmag", "SpType","SpType_ref","E(B-V)","Distance","Z","Lower_log(NHI)","log(NHI)","Flag_log(NHI)","Upper_log(NHI)","ref_log(NHI)","Lower_log(NH2)","log(NH2)","Upper_log(NH2)","ref_log(NH2)"]
    
    rows_to_skip=np.linspace(0,73,74)
    colspecs=[(0,6),(10,13),(15,31),(32,38),(39,46),(46,51),(52,68),(69,75),(76,81),(82,87),(88,93),(94,99),(102,107),(108,109),(110,115),(118,124),(126,131),(132,137),(138,143),(144,150)]

    file_name=str("Jenkins2009_data.txt")
    jen_data=pd.read_fwf(file_name,skiprows=rows_to_skip,header=None,colspecs=colspecs,names=Headers,engine='python')
    return(jen_data)


