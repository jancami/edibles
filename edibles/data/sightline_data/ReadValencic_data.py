import pandas as pd
def CreateInputFile(Parameter):
    Headers=["Object","E(B-V)","Sigma_E(B-V)","R(V)","Sigma_R(V)","A(V)","Sigma_A(V)","Distance","MinDist","MaxDist"]
    file_name=str("Valencic2004_datafile4.txt")
    val_data=pd.read_csv(file_name,sep='        |    '   ,skiprows=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18],header=None,names=Headers,engine='python')
    if "Sigma_"+Parameter in Headers:
        df=val_data[["Object",Parameter,"Sigma_"+Parameter]].copy()
    else:
        df=val_data[["Object",Parameter]].copy()
    
        
        
    df.to_csv("Input"+Parameter+"_Val.csv",index=False)
for param in ("E(B-V)","R(V)","A(V)"):
    CreateInputFile(param)
