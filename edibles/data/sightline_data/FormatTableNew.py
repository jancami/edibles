import pandas as pd
import numpy as np
import glob

def FormatAv():
    ##Find All Input CSV files that contain A(v)
    files=glob.glob("InputA(V)*.csv")
    print(files)

FormatAv()
    
