
import glob
from SRC.AvgSpectrum import CreateAverageSpectrum
from SRC.Cont_Fit import Cont_Fit
from SRC.Data_Fit_Pathways import Data_Fit
Sightlines=['HD144470','HD147165','HD147683','HD149757','HD166937','HD170740','HD184915','HD185418','HD185859','HD203532','HD23180', 'HD24398']
DIBs=[6614,6379,5797,6196]


LMFIT=False
print("Continuum fitting each observation")
Cont_Fit(DIBs,Sightlines)
for DIB in DIBs:
    print(DIB)
    print('Creating average spectrums')
    for Sightline in Sightlines:

        files_list=glob.glob('Data/ContinuumNorm/'+str(DIB)+'/'+str(DIB)+'_'+Sightline+'*.csv')
      
        avg_spectrum=CreateAverageSpectrum(DIB,Sightline,files_list)
        
        if LMFIT==True:
            if DIB==6614:
                peak_type='triple'
            else:
                if DIB==6196:
                    print("6196 fitting not optimized. Will fail.")
                peak_type='double'

            Data=Data_Fit(DIB,Sightline,avg_spectrum[0],avg_spectrum[1],peak_type,FitReport=False,ShowPlot=False)


