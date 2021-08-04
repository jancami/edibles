import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import lmfit
from lmfit.model import save_modelresult
from edibles.utils.simulations.SRC.Functions import Signal_Noise_Calculator, WavelengthToWavenumber
from scipy.interpolate import griddata
from edibles.utils.simulations.SimulatedContour import Simulated_Contour as sim
from timeit import default_timer as timer
import glob
from collections import OrderedDict
class Fitted_Contours:

    def __init__(self, DIB, Target, Initial_Params, J_Limit, Q_Branch=True):
        filename = 'SpectraData/'+str(DIB)+'/'+str(Target)+'_avg_spectra.csv'
        df = pd.read_csv(filename)
        self.DIB = DIB
        self.Target = Target
        self.wavelength = df['Wavelength'].to_numpy()
        self.flux = df['Flux'].to_numpy()
        lambda_guess=df["Wavelength"].iloc[df["Flux"].idxmin]
        Initial_Params.add('Lambda_0',lambda_guess,min=lambda_guess-0.45,max=lambda_guess+0.1)
        self.parms = Initial_Params
        self.J_Limit = J_Limit
        self.Q_Branch = Q_Branch
        
        filename="LocationsOfWings.csv"
        df2=pd.read_csv(filename)
        print(df2["Target"].dtype)
        if df2["Target"].str.contains(str(self.Target)).any():
            self.wing_flag=True
            wing_cutoff= df2["Wing-Coordinates"].loc[df2["Target"] == self.Target]
            self.wing_cutoff=wing_cutoff.iloc[0]
            first_cutoff=df2["Start-Coordinates"].loc[df2["Target"] == self.Target]
            self.first_cutoff=first_cutoff.iloc[0]
        else:
            self.wing_flag=False

    def calc_uncertainty(self):

        cont_x1  =  np.array(self.wavelength[-20:])
        cont_y1  =  np.array(self.flux[-20:])
        SN1,Fit1  =  Signal_Noise_Calculator(cont_x1,cont_y1)

        cont_x2  =  np.array(self.wavelength[:15])
        cont_y2  =  np.array(self.flux[:15])
        SN2,Fit2  =  Signal_Noise_Calculator(cont_x2,cont_y2)

        #Average of the two S/N
        Signal_Noise  =  0.5*(SN1+SN2)
        #Create an array of flux uncertainty values based on calculated S/N
        New_Noise  =  self.flux/Signal_Noise

        self.uncertainty_array  =  np.full(self.wavelength.shape,np.mean(New_Noise),dtype  =  float)
    
    


    def _residual(self, params, Data_X, Data_Y, Yerr, Jlimit=50, Symmetry_Type=None, Sightline=None, Q_Branch=False):
        

        # read in the dict of parameters. Apply appropriate constraints depending on symmetry type

        parvals = params.valuesdict()
        A = parvals['B']
        B = parvals['B']
        C = parvals['B']
        Q_scale = parvals['Q_Scale']
        PR_scale = parvals['PR_Scale']
        T = parvals['T']
        delta = parvals['Delta']
        lambda0 = parvals['Lambda_0']
        Y_Scale = parvals['Y_Scale']

        build = sim(A=A, B=B, C=C, Delta_A=delta*A, Delta_B=delta*B, Delta_C=delta*C, Trot=T, Jlimit=Jlimit,
                    Target=Sightline, lambda0=lambda0, Q_Branch=self.Q_Branch, Q_scale=Q_scale, PR_scale=PR_scale)
        old_x = np.asarray(build[0])

        old_y = np.asarray(build[1])

        new_y = griddata(old_x, old_y, Data_X, method='cubic', fill_value=1.0)
        new_y = (new_y/np.max(new_y*Y_Scale))*Y_Scale
        if self.wing_flag:
            arg=self._restrict_range(Data_X, new_y, Data_Y, Yerr)
            if self.uncertainty_array is None:
                residual = arg[0]-arg[1]
            else:
                residual = (arg[0]-arg[1])/arg[2]
        else:
            if self.uncertainty_array is None:
                residual = new_y-Data_Y
            else:
                residual = (new_y-Data_Y)/Yerr
            
        print(A, B, C, Q_scale, T, delta, lambda0)

#        plt.plot(Data_X,new_y,'c-')
#        plt.plot(Data_X,Data_Y,'b-')
#        plt.axvline(x=Data_X[arg[3]])
#        plt.axvline(x=Data_X[arg[4]])
#        plt.show()
        return(residual)
        
    def _restrict_range(self, x_values, model_y, data_y, data_error):
        #find index of first value greater than wing cutoff
        start_search=np.argmin(abs(x_values-self.first_cutoff))
        wing_search=np.argmin(abs(x_values-self.wing_cutoff))
        cut_new_y=model_y[start_search:wing_search]
        cut_model_y=data_y[start_search:wing_search]
        cut_error=data_error[start_search:wing_search]
        
        return(cut_new_y,cut_model_y, cut_error,wing_search, start_search)
        
        
    

        
    def minimize_contour(self,verbose=True):
        
        self.results=lmfit.minimize(self._residual,params=self.parms,kws={'Data_X': self.wavelength, 'Data_Y':self.flux,'Yerr':self.uncertainty_array, 'Jlimit': self.J_Limit, 'Sightline':self.Target ,'Q_Branch':self.Q_Branch}, method='nelder', nan_policy='omit',options={'disp':True,'fatol':0.001,'xatol':0.001})
        
        f = open('SavedModels/FinalFits2/'+str(self.DIB)+'_'+str(self.Target)+'_Lambda_'+str(round(self.results.params['Lambda_0'].value,2))+'FixedJ.txt', "w")
        f.write(str(self.results.params.valuesdict()))
        f.close()

        if verbose == True:
            print(self.Target)
            print(self.results.status)
            print(lmfit.fit_report(self.results))

    def plot_results(self):

        B = self.results.params['B'].value
        delta = self.results.params['Delta'].value
        T_r = self.results.params['T'].value
        Q_s = self.results.params['Q_Scale'].value
        Y_s = self.results.params['Y_Scale'].value
        P_s = self.results.params['PR_Scale'].value
        DIB_l = self.results.params['Lambda_0'].value

        build = sim(A=B, B=B, C=B, Delta_A=delta*B, Delta_B=delta*B, Delta_C=delta*B, Trot=T_r,
                    Jlimit=self.J_Limit, Q_scale=Q_s, PR_scale=P_s, Target=self.Target, lambda0=DIB_l, Q_Branch=True)

        new_y = griddata(build[0], build[1], self.wavelength, method='cubic', fill_value=1.0)
        new_y = (new_y/np.max(new_y*Y_s))*Y_s

        plt.title('Best Fit'+str(self.results.params))
        plt.plot(self.wavelength, self.flux, label='Original Data')
        plt.plot(self.wavelength, new_y, label='Best Fit')
        plt.legend()
        plt.savefig('SavedModels/FinalFits2/'+str(self.DIB)+'_'+str(self.Target)+'_Lambda_'+str(round(self.results.params['Lambda_0'].value,2))+'FixedJ.pdf')
        #plt.show()
        plt.close()
################################################################################


#####
def Fit_Contour(DIB, Target, Initial_Params, J_Limit, Q_Branch=True):
    build = Fitted_Contours(DIB, Target, Initial_Params, J_Limit, Q_Branch)
    build.calc_uncertainty()
    build.minimize_contour()
    build.plot_results()
    print('hi')


################################################################################

#Sightlines=['HD144470','HD147165','HD147683','HD149757','HD166937','HD170740','HD184915','HD185418','HD185859','HD203532','HD23180', 'HD24398']
Sightlines=['HD147165']
DIB=6614
################################################################################
##Create Initial Params##

p=lmfit.Parameters()
p.add('Q_Scale',0.3,min=0.0)
p.add('PR_Scale',0.6,min=0.0)
p.add('Y_Scale',0.005,min=0)
p.add('T', 25,min=3, max=50)
p.add('Delta', 0,min=-0.05, max=0.05)
p.add('B', 20e-3, min=10e-4, max=10e-1)

for Target in Sightlines:

#    file_name=glob.glob('SavedModels/FinalChoices/'+str(DIB)+'_'+str(Target)+'_Lambda_*.txt')[0]
#    file=open(file_name,'r')
#
#    contents=file.read()
#    saved_cont=dict(eval(contents))
#    B=saved_cont.get('B')
#    T=saved_cont.get('T')
#    Q_scale=saved_cont.get('Q_Scale')
#    PR_scale=saved_cont.get('PR_Scale')
#    delta=saved_cont.get('Delta')
#    Y_Scale=saved_cont.get('Y_Scale')
#    lambda0=saved_cont.get('Lambda_0')
    

    start=timer()
    print(start)
    print(Target)


    Fit_Contour(6614,Target,p,100)
    end=timer()
    print(end)
    print(end-start,' seconds')

