import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy import constants as const
from astropy import units as u
from astropy.convolution import Gaussian1DKernel, convolve
from astropy.modeling.models import Voigt1D
pd.options.mode.chained_assignment = None
c=const.c.to('km/s')
kms=(u.km / u.s)
def WavelengthToWavenumber(values):
    wavenumbers=1/(values*(10**-8))
    return(wavenumbers)
class Rotational_Energies:

    def __init__(self,A,B,C, Name,Q_scale,PR_scale):
        #print('init')
        self.A=A
        self.B=B
        self.C=C
        self._determine_symmetry_type()
        self.name=Name
        self.Q_scale=Q_scale
        self.PR_scale=PR_scale
        
    def _determine_symmetry_type(self):
        self.flag=False
        if self.A==self.B==self.C:
            self.symmetry_type='spherical'
        elif self.B==self.C and self.A>self.B:
            self.symmetry_type='symmetric_prolate'
        elif self.B==self.A and self.B>self.C:
            self.symmetry_type='symmetric_oblate'
        else:
            self.symmetry_type='asymmetric'
            print('Asymmetric tops are too hard for me to deal with right now!')
            self.flag=True
        print('Symmetry Type: ',self.symmetry_type)
    
    def rotational_energies(self,Jlimit):
        #print("Calculating rotational energy values")
        self.Jlimit=Jlimit
        if self.symmetry_type=='spherical':
            df=pd.DataFrame(columns=["J","E"])
            J_vals=list(range(0,self.Jlimit+1))
            
            for J in J_vals:
                
                E=self.B*J*(J+1)
                df=df.append({"J":J, "E":E},ignore_index=True)
                self.K=pd.Series(np.nan)
        else:
            df=pd.DataFrame(columns=["J","K"])
            J_vals=list(range(0,self.Jlimit+1))
            
            for J in J_vals:
                for K in list(range(-J,J+1)):
                    
                    df=df.append({"J" : J, "K": K}, ignore_index=True)
            if self.symmetry_type=='symmetric_prolate':
                df["E"]=self.B*df["J"]*(df["J"]+1)+(self.A-self.B)*df["K"]**2
            elif self.symmetry_type=='symmetric_oblate':
                df["E"]=self.B*df["J"]*(df["J"]+1)+(self.C-self.B)*df["K"]**2
            self.K=df["K"]
        self.J=df["J"]
        
        self.E=df["E"]
        
    def boltzmann(self,T):
        #print("Calculating state population")
        h=const.h.value
        c=const.c.to('cm/s').value
        k=const.k_B.value
        self.T=T
        exponent=(-h*c*self.B*self.J*(self.J+1)*(1/(k*T))).astype('float64')
        
        self.population=(2*self.J+1)*np.exp(exponent)
        #renormalize
        norm_pop=self.population/(np.sum(self.population))
        self.population=pd.Series(norm_pop)
        
        
        
    def allowed_combinations(self,Jup,Kup,Eup,Q_Branch=False):
        
        #determines which combinations are allowed based on J and K selection rules.
        df=pd.concat([self.J,self.K,self.E,self.population],axis=1)
        df.columns=["J","K","E","nJ"]
        df2=pd.DataFrame({"J'": Jup,"K'": Kup,"E'": Eup})
        E_list,E_prime_list,J_list,J_prime_list,K_list, K_prime_list,nJ_list=[],[],[],[],[],[],[]
        if self.symmetry_type=='spherical':
            if Q_Branch==False:
                DeltaJ=[-1,1]
            else:
                DeltaJ=[-1,0,1]
            #print("Starting search for allowed transitions")
            for i in range(len(df.index)):
                
                Jupp=df2["J'"].iloc[i]
                Eupp=df2["E'"].iloc[i]
                Kupp=df2["K'"].iloc[i]
                
                allowedJ=[Jupp-DelJ for DelJ in DeltaJ]
                df4=pd.DataFrame(df.loc[df.J.isin(allowedJ)])
                df4["J'"]=Jupp
                df4["E'"]=Eupp
                df4["K'"]=Kupp
                E_list.extend(df4["E"].values)
                J_list.extend(df4["J"].values)
                nJ_list.extend(df4["nJ"].values)
                E_prime_list.extend(df4["E'"].values)
                J_prime_list.extend(df4["J'"].values)
                K_list.extend(df4["K"].values)
                K_prime_list.extend(df4["K'"].values)
                

        elif 'symmetric' in self.symmetry_type:
            #print("Starting search for allowed transitions")
            for i in range(len(df.index)):

                Jupp=df2["J'"].iloc[i]
                Eupp=df2["E'"].iloc[i]
                Kupp=df2["K'"].iloc[i]
                DeltaK=[0]
                if df["K"].iloc[i]==0:
                    DeltaJ=[-1,1]
                else:
                    DeltaJ=[-1,0,1]
                    
                allowedJ=[Jupp-DelJ for DelJ in DeltaJ]
                allowedK=[Kupp-DelK for DelK in DeltaK]
                df4=pd.DataFrame(df.loc[df.J.isin(allowedJ) & df.K.isin(allowedK)])
                df4["J'"]=Jupp
                df4["E'"]=Eupp
                df4["K'"]=Kupp
                E_list.extend(df4["E"].values)
                J_list.extend(df4["J"].values)
                nJ_list.extend(df4["nJ"].values)
                E_prime_list.extend(df4["E'"].values)
                J_prime_list.extend(df4["J'"].values)
                K_list.extend(df4["K"].values)
                K_prime_list.extend(df4["K'"].values)
                
        else:
           print('symmetry type not yet available')
        
        
            
        df3=pd.DataFrame({"E": E_list, "E'": E_prime_list, "J": J_list, "J'": J_prime_list, "K": K_list, "K'": K_prime_list, "nJ": nJ_list})
        self.allowed_combo_data=df3
        
            
                


        
    def transition_freq_and_pop(self):
        if self.allowed_combo_data is None:
            print("Do not have relevant data")
        else:
            self.transition_freqs=self.allowed_combo_data["E'"]-self.allowed_combo_data["E"]
            self.allowed_combo_data['DeltaJ']=self.allowed_combo_data["J'"]-self.allowed_combo_data["J"]
           
            self.allowed_combo_data['nJ'].loc[self.allowed_combo_data['DeltaJ']==0] = self.allowed_combo_data["nJ"]*self.Q_scale
            self.allowed_combo_data['nJ'].loc[self.allowed_combo_data['DeltaJ']!=0] = self.allowed_combo_data["nJ"]*self.PR_scale
            
            
            
            self.transition_intensity=self.allowed_combo_data["nJ"]
            
            #find highest populated state
            max_pop_idx=self.transition_intensity.astype('float64').idxmax()
            self.highest_pop_state=self.allowed_combo_data["J"].iloc[max_pop_idx]
            #print("Jmax=" ,self.highest_pop_state)
        
        
    def plot_transitions(self):
        #print("Plotting the resulting transitions")
        plt.vlines(x=self.transition_freqs,ymax=self.transition_intensity,ymin=0,color='black',label='transitions')
  
        df=pd.DataFrame({"Trans_Freqs": self.transition_freqs, "Trans_Intens": self. transition_intensity, "Delta J": self.allowed_combo_data["J'"]-self.allowed_combo_data["J"]})
        df=df.sort_values(by=["Trans_Freqs"])
        self.sorted_freqs=df["Trans_Freqs"]
        self.sorted_intens=df["Trans_Intens"]
        
        plt.xlabel("Transition Frequency (1/cm)")
        plt.ylabel("Intensity")
        plt.title(str(self.name))
        plt.savefig(str(self.name)+'_'+self.symmetry_type+".pdf")


            
      
    def _rebin_data(self,X,Y,bin_size,Verbose=False):
        rebinned_x,rebinned_y=[],[]
        bins=np.arange(np.float(np.min(X))-5,np.float(np.max(X))+5, bin_size )
        
        inds=np.digitize(X,bins)

        df2=pd.DataFrame({"X":X,"Y": Y, "bin": inds})

        for i in range(1,len(bins)):

            search=df2.loc[df2["bin"]==i]

            if len(search.index)!=0:
                center_wave=bins[i]-bin_size/2
                sum_y=search["Y"].sum()

                rebinned_x.append(center_wave)
                rebinned_y.append(sum_y/bin_size)
                if Verbose==True:
                    print('Bin:',i,'from', bins[i-1], 'to', bins[i])
                    print('Wavenlength:', bins[i]-bin_size/2)
                    print(search)
                    print('Sum of Y:', sum_y)
                    print('Y To Plot:', sum_y/bin_size)
                    print('---------')
            else:
                center_wave=bins[i]-bin_size/2


                rebinned_x.append(center_wave)
                rebinned_y.append(0)
        
    
        return(rebinned_x,rebinned_y)

      
      
    def _CreateWave(self,X,Y,lambda0):
        c=const.c.to('km/s')
        kms=(u.km / u.s)
        
    #convert X into velocity space. del_sigma is measured in 1/cm (built into rotational_energies class) and lambda0 is measured in AA.
        V=[((c*(del_sigma*u.k*lambda0*u.AA)).to(u.km / u.s)).value for del_sigma in X]
        #################################################################
        #resample V and Y at 0.1km/s intervals
        bin_size=0.1
        X_sampled,Y_sampled=self._rebin_data(V,Y,bin_size)


        #################################################
        #Convert velocities back into wavelengths

        Final_X=[((lambda0*u.AA*(1+(V*kms/c))).to(u.AA)).value for V in X_sampled]
        
        return(Final_X,Y_sampled)
        
    def apply_voigt(self,lambda0,show_figure=False):
        from edibles.utils.voigt_profile import voigt_optical_depth
        X=self.transition_freqs
        Y=self.transition_intensity
        Wave,Intensity=self._CreateWave(X,Y,lambda0)
        final_tau=np.zeros(len(Wave))
        for i in range(len(Wave)):
            tau=voigt_optical_depth(Wave, lambda0=Wave[i],b=1,N=Intensity[i]*10**10,f=1,gamma=1e7,v_rad=0)
            final_tau=final_tau+tau
        self.spectrax=Wave
        self.spectray=final_tau
        if show_figure==True:
            plt.plot(Wave,final_tau,'k-',label='voigt profile applied')
            plt.xlabel('Wavelength $\mathrm{\AA}$')
            plt.ylabel('Tau')
            plt.show()
            
    def apply_radiative_transfer(self,show_figure=False):
        self.simple_rt_y=1-self.spectray
        self.full_rt_y=np.exp(-self.spectray)
        if show_figure==True:
            plt.plot(self.spectrax,self.simple_rt_y,'k-',label='after radiative transfer step')
            
            plt.xlabel('Wavelength $\mathrm{\AA}$')
            plt.ylabel('Optical Depth')
            plt.show()
        return()
        
    def smooth_spectra(self, lambda0,show_figure=False):
    
       
        dx=np.asarray(self.spectrax[1:])-np.asarray(self.spectrax[0:-1])
        d_lambda=lambda0/80000
        sigma_a=d_lambda/2.355
        sigma_p=sigma_a/dx[0]
        
        kernal=Gaussian1DKernel(sigma_p)
        convolved_y=convolve(self.full_rt_y,kernal,boundary='extend')
        
        if show_figure==True:
            plt.plot(self.spectrax, convolved_y/np.max(convolved_y), 'c-',label='after smoothing')
            plt.xlabel('Wavelength $\mathrm{\AA}$')
            plt.ylabel('Normalized Intensity')

        self.final_y=convolved_y/np.max(convolved_y)

        return()
        

        