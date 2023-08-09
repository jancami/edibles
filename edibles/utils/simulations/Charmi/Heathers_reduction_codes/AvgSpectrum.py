import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import re
import seaborn as sns
from SRC.Functions import Signal_Noise_Calculator, weighted_average, WavelengthToWavenumber, WavenumberToWavelength
def CreateAverageSpectrum(DIB,Sightline,file_list):
    cmap = sns.color_palette("colorblind", as_cmap=True)
    h=0
    j=0
    df=pd.DataFrame()
    fig,ax=plt.subplots(constrained_layout=True)
    for file in file_list:
        if j==6:
            j=j+1
        
        
        m = re.search('_(.+?).csv', file).group(1)
        target_date=m[-10:]
        
        data=pd.read_csv(file,index_col=0)
        DIB_wavelength=data["Wavelength"]
        DIB_flux=data["Flux"]



        cont_x1=np.array(DIB_wavelength[-20:])
        cont_y1=np.array(DIB_flux[-20:])
        SN1,Fit1=Signal_Noise_Calculator(cont_x1,cont_y1)


        cont_x2=np.array(DIB_wavelength[:15])
        cont_y2=np.array(DIB_flux[:15])
        SN2,Fit2=Signal_Noise_Calculator(cont_x2,cont_y2)

        Signal_Noise=0.5*(SN1+SN2)
        New_Noise=DIB_flux/Signal_Noise
        
        uncertainty=np.full(DIB_wavelength.shape,np.mean(New_Noise),dtype=float)
        print(uncertainty)
        #plot spectra and RMS envelope
        
        ax.plot(DIB_wavelength,DIB_flux+h,color=cmap[1*j], label=target_date+" S/N: "+format(Signal_Noise,'.2f'))
        
        ax.plot(cont_x1,Fit1+h,'k-')
        ax.plot(cont_x1,Fit1+h-uncertainty[0],'k--')
        ax.plot(cont_x1,Fit1+h+uncertainty[0],'k--')
        ax.plot(cont_x2,Fit2+h,'k-')
        ax.plot(cont_x2,Fit2+h-uncertainty[0],'k--')
        ax.plot(cont_x2,Fit2+h+uncertainty[0],'k--')



        h=h+0.05
        df[target_date+"_data"]=DIB_flux
        df[target_date+"_error"]=uncertainty
        j=j+2

    weig_avg=[]
    weig_avg_err=[]
    print(df)
    for i in range(len(df.index)):
        values=df[df.columns[::2]].iloc[i].values
        error=df[df.columns[1::2]].iloc[i].values
        avg_weig=weighted_average(values,error)
        weig_avg.append(avg_weig[0])
        weig_avg_err.append(avg_weig[1])

    df["Weighted_Average"]=weig_avg
    df["Weighted_Average_Error"]=weig_avg_err

    cont_x1=np.array(DIB_wavelength[-20:])
    cont_y1=np.array(df["Weighted_Average"][-20:])
    SN1,Fit1=Signal_Noise_Calculator(cont_x1,cont_y1)


    cont_x2=np.array(DIB_wavelength[:15])
    cont_y2=np.array(df["Weighted_Average"][:15])
    SN2,Fit2=Signal_Noise_Calculator(cont_x2,cont_y2)

    Signal_Noise=0.5*(SN1+SN2)
    print(Signal_Noise)
    New_Noise=DIB_flux/Signal_Noise
    uncertainty=np.full(data["Wavelength"].shape,np.mean(New_Noise),dtype=float)


    ax.plot(DIB_wavelength,df["Weighted_Average"]+h,label='Weighted Average'+" S/N: "+format(Signal_Noise,'.2f'),linewidth=3,color='red')
    
    ax.plot(cont_x1,Fit1+h,'k-')
    ax.plot(cont_x1,Fit1+h-uncertainty[0],'k--')
    ax.plot(cont_x1,Fit1+h+uncertainty[0],'k--')
    ax.plot(cont_x2,Fit2+h,'k-')
    ax.plot(cont_x2,Fit2+h-uncertainty[0],'k--')
    ax.plot(cont_x2,Fit2+h+uncertainty[0],'k--')
    ax.legend(loc='lower right', fontsize='x-small')
    ax.set_xlabel('Wavelength ($\mathrm{\AA}$)')
    ax.set_ylabel('Normalized Flux')
    secax=ax.secondary_xaxis('top',functions=(WavelengthToWavenumber,WavenumberToWavelength))
    secax.set_xlabel(r'Wavenumber (cm$^{-1}$)' )
    plt.title(Sightline)
    plt.savefig('Model_Images/Comps/'+str(DIB)+'/'+str(DIB)+'_'+Sightline+'_WeightedAverage.pdf',bbox_inches='tight')
    #plt.show()
    plt.close()
    final=pd.DataFrame({"Wavelength":DIB_wavelength.to_numpy(), "Flux": df["Weighted_Average"].to_numpy()})
    final.to_csv("Data/AverageSpectraData/"+str(DIB)+"/"+Sightline+"_avg_spectra.csv")
    return(DIB_wavelength.to_numpy(),df["Weighted_Average"].to_numpy())