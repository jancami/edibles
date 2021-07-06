#import the stuff from my class
from RotationalEnergies import Rotational_Energies
import matplotlib.pyplot as plt
import pandas as pd
def Simulated_Contour(A,Delta_A,B,Delta_B,C,Delta_C,Trot,Jlimit, Name,Q_scale=1,PR_scale=1, Q_Branch=True,lambda0=0):
    re_low = Rotational_Energies(A=A,B=B,C=C, Name=Name,Q_scale=Q_scale,PR_scale=PR_scale)
    if re_low.flag:
        print("Can't deal with this molecule yet")
        return([0],[0])
    else:
        re_low.rotational_energies(Jlimit=Jlimit)
        re_low.boltzmann(T=Trot)
        re_up = Rotational_Energies(A=A+Delta_A,B=B+Delta_B,C=C+Delta_C,Name=Name,Q_scale=Q_scale,PR_scale=PR_scale)


        
        re_up.rotational_energies(Jlimit=Jlimit)
        re_low.allowed_combinations(Jup=re_up.J,Kup=re_up.K,Eup=re_up.E, Q_Branch=Q_Branch)
        re_low.transition_freq_and_pop()
    #    plot1=plt.figure(1)
    #    re_low.plot_transitions()
    #    plt.legend()
    #    plot2=plt.figure(2)
        re_low.apply_voigt(lambda0=lambda0,show_figure=False)
    #    plt.legend()
    #    plot3=plt.figure(3)
        re_low.apply_radiative_transfer(show_figure=False)
    #    plt.legend()
    #    plot4=plt.figure(4)
        re_low.smooth_spectra(lambda0=lambda0,show_figure=False)
    #    plt.legend()
    #    plt.show()
    
   
        return(re_low.spectrax,re_low.final_y)
    
    
    
if __name__ == "__main__":


    sim=Simulated_Contour(A=42e-3,B=42e-3,C=42e-3,Delta_A=42e-3*(0.001),Delta_B=42e-3*(0.001),Delta_C=42e-3*(0.001),Trot=15,Jlimit=50, Name='Test',lambda0=6614, Q_Branch=True)
    #
    plt.plot(sim[0],sim[1],'k-')
    plt.xlabel(r'Wavelength ($\mathrm{\AA}$)')
    plt.ylabel('Normalized Intensity')
    plt.title('Simulated Spectrum')
    plt.show()
