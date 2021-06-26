"""A one line summary of the module or program, terminated by a period.

Leave one blank line.  The rest of this docstring should contain an
overall description of the module or program.  Optionally, it may also
contain a brief description of exported classes and functions and/or usage
examples.
"""

#import the stuff from my class
from RotationalEnergies import Rotational_Energies
import matplotlib.pyplot as plt
import pandas as pd

def Simulated_Contour(A,Delta_A,B,Delta_B,C,Delta_C,Trot,Jlimit, Name,Q_scale=1,PR_scale=1, Q_Branch=True,lambda0=0):
    """Summary of function...
    
    Longer information...

    Args: 
        A (TYPE): DESCRIPION
        Delta_A (TYPE): DESCRIPION
        B (TYPE): DESCRIPION
        Delta_B (TYPE): DESCRIPION
        C (TYPE): DESCRIPION
        Delta_C (TYPE): DESCRIPION
        Trot (TYPE): DESCRIPION
        Jlimit (TYPE): DESCRIPION
        Name (TYPE): DESCRIPION
        Q_scale (TYPE): DESCRIPION. Optional; the default is 1.
        PR_scale (TYPE): DESCRIPION. Optional; the default is 1.
        Q_Branch (TYPE): DESCRIPION. Optional; the default is True.
        lambda0 (TYPE): DESCRIPION. Optional; the default is 0.

    Returns: 
        re_low.spectrax (TYPE): DESCRIPTION
        re_low.final_y (TYPE): DESCRIPTION

    """
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
