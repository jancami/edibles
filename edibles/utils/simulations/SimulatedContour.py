"""Script to generate a simulated contour of a molecule with known rotational constants."""

from edibles.utils.simulations.RotationalEnergies import Rotational_Energies
import matplotlib.pyplot as plt


def Simulated_Contour(A, Delta_A, B, Delta_B, C, Delta_C, Trot, Jlimit, Target, Q_scale=1,
                      PR_scale=1, Q_Branch=True, lambda0=0):
    """Generate simulated contour.

    Args:
        A (float): Rotational constant of the first rotational axis.
        Delta_A (float): Difference between the values of the first rotational constant of the 
            upper and lower states.
        B (float): Constant of the second axis.
        Delta_B (float): Difference of the second constant.
        C (float): Constant of the third axis.
        Delta_C (float): Difference of the third constant.
        Trot (float): Temperature (Kelvin degrees).
        Jlimit (int): Upper bound of the first rotational quantum number J.
        Target (str): Name of the target sightline.
        Q_scale (float, optional): Scale of the Q-branch. Default to 1.
        PR_scale (float, optional): Scale of the P-branch and R-branch. Default to 1.
        Q_Branch (bool, optional): Wheter to consideror not the Q-branch. Default to True.
        lambda0 (float, optional): Center wavelength of DIB (Angstrom). Default to 0.

    Returns:
        re_low.spectrax (1darray): Resulting spectrum.
        re_low.final_y (1darray): Intensity of spectrum.
    """
    # Generate class object.
    re_low = Rotational_Energies(A=A, B=B, C=C, Target=Target,
                                 Q_scale=Q_scale, PR_scale=PR_scale)

    # Check for available symmetries.
    if re_low.flag:
        print("Can't deal with this molecule yet")
        return([0], [0])
    else:
        # Get rotational energies.
        re_low.rotational_energies(Jlimit=Jlimit)

        # Get energies populations.
        re_low.boltzmann(T=Trot)

        # Define new class to transitionate.
        re_up = Rotational_Energies(A=A+Delta_A, B=B+Delta_B, C=C+Delta_C,
                                    Target=Target, Q_scale=Q_scale, PR_scale=PR_scale)

        # Get rotational energies of new class
        re_up.rotational_energies(Jlimit=Jlimit)

        # Get allowed combinations between the two clases (ie states).
        re_low.allowed_combinations(Jup=re_up.J, Kup=re_up.K,
                                    Eup=re_up.E, Q_Branch=Q_Branch)

        # Get transition frequencies and populations.
        re_low.transition_freq_and_pop()
    #    plot1=plt.figure(1)
    #    re_low.plot_transitions()
    #    plt.legend()
    #    plot2=plt.figure(2)

        # Apply voigt profile.
        re_low.apply_voigt(lambda0=lambda0, show_figure=False)
    #    plt.legend()
    #    plot3=plt.figure(3)

        # Apply radiative transfer
        re_low.apply_radiative_transfer(show_figure=False)
    #    plt.legend()
    #    plot4=plt.figure(4)

        # Apply 1D Gaussian kernel.
        re_low.smooth_spectra(lambda0=lambda0, show_figure=False)
    #    plt.legend()
    #    plt.show()

        return(re_low.spectrax, re_low.final_y)


if __name__ == "__main__":

    # Perform simulation.
    sim = Simulated_Contour(A=42e-3, B=42e-3, C=42e-3, Delta_A=42e-3*(0.001),
                            Delta_B=42e-3*(0.001), Delta_C=42e-3*(0.001), Trot=15,
                            Jlimit=50, Target='Test', lambda0=6614, Q_Branch=True)

    # Plot result.
    plt.plot(sim[0], sim[1], 'k-')
    plt.xlabel(r'Wavelength ($\mathrm{\AA}$)')
    plt.ylabel('Normalized Intensity')
    plt.title('Simulated Spectrum')
    plt.show()
