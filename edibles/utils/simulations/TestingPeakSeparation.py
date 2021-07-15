"""Test the branches peak separations.

Test the peaks separations between the simulations results and the theroetical values
for a molecule with a spherical symmetry. Constant B and delta values are considered,
but the test is performed for several T values. Plots with the error and the relative error
are generated.
"""

from edibles.utils.simulations.SimulatedContour import Simulated_Contour as sim
from edibles.utils.simulations.RotationalEnergies import WavelengthToWavenumber
from astropy import constants as const

import numpy as np
from scipy.signal import argrelextrema
import matplotlib.pyplot as plt

colormap = plt.cm.cool
aux = (const.k_B/(const.h*const.c.to('cm/s'))).value


def dist_RQ(T, B, delta):
    """Theoretical distance between R and Q branches for a molecule with given T, B and delta."""
    return (B + B*delta)*((2*aux*T/B)**(1/2) + 1)


def dist_QP(T, B, delta):
    """Theoretical distance between Q and P branches for a molecule with given T, B and delta."""
    return (B + B*delta)*((2*aux*T/B)**(1/2) - 1)


def dist_RP(T, B, delta):
    """Theoretical distance between R and P branches for a molecule with given T, B and delta."""
    return 2*(B + B*delta)*(2*aux*T/B)**(1/2)


def test_peak_separation(B, delta, Jlimit, N, Sightline, lambda0, Q_Branch, Q_scale):
    """Test branches peak separations between the simulation results and the theroetical values.

    The test is performed for a molecule with known B and delta values, and is tested it with
    several values of T. A plot of the simulation's profiles is generated.
    Also a plot with the error and relative error.

    Args:
        B (float): Value of rotational constant.
        delta (float): Percentage of change of the rotational constant to the upper states.
        Jlimit (int): Upper limit to consider of the rotational quantum number J.
        N (int): Number of splits in the parameter range for the temperature (10k to 90k).
        Sightline (str): Target sighline.
        lambda0 (float): Center wavelength of profile.
        Q_Branch (bool): Wheter or not to consider the Q branch.
        Q_scale (float): Scale for the Q branch (if considered).
    """
    # Only for spherical tops.
    A = C = B

    # Temperature values to try
    T_values = np.linspace(10, 90, N)

    # list to store results.
    x_vals = []
    y_vals = []

    # Run simulations
    for T in T_values:
        # Run simulation
        build = sim(A=A, B=B, C=C, Delta_A=delta*A, Delta_B=delta*B, Delta_C=delta*C,
                    Trot=T, Jlimit=Jlimit, Target=Sightline, lambda0=lambda0,
                    Q_Branch=Q_Branch, Q_scale=Q_scale)

        # Save and convert wavelength to wavenumbers.
        x_vals.append(WavelengthToWavenumber(np.asarray(build[0])) -
                      WavelengthToWavenumber(lambda0))
        y_vals.append(build[1])

    # List to save local minim.
    local_min = []

    for i in range(N):
        # Save local minima per simulation
        local_min.append(argrelextrema(y_vals[i], np.less))

    # Create colour for profiles.
    colours = [colormap(i) for i in np.linspace(0, 0.9, N)]

    for i in range(N):
        # Plot profiles.
        plt.plot(x_vals[i], y_vals[i], label=f'T = {T_values[i]}', alpha=0.9, color=colours[i])
        # Scatter minima points.
        plt.scatter(x_vals[i][local_min[i][0]], y_vals[i][local_min[i][0]], color=colours[i])

    # Plot details.
    plt.ylim((0.8, 1.01))
    plt.legend()
    plt.xlabel(r'Wavenumber (cm$^{-1}$)')
    plt.ylabel('Normalized Intensity')
    plt.title(f'Testing PQR peak separation. B = {B}, delta = {delta}')
    plt.grid()
    plt.show()

    # Lists to save error and realative errors for each pair of branches per simulation.
    err_RQ, err_QP, err_RP = [], [], []
    rel_RQ, rel_QP, rel_RP = [], [], []

    # For every simulation.
    for i in range(N):

        # Compute simulation peak position.
        v_P = x_vals[i][local_min[i][0][2]]
        v_Q = x_vals[i][local_min[i][0][1]]
        v_R = x_vals[i][local_min[i][0][0]]

        # Simulation's peaks distance.
        sim_RQ = v_R - v_Q
        sim_QP = v_Q - v_P
        sim_RP = v_R - v_P

        # Theoretical values of distance.
        the_RQ = dist_RQ(T_values[i], B, delta)
        the_QP = dist_QP(T_values[i], B, delta)
        the_RP = dist_RP(T_values[i], B, delta)

        # Error.
        err_RQ.append(the_RQ-sim_RQ)
        err_QP.append(the_QP-sim_QP)
        err_RP.append(the_RP-sim_RP)

        # Relative error.
        rel_RQ.append((the_RQ-sim_RQ)/the_RQ*100)
        rel_QP.append((the_QP-sim_QP)/the_QP*100)
        rel_RP.append((the_RP-sim_RP)/the_RP*100)

    # Plot of errors.
    fig, axs = plt.subplots(2, 1, sharex=True)

    axs[0].plot(T_values, err_RQ, label='R-Q')
    axs[0].plot(T_values, err_QP, label='Q-P')
    axs[0].plot(T_values, err_RP, label='R-P')

    axs[0].legend()
    axs[0].set_ylabel('Error')
    axs[0].set_title('Theory - Simulation')
    axs[0].grid()

    axs[1].plot(T_values, rel_RQ, label='R-Q')
    axs[1].plot(T_values, rel_QP, label='Q-P')
    axs[1].plot(T_values, rel_RP, label='R-P')

    axs[1].legend()
    axs[1].set_ylabel('Relative error (%)')
    axs[1].grid()

    plt.xlabel('Temperature (K)')
    plt.suptitle(f'Testing PQR peak separation. B = {B}, delta = {delta}')
    plt.show()


if __name__ == "__main__":

    test_peak_separation(B=0.02, delta=0.005, Jlimit=200, N=9, Sightline='Testing P-R distance',
                         lambda0=6614, Q_Branch=True, Q_scale=0.5)
