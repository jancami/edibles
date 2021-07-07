# import statements
from edibles.utils.simulations.SimulatedContour import Simulated_Contour as sim
from edibles.utils.simulations.RotationalEnergies import WavelengthToWavenumber

# Import some modules
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
import pickle
from scipy.signal import argrelextrema

# Colors and style for plots
colormap = plt.cm.cool
plt.rcParams.update({'font.size': 10})
plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]


def animation(simulations, param_name, param_values, SymmetryType,
              path='Test_data', fps=5, plot_min=False):
    """Animate a series of simulations.

    Creates a GIF with the result of different profiles simulations obtained by varying
    some parameter(s).

    Args:
        simulations (list): List with simulations data.
        param_name (str): Name of varying variable(s).
        param_values (list): Values of parameters of each simulation.
        path (str): Directory to save image.
        fps (int, optional): Frames per second of animation. Defaults to 5.
        plot_min (bool, optional): To plot local minimum values or not. Defaults to False.
    """
    # Fig object and size.
    fig, ax = plt.subplots(figsize=(6, 6))

    # Find interval for x_aix values
    x_min = min([min(simulation[0]) for simulation in simulations])
    x_max = max([max(simulation[0]) for simulation in simulations])

    # Set axes limits.
    ax.set_xlim((x_min-1, x_max+1))

    y_min = min([min(simulation[1]) for simulation in simulations])
    ax.set_ylim(y_min, 1)

    # add grid to plot.
    ax.grid()

    # Plot to update in simulation.
    line, = ax.plot([], [], lw=2, alpha=1, color='black')

    # If plotting local minimum
    if plot_min:
        global vlines
        vlines = []
        scat = ax.scatter([0, 0, 0], [0, 0, 0])

    # Function to initializate animation.
    def init():
        line.set_data([0], [0])
        return (line,)

    # Main animation function.
    def animate(i):

        # Remove past local minimum plots
        if plot_min:
            global vlines
            for vline in vlines:
                vline.remove()
            vlines = []

        # Update plot.
        line.set_data(simulations[i][0], simulations[i][1])
        ax.set_title(f'B: {param_values[i][0]:.3f} T: {param_values[i][1]:.3f} ' +
                     f'Delta: {param_values[i][2]:.3f}')

        # Update local minimum markers.
        if plot_min:
            # Find and plot local minimum (PQR branches peaks)
            index = argrelextrema(simulations[i][1], np.less)
            scat.set_offsets(np.array([simulations[i][0][index], simulations[i][1][index]]).T)

            # Save vlines.
            for j in list(index)[0]:
                vlines.append(plt.axvline(x=simulations[i][0][j], ymin=0, ymax=1))

            # Return animation.
            return (line, scat)

        else:
            # Return animation.
            return (line, )

    # Animate
    anim = FuncAnimation(fig, animate, init_func=init,
                         frames=len(simulations), interval=50, blit=True)

    # Axes labels.
    ax.set_xlabel(r'Wavenumber ($cm^{-1}$)')
    ax.set_ylabel('Normalized Intensity')

    # Save animation.
    anim.save(f"{path}/"+SymmetryType+"_Changing"+param_name+'.gif',
              writer='imagemagick', fps=fps)

    plt.close()


def plot_simulations(simulations, param_name, param_values, SymmetryType, path='Test_data'):
    """Plot a set of simulations.

    Plot in a single graph a set of simulations obtained by varying some
    parameter(s).

    Args:
        simulations (list): List with simulations data.
        param_name (str): Name of varying variable(s).
        param_values (list): Values of parameters of each simulation.
        path (str): Directory to save image.
    """
    # Colours for plots
    colours = [colormap(i) for i in np.linspace(0, 0.9, len(param_values))]

    # Plot all the profiles.
    for i, simulation in enumerate(simulations):
        plt.plot(simulation[0], simulation[1], ls='-',
                 label=f'B: {param_values[i][0]:.4f} ' +
                 f'T: {param_values[i][1]:.1f} ' +
                 f'Delta: {param_values[i][2]:.3f}', color=colours[i])

    # Plot style
    plt.xlabel(r'Wavenumber ($cm^{-1}$)')
    plt.ylabel('Normalized Intensity')
    plt.legend(fontsize='x-small')
    plt.title('Parameter Survey - Changing '+param_name)
    plt.savefig(f"{path}/"+str(SymmetryType)+"_Changing"+param_name+".pdf", dpi=300)
    plt.show()
    plt.close()


def one_variable_survey(path, B_values, T_values, delta_values, Q_Branch=True,
                        B_default=0.02, T_default=10, delta_default=0):
    """Perform a parameter survey over individual parameters in a Spherical Symmetry.

    Generates plots and GIFs of surveys performed over individual parameters.
    (T, B and delta).

    Args:
        B_values (1darray): Array with values of B.
        T_values (1darray): Array with values of T.
        delta_values (1darray): Array with values of deltaB.
        Q_Branch (bool, optional): To consider or not the Q-branch. Defaults to True.
        B_default (flota, optional): Default value of B when constant. Defaults to 0.02
        T_default (flota, optional): Defaults to 10
        B_default (flota, optional): Defaults to 0.0
    """
    # Declare constants
    SymmetryType = 'Spherical'
    Jlimit = 200
    lambda0 = 6614
    Sightline = 'ParameterSurvey'
    Q_scale_init = 1

    # Dictionary to save simulations per parameter
    param_simulations = {'B': [], 'T': [], 'delta': []}

    # Dictionary to save parameters used for simulation
    parameters = {'B': [], 'T': [], 'delta': []}

    # Dictionary with parameter space
    param_values = {'B': B_values, 'T': T_values, 'delta': delta_values}

    # Varying all the parameters.
    for param in param_values.keys():

        # Iterate over values of parameter.
        for i in range(len(param_values[param])):

            # Initializate Q_scale
            Q_scale = Q_scale_init

            # Parameter values.
            if param == 'B':
                B = B_values[i]
                A = C = B
                T = T_default
                delta = delta_default

            elif param == 'T':
                B = B_default
                A = C = B
                T = T_values[i]
                delta = delta_default

            elif param == 'delta':
                B = B_default
                A = C = B
                T = T_default
                delta = delta_values[i]

            # Run simulation
            build = sim(A=A, B=B, C=C, Delta_A=delta*A, Delta_B=delta*B, Delta_C=delta*C,
                        Trot=T, Jlimit=Jlimit, Target=Sightline, lambda0=lambda0,
                        Q_Branch=Q_Branch, Q_scale=Q_scale)

            # Obtain wavenumbers.
            x_vals = WavelengthToWavenumber(np.asarray(build[0]))
            y_vals = build[1]

            # Save results and parameters.
            param_simulations[param].append([x_vals, y_vals])
            parameters[param].append([B, T, delta])

        plot_simulations(param_simulations[param], param,
                         parameters[param], SymmetryType, path=path)
        animation(param_simulations[param], param, parameters[param],
                  SymmetryType, path=path)


if __name__ == "__main__":

    # Define parameter ranges.
    T_values = np.linspace(10, 100, 10)
    B_values = 10**np.linspace(-4, -1, 10)
    delta_values = np.linspace(0, 0.05, 10)

    # Run simulations with and without the Q-branch
    one_variable_survey('Parameter_survey/QTrue', B_values, T_values, delta_values)
    one_variable_survey('Parameter_survey/QFalse',  B_values, T_values,
                        delta_values, Q_Branch=False)
