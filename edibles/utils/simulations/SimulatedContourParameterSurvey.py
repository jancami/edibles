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


def animation(simulations, param_name, param_values, SymmetryType, lambda0,
              path='Test_data', fps=5, plot_min=False):
    """Animate a series of simulations.

    Creates a GIF with the result of different profiles simulations obtained by varying
    some parameter(s).

    Args:
        simulations (list): List with simulations data.
        param_name (str): Name of varying variable(s).
        param_values (list): Values of parameters of each simulation.
        lambda0 (float): Center wavelength of profile.
        path (str): Directory to save image.
        fps (int, optional): Frames per second of animation. Defaults to 5.
        plot_min (bool, optional): To plot local minimum values or not. Defaults to False.
    """
    # Fig object and size.
    fig, ax = plt.subplots(figsize=(6, 6))

    # Find interval for x_aix values
    # x_min = min([min(simulation[0]) for simulation in simulations])
    # x_max = max([max(simulation[0]) for simulation in simulations])
    # ax.set_xlim((x_min-1, x_max+1))

    # Set axes limits.
    a = WavelengthToWavenumber(lambda0)
    ax.set_xlim((a-10, a+10))

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


def plot_simulations(simulations, param_name, param_values,
                     SymmetryType, lambda0, path='Test_data'):
    """Plot a set of simulations.

    Plot in a single graph a set of simulations obtained by varying some
    parameter(s).

    Args:
        simulations (list): List with simulations data.
        param_name (str): Name of varying variable(s).
        param_values (list): Values of parameters of each simulation.
        lambda0 (float): Center wavelength of profile.
        path (str): Directory to save image.
    """
    # Colours for plots
    colours = [colormap(i) for i in np.linspace(0, 0.9, len(param_values))]

    # Plot all the profiles.
    for i, simulation in enumerate(simulations):
        plt.plot(simulation[0], simulation[1], ls='-', alpha=0.7,
                 label=f'B: {param_values[i][0]:.4f} ' +
                 f'T: {param_values[i][1]:.1f} ' +
                 f'Delta: {param_values[i][2]:.3f}', color=colours[i])

    # Plot style
    a = WavelengthToWavenumber(lambda0)
    plt.xlim((a-10, a+10))

    plt.xlabel(r'Wavenumber ($cm^{-1}$)')
    plt.ylabel('Normalized Intensity')
    plt.legend(fontsize='x-small')
    plt.title('Parameter Survey - Changing '+param_name)
    plt.savefig(f"{path}/"+str(SymmetryType)+"_Changing"+param_name+".pdf", dpi=300)
    plt.show()
    plt.close()


def one_variable_survey(path, B_values, T_values, delta_values, Q_Branch=True,
                        B_default=0.02, T_default=10, delta_default=0, save=False, load=False):
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

    # Dictionary to save simulations per parameter varying.
    param_simulations = {'B': [], 'T': [], 'delta': []}

    # Dictionary to save parameters used for simulation
    parameters = {'B': [], 'T': [], 'delta': []}

    # Dictionary with parameter space
    param_values = {'B': B_values, 'T': T_values, 'delta': delta_values}

    # Varying all the parameters.
    for param in param_values.keys():

        # Load pre-calculated simulations
        if load:
            with open(f"{path}/{SymmetryType}_Changing{param}_.bin", "rb") as data:
                # Extract data
                build = pickle.load(data)
                param_simulations[param] = build[0]
                parameters[param] = build[1]

        # Run simulation.
        else:
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

        # Generates GIF and plot.
        plot_simulations(param_simulations[param], param,
                         parameters[param], SymmetryType, lambda0, path=path)
        animation(param_simulations[param], param, parameters[param],
                  SymmetryType, lambda0, path=path)

        # Save simulation to future reference.
        if save:
            with open(f"{path}/{SymmetryType}_Changing{param}_.bin", "wb") as output:
                pickle.dump([param_simulations[param], parameters[param]], output)


def assign_two_variables(value1, value2, name1, name2,
                         B_default, T_default, delta_default):
    """Assign parameters values in the two-parameter survey.

    Args:
        value1 (float): Value of the first variable.
        value2 (float): Value of the second variable.
        name1 (str): Name of first variable (B, T or delta).
        name2 (str): Name of second variable.
        B_default (float): Default value of B when constant.
        T_default (float): Default value of T.
        delta_default (float): Default value of delta.

    Returns:
        B (float): Final value of B
        T (float): Value of T.
        delta (float): Value of delta.
    """
    if name1 == 'B':
        B = value1
        if name2 == 'T':
            T = value2
            delta = delta_default

        elif name2 == 'delta':
            delta = value2
            T = T_default

    elif name1 == 'T':
        T = value1
        if name2 == 'B':
            B = value2
            delta = delta_default

        elif name2 == 'delta':
            delta = value2
            B = B_default

    elif name1 == 'delta':
        delta = value1
        if name2 == 'B':
            B = value2
            T = T_default

        elif name2 == 'T':
            T = value2
            B = B_default

    return B, T, delta


def two_variable_survey(path, B_values, T_values, delta_values, Q_Branch=True,
                        B_default=0.02, T_default=10, delta_default=0, save=False, load=False):

    # Declare constants
    SymmetryType = 'Spherical'
    Jlimit = 200
    lambda0 = 6614
    Sightline = 'ParameterSurvey'
    Q_scale = 1

    # Dictionary to save simulations
    param_simulations = {'BT': {},
                         'Tdelta': {},
                         'Bdelta': {}}

    # Dictionary to save parameters used in simulations.
    parameters = {'BT': {},
                  'Tdelta': {},
                  'Bdelta': {}}

    # Dictionary with parameter space
    param_values = {'BT': {'B': B_values, 'T': T_values},
                    'Tdelta': {'T': T_values, 'delta': delta_values},
                    'Bdelta': {'B': B_values, 'delta': delta_values}}

    for combinations in param_values.keys():
        # Load pre-calculated simulations
        if load:
            with open(f"{path}/{SymmetryType}_Changing{combinations}_.bin", "rb") as data:
                # Extract data
                build = pickle.load(data)
                param_simulations[combinations] = build[0]
                parameters[combinations] = build[1]

        # Run simulation.
        else:

            variables = list(param_values[combinations].keys())

            # Simulation
            for value1 in param_values[combinations][variables[0]]:

                param_simulations[combinations][value1] = {}
                parameters[combinations][value1] = {}
                for value2 in param_values[combinations][variables[1]]:

                    B, T, delta = assign_two_variables(value1, value2, variables[0], variables[1],
                                                       B_default, T_default, delta_default)
                    A = C = B

                    build = sim(A=A, B=B, C=C, Delta_A=delta*A, Delta_B=delta*B, Delta_C=delta*C,
                                Trot=T, Jlimit=Jlimit, Target=Sightline, lambda0=lambda0,
                                Q_Branch=Q_Branch, Q_scale=Q_scale)

                    # Obtain wavenumbers.
                    x_vals = WavelengthToWavenumber(np.asarray(build[0]))
                    y_vals = build[1]

                    # Save results and parameters.
                    param_simulations[combinations][value1][value2] = {
                        'simulation': [x_vals, y_vals], 'parameters': [B, T, delta]}
                    parameters[combinations][value1][value2] = [B, T, delta]

        # Save simulation to future reference.
        if save:
            with open(f"{path}/{SymmetryType}_Changing{combinations}_.bin", "wb") as output:
                pickle.dump([param_simulations[combinations], parameters[combinations]], output)


if __name__ == "__main__":

    # Define parameter ranges.
    T_values = np.linspace(10, 100, 3)
    B_values = 10**np.linspace(-4, -1, 3)
    delta_values = np.linspace(0, 0.05, 3)

    # Run simulations with and without the Q-branch
    # one_variable_survey('Parameter_survey/QTrue', B_values, T_values, delta_values, load=True)
    # one_variable_survey('Parameter_survey/QFalse',  B_values, T_values, delta_values, Q_Branch=False, load=True)

    two_variable_survey('Parameter_survey/2D/QTrue', B_values, T_values, delta_values, save=True)
