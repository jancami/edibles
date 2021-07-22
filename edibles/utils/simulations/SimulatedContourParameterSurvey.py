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


def animation(simulations, param_name, param_values,
              lambda0, path='Test_data', fps=5, plot_min=False):
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
    ax.set_xlim((-10, +10))

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
        ax.set_title(f'A: {param_values[i]["A"]:.4f} B: {param_values[i]["B"]:.4f}\
                     C: {param_values[i]["C"]:.4f} T: {param_values[i]["T"]:.2f}\
                     Delta: {param_values[i]["delta"]:.4f}')

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
    anim.save(f"{path}/_Changing"+param_name+'.gif',
              writer='imagemagick', fps=fps)

    plt.close()


def plot_simulations(simulations, param_name, param_values,
                     lambda0, path='Test_data', save=False):
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
                 label=f'{param_name} = {param_values[i][param_name]:.5f}', color=colours[i])

    # Plot style
    plt.xlim((-10, 10))

    plt.xlabel(r'Wavenumber ($cm^{-1}$)')
    plt.ylabel('Normalized Intensity')
    plt.legend(fontsize='x-small')
    plt.title(f'Changing {param_name}. Initial values A = {param_values[0]["A"]},\
              B = {param_values[0]["B"]}, C = {param_values[0]["C"]}, T = {param_values[0]["T"]},\
              delta = {param_values[0]["delta"]}')

    if save:
        plt.savefig(f"{path}/_Changing"+param_name+".pdf", dpi=300)
    plt.show()
    plt.close()


def one_variable_survey(path, A_values=None, B_values=None, C_values=None, T_values=None,
                        delta_values=None, Q_Branch=True, B_default=0.02, T_default=10,
                        delta_default=0.001, save=False, load=False, save_fig=False, anim=False):
    """Perform a parameter survey over individual parameters in a Spherical Symmetry.

    Generates plots and GIFs of surveys performed over individual parameters.
    (T, A, B, C, and delta).

    Args:
        A_values (1darray): Array with values of A. Ignored if empty. Defaults to None.
        B_values (1darray): Array with values of B. Ignored if empty. Defaults to None.
        C_values (1darray): Array with values of C. Ignored if empty. Defaults to None.
        T_values (1darray): Array with values of T. Ignored if empty. Defaults to None.
        delta_values (1darray): Array with values of deltaB. Ignored if empty.  Defaults to None.
        Q_Branch (bool, optional): To consider or not the Q-branch. Defaults to True.
        B_default (float, optional): Default value of B when constant. Defaults to 0.02
        T_default (float, optional): Defaults to 10
        B_default (float, optional): Defaults to 0.0
        save_fig (bool, optional): Defaults to False.
        anim (bool, optional): Create GIF or not. Defaults to False.
    """
    # Declare constants
    Jlimit = 200
    lambda0 = 6614
    Sightline = 'ParameterSurvey'
    Q_scale_init = 1

    # Dictionary to save simulations per parameter varying.
    param_simulations = {'A': [], 'B': [], 'C': [], 'T': [], 'delta': []}

    # Dictionary to save parameters used for simulation
    parameters = {'A': [], 'B': [], 'C': [], 'T': [], 'delta': []}

    # Dictionary with parameter space
    param_values = {'A': A_values, 'B': B_values, 'C': C_values,
                    'T': T_values, 'delta': delta_values}

    # Varying all the parameters.
    for param in param_values.keys():
        if param_values[param] is not None:
            # Load pre-calculated simulations
            if load:
                with open(f"{path}/_Changing{param}_.bin", "rb") as data:
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
                    if param == 'A':
                        B = B_default
                        C = B
                        A = A_values[i]
                        T = T_default
                        delta = delta_default

                    if param == 'B':
                        B = B_values[i]
                        A = C = B
                        T = T_default
                        delta = delta_default

                    if param == 'C':
                        B = B_default
                        A = B
                        C = C_values[i]
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
                    x_vals = WavelengthToWavenumber(np.asarray(build[0])) - \
                        WavelengthToWavenumber(lambda0)
                    y_vals = build[1]

                    # Save results and parameters.
                    param_simulations[param].append([x_vals, y_vals])
                    parameters[param].append({'A': A, 'B': B, 'C': C, 'T': T, 'delta': delta})

            # Generates GIF and plot.
            plot_simulations(param_simulations[param], param, parameters[param],
                             lambda0, path=path, save=save_fig)

            if anim:
                animation(param_simulations[param], param, parameters[param],
                          lambda0, path=path)

            # Save simulation to future reference.
            if save:
                with open(f"{path}/_Changing{param}_.bin", "wb") as output:
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
                        B_default=0.02, T_default=10, delta_default=0, save=False,
                        load=False, plot_wavenumber=True, Q_scale=1):

    # Declare constants
    SymmetryType = 'Spherical'
    Jlimit = 200
    lambda0 = 6614
    Sightline = 'ParameterSurvey'

    # Dictionary to save simulations
    param_simulations = {'BT': {},
                         'Tdelta': {},
                         'Bdelta': {}}

    # Dictionary to save parameters used in simulations.
    parameters = {'BT': {},
                  'Tdelta': {},
                  'Bdelta': {}}

    # Dictionary with parameter space
    param_values = {'BT': [{'B': B_values, 'T': T_values}, {'delta': delta_default}],
                    'Tdelta': [{'T': T_values, 'delta': delta_values}, {'B': B_default}],
                    'Bdelta': [{'B': B_values, 'delta': delta_values}, {'T': T_default}]}

    # Load pre-calculated simulations
    if load:
        with open(f"{path}/{SymmetryType}_2D_survey.bin", "rb") as data:
            # Extract data
            build = pickle.load(data)
            param_simulations = build[0]
            parameters = build[1]

# Run simulation.
    else:

        for combinations in param_values.keys():

            parameter_space = param_values[combinations][0]
            variables = list(parameter_space.keys())

            # Simulation
            for value1 in parameter_space[variables[0]]:

                param_simulations[combinations][value1] = {}
                parameters[combinations][value1] = {}

                for value2 in parameter_space[variables[1]]:

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
                    param_simulations[combinations][value1][value2] = [x_vals, y_vals]
                    parameters[combinations][value1][value2] = [B, T, delta]

    # Save simulation to future reference.
    if save:
        with open(f"{path}/{SymmetryType}_2D_survey.bin", "wb") as output:
            pickle.dump([param_simulations, parameters], output)

    plot_grid(param_simulations, parameters, param_values, lambda0,
              path, SymmetryType, plot_wavenumber)


def plot_grid(simulations, parameters, param_values,
              lambda0, path, SymmetryType, plot_wavenumber):

    a = WavelengthToWavenumber(lambda0)

    for combinations in param_values.keys():

        parameter_space = param_values[combinations][0]
        variables = list(parameter_space.keys())

        rows = len(list(parameter_space[variables[0]]))
        cols = len(list(parameter_space[variables[1]]))

        fig, axs = plt.subplots(nrows=rows, ncols=cols, sharex=True,
                                sharey=True, figsize=(rows*2, cols*2))

        for i, value1 in enumerate(parameter_space[variables[0]]):
            for j, value2 in enumerate(parameter_space[variables[1]]):
                x = simulations[combinations][value1][value2][0]
                y = simulations[combinations][value1][value2][1]

                if plot_wavenumber:
                    axs[i, j].set_xlim((-5, +5))
                    x -= a

                else:
                    x = WavelengthToWavenumber(x)
                    axs[i, j].set_xlim((lambda0-2, lambda0+2))

                axs[i, j].plot(x, y)
                axs[i, j].grid(alpha=0.5)

                if i == rows-1:
                    if plot_wavenumber:
                        axs[i, j].set_xlabel(r'Wavenumber (cm$^{-1}$)')
                    else:
                        axs[i, j].set_xlabel('Wavelength $\mathrm{\AA}$')

                if i == 0:
                    axs[i, j].set_title(f'{variables[1]} = {value2:.4f}')

                if j == cols-1:
                    axs[i, j].set_ylabel(f'{variables[0]} = {value1:.4f}')
                    axs[i, j].yaxis.set_label_position("right")

                if j == 0:
                    axs[i, j].set_ylabel('Normalized Intensity')

        constant = list(param_values[combinations][1].keys())[0]
        value = list(param_values[combinations][1].values())[0]

        plt.suptitle(f'2D parameter survey with {variables[0]} and {variables[1]}. ' +
                     f'Constant {constant} = {value:.4f}')
        plt.tight_layout()
        plt.savefig(f"{path}/"+str(SymmetryType)+"_Changing"+combinations+".pdf", dpi=300)
        plt.show()


if __name__ == "__main__":

    if False:

        # Define parameter ranges.
        T_values = np.linspace(10, 100, 5)
        B_values = 10**np.linspace(-4, -1, 5)
        delta_values = np.linspace(0, 0.05, 5)

        # Default values when constant.
        T_default = np.linspace(10, 100, 3)
        B_default = 10**np.linspace(-4, -1, 3)
        delta_default = np.linspace(0, 0.05, 3)

        # Run 2D surveys with three different default values
        for i in range(3):
            print(i, 1)
            # Simulation without Q branch.
            two_variable_survey(path=f'Parameter_survey/2D/QFalse/default{i+1}', B_values=B_values,
                                T_values=T_values, delta_values=delta_values, Q_Branch=False,
                                B_default=B_default[i], T_default=T_default[i],
                                delta_default=delta_default[i], load=True)
            print(i, 2)
            # Simulation with Q branch
            two_variable_survey(path=f'Parameter_survey/2D/QTrue/default{i+1}', B_values=B_values,
                                T_values=T_values, delta_values=delta_values, Q_Branch=True,
                                B_default=B_default[i], T_default=T_default[i],
                                delta_default=delta_default[i], load=True)
            print(i, 3)
            # Simultion with scaled Q branch
            two_variable_survey(path=f'Parameter_survey/2D/QScale/default{i+1}', B_values=B_values,
                                T_values=T_values, delta_values=delta_values, Q_Branch=True,
                                B_default=B_default[i], T_default=T_default[i],
                                delta_default=delta_default[i], load=True, Q_scale=0.1)

    if True:

        # Define parameter ranges.
        A_values = np.arange(0.05, 0.5, 0.05)
        one_variable_survey(path='Parameter_survey', A_values=A_values,
                            save=True, save_fig=True, anim=True)

        # Define parameter ranges.
        C_values = np.arange(0.0001, 0.015, 0.002)
        one_variable_survey(path='Parameter_survey', C_values=C_values,
                            save=True, save_fig=True, anim=True)
