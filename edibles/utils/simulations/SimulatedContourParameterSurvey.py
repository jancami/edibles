"""Script to perform paramter survey of simulations.

This script provides methods to perform one and two-parameter surveys. It's possible
to vary over parameter space of rotational constants (A, B, C), change in rotational
constants (deltaA, deltaB, deltaC), and temperature. It also includes some methods
and options to generate visualizations.
"""

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
plt.rcParams['text.latex.preamble'] = r"\usepackage{amsmath}"


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
    ax.set_xlim((-4, +4))

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
        ax.set_title(f'A: {param_values[i]["A"]:.4f} B:{param_values[i]["B"]:.4f}' +
                     f'C: {param_values[i]["C"]:.4f} T:{param_values[i]["T"]:.2f}\n' +
                     f'deltaA: {param_values[i]["deltaA"]:.5f}, ' +
                     f'deltaB: {param_values[i]["deltaB"]:.5f}, ' +
                     f'deltaC: {param_values[i]["deltaC"]:.5f}')

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
    anim = FuncAnimation(fig, animate, init_func=init, repeat_delay=1000,
                         frames=len(simulations), interval=500, blit=True)

    # Axes labels.
    ax.set_xlabel(r'Wavenumber ($cm^{-1}$)')
    ax.set_ylabel('Normalized Intensity')

    # Save animation.
    anim.save(f"{path}/_Changing"+param_name+'.gif',
              writer='imagemagick', fps=fps)

    plt.close()


def plot_simulations(simulations, param_name, param_values,
                     lambda0, path='Test_data', save=False, SymmetryType='Prolate'):
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
    plt.xlim((-4, 4))

    plt.xlabel(r'Wavenumber ($cm^{-1}$)')
    plt.ylabel('Normalized Intensity')
    plt.legend(fontsize='x-small')
    title = (f'Changing {param_name}. Initial values A:{param_values[0]["A"]:.4f}\n' +
             f'B:{param_values[0]["B"]:.4f}, C:{param_values[0]["C"]:.4f}, ' +
             f'T:{param_values[0]["T"]:.2f}\ndeltaA = {param_values[0]["deltaA"]:.5f},' +
             f' deltaB = {param_values[0]["deltaB"]:.4f}, ' +
             f'deltaC = {param_values[0]["deltaC"]:.4f}')
    plt.title(title)
    plt.tight_layout()

    if save:
        plt.savefig(f"{path}/{SymmetryType}_Changing"+param_name+".pdf", dpi=300)
    # plt.show()
    plt.close()


def one_variable_survey(path, A_values=None, B_values=None, C_values=None, T_values=None,
                        delta_A_values=None, delta_B_values=None, delta_C_values=None,
                        Q_Branch=True, A_default=0.02, B_default=0.02, C_default=0.02,
                        T_default=10, delta_A_default=0.001, delta_B_default=0.001,
                        delta_C_default=0.001, Jlimit=200, lambda0=6614, save=False, load=False,
                        save_fig=False, anim=False, SymmetryType='Prolate'):
    """Perform a parameter survey over individual parameters in a Spherical Symmetry.

    Generates plots and GIFs of surveys performed over individual parameters.
    (T, A, B, C, delta_A, delta_B and delta_C).

    Args:
        A_values (1darray): Array with values of A. Ignored if empty. Defaults to None.
        B_values (1darray): Values of B. Ignored if empty. Defaults to None.
        C_values (1darray): Values of C. Ignored if empty. Defaults to None.
        T_values (1darray): Values of T. Ignored if empty. Defaults to None.
        delta_A_values (1darray): Array with values of deltaA. Ignored if empty.  Defaults to None.
        delta_B_values (1darray): Values of deltaB. Ignored if empty.  Defaults to None.
        delta_C_values (1darray): Values of deltaC. Ignored if empty.  Defaults to None.
        Q_Branch (bool, optional): Default to False. This parameter only
                affects linear/spherical tops. When True, the perpendicular
                band will be computed (it has a Q-branch). When False, then the
                parallel band will be computed (without Q-branch).
        A_default (float, optional): Default value of A when constant. Defaults to 0.02.
        B_default (float, optional): Default value of B. Defaults to 0.02.
        C_default (float, optional): Default value of C. Defaults to 0.02
        T_default (float, optional): Defaults to 10
        delta_A_default (float, optional): Defaults to 0.001
        delta_B_default (float, optional): Defaults to 0.001
        delta_C_default (float, optional): Defaults to 0.001
        Jlimit (float, optional): Defaults to 200. How many values of J to consider.
        lambda0 (float, optional): Defaults to 6614. Center wavelength (in A) of profile
        save_fig (bool, optional): Defaults to False.
        anim (bool, optional): Create GIF or not. Defaults to False.
        SymmetryType (str, optional): Type of symmetry of molecule. Options Oblate, Prolate
            or Spherical. Defaults to Prolate.
    """
    # Declare constants

    Sightline = 'ParameterSurvey'
    Q_scale_init = 1

    # Dictionary to save simulations per parameter varying.
    param_simulations = {'A': [], 'B': [], 'C': [], 'T': [],
                         'deltaA': [], 'deltaB': [], 'deltaC': []}

    # Dictionary to save parameters used for simulation
    parameters = {'A': [], 'B': [], 'C': [], 'T': [],
                  'deltaA': [], 'deltaB': [], 'deltaC': []}

    # Dictionary with parameter space
    param_values = {'A': A_values, 'B': B_values, 'C': C_values, 'T': T_values,
                    'deltaA': delta_A_values, 'deltaB': delta_B_values,
                    'deltaC': delta_C_values}

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

                    # Set default parameters values
                    A = A_default
                    B = B_default
                    C = C_default
                    T = T_default
                    delta_A = delta_A_default
                    delta_B = delta_B_default
                    delta_C = delta_C_default

                    # Set changing parameter value
                    if param == 'A':
                        A = A_values[i]

                    elif param == 'B':
                        B = B_values[i]

                    elif param == 'C':
                        C = C_values[i]

                    elif param == 'T':
                        T = T_values[i]

                    elif param == 'deltaA':
                        delta_A = delta_A_values[i]

                    elif param == 'deltaB':
                        delta_B = delta_B_values[i]

                    elif param == 'deltaC':
                        delta_C = delta_C_values[i]

                    if SymmetryType == 'Prolate':
                        C = B
                        delta_C = delta_B

                    elif SymmetryType == 'Oblate':
                        A = B
                        delta_A = delta_B

                    elif SymmetryType == 'Spherical':
                        A = C = B
                        delta_A = delta_A = delta_B

                    # Run simulation
                    build = sim(A=A, B=B, C=C, Delta_A=delta_A*A, Delta_B=delta_B*B,
                                Delta_C=delta_C*C, Trot=T, Jlimit=Jlimit, Target=Sightline,
                                lambda0=lambda0, Q_Branch=Q_Branch, Q_scale=Q_scale)

                    # Obtain wavenumbers.
                    x_vals = WavelengthToWavenumber(np.asarray(build[0])) - \
                        WavelengthToWavenumber(lambda0)
                    y_vals = build[1]

                    # Save results and parameters.
                    param_simulations[param].append([x_vals, y_vals])
                    parameters[param].append({'A': A, 'B': B, 'C': C, 'T': T, 'deltaA': delta_A,
                                              'deltaB': delta_B, 'deltaC': delta_C})

            # Save simulation to future reference.
            if save:
                with open(f"{path}/_Changing{param}_.bin", "wb") as output:
                    pickle.dump([param_simulations[param], parameters[param]], output)

            # Generates GIF and plot.
            plot_simulations(param_simulations[param], param, parameters[param],
                             lambda0, path=path, save=save_fig, SymmetryType=SymmetryType)

            if anim:
                animation(param_simulations[param], param, parameters[param],
                          lambda0, path=path)


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
    """Perform a two-variable survey on spherical/linear tops. All combinations.

    This method works for spherical/linear tops. It performs six two-variable surveys
    corresponding to the combinations of varying B, deltaB, and T. The values of each
    variable must be provided. Default values of variables to use when not varying are
    optional.

    Args:
        path (str): Directory for saving simulations and figures.
        B_values (1darray): Parameter space values of rotational constant.
        T_values (1darray): Values of rotational temperature.
        delta_values (1darray): Values of the change in the rotational constant
            to the upper level.
        Q_Branch (bool, optional): Specify if consider the Q_branch or not.
            Defaults to True.
        B_default (float, optional): Default value of B when constant. Defaults
            to 0.02.
        T_default (float, optional): Default value of T. Defaults to 10.
        delta_default (float, optional): Default value of deltaB. Defaults to 0.
        save (bool, optional): If True, the survey simulations are saved on the
            provided path. Defaults to False.
        load (bool, optional): If True, the simulations are loaded instead of
            computed. A previous executions whit save = True is necessary.
            Defaults to False.
        plot_wavenumber (bool, optional): If True, the plots are created with
            wavenumbers instead of wavelength. Defaults to True.
        Q_scale (float, optional): Scale factor for the Q_branch. Defaults to 1.

    Returns:
        None.
    """
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

    # Save simulation for future reference.
    if save:
        with open(f"{path}/{SymmetryType}_2D_survey.bin", "wb") as output:
            pickle.dump([param_simulations, parameters], output)

    # Create visualizations of general two-paramter survey
    plot_grid(param_simulations, parameters, param_values, lambda0,
              path, SymmetryType, plot_wavenumber)


def plot_grid(simulations, parameters, param_values,
              lambda0, path, SymmetryType, plot_wavenumber):
    """Create grid plots of two-variable surveys.

    The dictionary structures and simulations needed to generate the grid plots are
    those generated by the two_variable_survey function. It can be used in past-saved
    surveys or new ones.

    Args:
        simulations (dict): Dictionary containing the simulations.
        parameters (dict): Dictionary of parameters used in each simulation.
        param_values (dict): Dictionary with parameter space used in the simulation.
        lambda0 (float): Center wavelength of interest (same as simulation).
        path (str): Directory for saving plots.
        SymmetryType (str): Symmetry type of simulations.
        plot_wavenumber (bool): If true, the plots will be generated in wavenumber space.

    Returns:
        None.
    """
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


def single_two_variable_survey(path, values1, values2, name1, name2, Q_Branch=True, A_default=0.02,
                               B_default=0.02, C_default=0.02, T_default=10, delta_A_default=0,
                               delta_B_default=0, delta_C_default=0, save=False,
                               load=False, plot_wavenumber=True, Q_scale=1, plot=True,
                               SymmetryType='Spherical'):
    """Perform a two-variable survey but only with two variables.

    Similar to the two_variable_survey function, but now for a single pair of parameters.
    It can perform surveys in prolate and oblate symmetric tops, not only with spherical/linear
    tops. Depending on the symmetry type, you can vary the parameters A, B, C, deltaA, deltaB,
    deltaC, and T. This method also creates a single grid plot of the survey.

    Args:
        path (str): Directory to save plot and simulation (if specified).
        values1 (1darray): Array with values of the first varying parameter.
        values2 (1darray): Array with values of the second varying parameter.
        name1 (str): Name of the first varying parameter. Options: A, B, C, T, delta_A,
            delta_B, delta_C.
        name2 (str): Name of the second varying parameter.
        Q_Branch (bool, optional): If True, the Q_branch is considered in simulations
            (only affects in Spherical tops). Defaults to True.
        A_default (float, optional): Default value for rotational constant A. Defaults to 0.02.
        B_default (float, optional): Default value for rotational constant B. Defaults to 0.02.
        C_default (float, optional): Default value for rotational constant C. Defaults to 0.02.
        T_default (float, optional): Default value for temperature T. Defaults to 10.
        delta_A_default (float, optional): Default value for the change in the rotational
            constant A, delta_A. Defaults to 0.
        delta_B_default (float, optional): Default value for delta_B. Defaults to 0.
        delta_C_default (float, optional): Default value for delta_C. Defaults to 0.
        save (bool, optional): If true, the resulting simulation is saved in the directory.
            Defaults to False.
        load (bool, optional): If True, the simulations are loaded instead of computed.
            A previous executions whit save = True is necessary. Defaults to False.
        plot_wavenumber (bool, optional): If True, the plots are created with wavenumbers
            instead of wavelength. Defaults to True.
        Q_scale (float, optional): Scale factor for the Q_branch. Defaults to 1.
        plot (bool, optional): If true, the grid plot of the survey is created. Defaults to True.
        SymmetryType (str, optional): Symmetry type to consider. It must agree with A, B,
            and C values. Defaults to 'Spherical'.

    Returns:
        None.
    """
    # Declare constants
    Jlimit = 200
    lambda0 = 6614
    a = WavelengthToWavenumber(lambda0)
    Sightline = 'ParameterSurvey'
    name = name1 + name2

    # Dictionary to save simulations
    param_simulations = {}

    # Dictionary to save parameters used in simulations.
    parameters = {}

    # Load pre-calculated simulations
    if load:
        with open(f"{path}/{name}_2D_single_survey.bin", "rb") as data:
            # Extract data
            build = pickle.load(data)
            param_simulations = build[0]
            parameters = build[1]

    # Run simulation.
    else:

        # Simulations
        for value1 in values1:
            param_simulations[value1] = {}
            parameters[value1] = {}
            for value2 in values2:

                # Assign default values.
                params = {'A': A_default, 'B': B_default, 'C': C_default, 'T': T_default,
                          'delta_A': delta_A_default, 'delta_B': delta_B_default,
                          'delta_C': delta_C_default}

                # Assign variable values.
                params[name1] = value1
                params[name2] = value2

                if SymmetryType == 'Prolate':
                    params['C'] = params['B']
                    params['delta_C'] = params['delta_B']

                elif SymmetryType == 'Oblate':
                    params['A'] = params['B']
                    params['delta_A'] = params['delta_B']

                elif SymmetryType == 'Spherical':
                    params['A'] = params['B']
                    params['delta_A'] = params['delta_B']
                    params['C'] = params['B']
                    params['delta_C'] = params['delta_B']

                build = sim(A=params['A'], B=params['B'], C=params['C'],
                            Delta_A=params['delta_A']*params['A'],
                            Delta_B=params['delta_B']*params['B'],
                            Delta_C=params['delta_C']*params['C'],
                            Trot=params['T'], Jlimit=Jlimit, Target=Sightline, lambda0=lambda0,
                            Q_Branch=Q_Branch, Q_scale=Q_scale)

                # Obtain wavenumbers.
                x_vals = WavelengthToWavenumber(np.asarray(build[0]))
                y_vals = build[1]

                # Save results and parameters.
                param_simulations[value1][value2] = [x_vals, y_vals]

    # Save simulation to future reference.
    if save:
        with open(f"{path}/{name}_2D_single_survey.bin", "wb") as output:
            pickle.dump([param_simulations, parameters], output)

    if plot:
        rows = len(values1)
        cols = len(values2)

        fig, axs = plt.subplots(nrows=rows, ncols=cols, sharex=True,
                                sharey=True, figsize=(rows*2, cols*2))

        for i, value1 in enumerate(values1):
            for j, value2 in enumerate(values2):

                x = param_simulations[value1][value2][0]
                y = param_simulations[value1][value2][1]

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
                    axs[i, j].set_title(f'{name2} = {value2:.4f}')

                if j == cols-1:
                    axs[i, j].set_ylabel(f'{name1} = {value1:.4f}')
                    axs[i, j].yaxis.set_label_position("right")

                if j == 0:
                    axs[i, j].set_ylabel('Normalized Intensity')

        plt.suptitle(f'2D parameter survey with {name1} and {name2}.\n' +
                     f'A = {A_default}, B = {B_default}, C = {C_default}, T = {T_default}\n' +
                     f'deltaA = {delta_A_default}, deltaB = {delta_B_default}, ' +
                     f'deltaC = {delta_C_default}')
        plt.tight_layout()
        plt.savefig(f"{path}/{name}_2D_single_survey.pdf", dpi=300)
        plt.show()


if __name__ == "__main__":

    # 2D parameter survey for spherical tops. Turn to True to run.
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

    # One variable survey for prolate and oblate symmetries.
    if False:

        # Define parameter ranges.
        delta_A_values = np.linspace(0, 0.5, 10)
        one_variable_survey(path='Parameter_survey', delta_A_values=delta_A_values,
                            A_default=0.1, load=True, save_fig=True, anim=True,
                            SymmetryType='Prolate')

        # Define parameter ranges.
        delta_C_values = np.linspace(0, 0.5, 10)
        one_variable_survey(path='Parameter_survey', delta_C_values=delta_C_values,
                            C_default=0.005, load=True, save_fig=True, anim=True,
                            SymmetryType='Oblate')

    # Single 2D variable survey for prolate and oblate tops.
    if True:

        # Define parameter ranges.
        B_values = [0.001, 0.005, 0.01, 0.05]
        A_values = [0.1, 0.5, 1, 5]
        single_two_variable_survey(path='Parameter_survey/single', values1=A_values,
                                   values2=B_values, name1='A', name2='B', delta_A_default=0.005,
                                   delta_B_default=0.005, delta_C_default=0.005, load=True,
                                   T_default=70, Q_scale=0.4, SymmetryType='Prolate')

        # Define parameter ranges.
        C_values = [0.001, 0.005, 0.01, 0.05]
        B_values = [0.1, 0.5, 1, 5]
        single_two_variable_survey(path='Parameter_survey/single', values1=A_values,
                                   values2=B_values, name1='A', name2='B', delta_A_default=0.005,
                                   delta_B_default=0.005, delta_C_default=0.005, save=True,
                                   T_default=70, Q_scale=0.4, SymmetryType='Oblate')
