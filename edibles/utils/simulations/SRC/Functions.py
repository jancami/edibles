from lmfit.models import LinearModel
import numpy as np
import pandas as pd
"""
Python script to hold all of the little functions that get used throughout this project"
"""


def Signal_Noise_Calculator(x_vals, y_vals):
    """
    Calculates the signal to noise ratio (S/N) of a inputted spectrum by fitting a linear model to the spectrum.
    Args:
        x_vals (1darray): Wavelength values.
        y_vals (1darray): Flux values

    Returns:
        SN (float): The S/N value for the input spectrum.
        y_fit (1darray): Array consisting of flux values which correspond to the best linear fit through the input spectrum.

    """
    # Fit a linear model through the inpit x_vals and y_vals
    lin_fit = LinearModel()
    lin_fit_pars = lin_fit.guess(y_vals, x=x_vals)
    lin_fit_out = lin_fit.fit(y_vals, lin_fit_pars, x=x_vals)
    y_fit = lin_fit_out.best_fit
    # Find the residuals of the best fit.
    resid = lin_fit_out.residual
    # Calculate the S/N using the Root Mean Squared Error(RMSE)
    resid_square = [r**2 for r in resid]
    resid_error = np.sqrt(np.sum(resid_square)/len(resid_square))
    SN = np.mean(y_vals)/resid_error

    return (SN, y_fit)


def weighted_average(values, uncertainties):
    """
    Calculates the weighted average, as well as the uncertainty of this average, of an input set of values
    Args:
        values (1darray): Values to have the average taken of.
        uncertanties (1darray): The uncertainty of the measurement of the values passed
    Returns:
        weighted_average (1d array): Weighted average of the values
        weighted_average_err (1d array): Uncertainty of the weighted average.

    """
    if isinstance(values, pd.Series):
        weights = 1/(uncertainties**2)
        values_times_weights = values*weights
        weighted_average = np.sum(values_times_weights)/np.sum(weights)
        weighted_average_err = 1/(np.sqrt(np.sum(weights)))
    else:
        weights = [1/sigma**2 for sigma in uncertainties]
        values_times_weights = [values[i]*weights[i] for i in range(len(values))]
        weighted_average = np.sum(values_times_weights)/np.sum(weights)
        weighted_average_err = 1/(np.sqrt(np.sum(weights)))

    return(weighted_average, weighted_average_err)


def calculate_average(values, uncertainties):
    """
    Calculates the average, as well as the uncertainty of this average, of an input set of values.
    Args:
        values (1darray or Series): Values to have the average taken of.
        uncertanties (1darray or Series): The uncertainty of the measurement of the values passed
    Returns:
        average (1darray or Series): Average of the values
        uncertainty (1darray or Series): Uncertainty of the average.

    """
    # If/Es = lse statement to determine if the values input were stored in arrays or series. The inner working of both components are the same, but the manipulation is different due to different types. Hope is to eventually transition to only using series in this project.
    if isinstance(values, pd.Series):

        average = values.sum()/values.size
        squares = uncertainties.pow(2)
        uncertainty = (1/uncertainties.size)*np.sqrt(squares.sum())

    else:
        average = np.sum(values)/len(values)
        squares = uncertainties**2
        uncertainty = (1/len(uncertainties))*np.sqrt(np.sum(squares))

    return(average, uncertainty)


def Sort_Points(x_values, y_values, error):
    """
    Function which combines three arrays and sorts by ascending values of the first.
    Args:
        x_values (1darray): Array of x-values. This array is the one which sorting will be based. Each index of this array is related to the same index value of the other two.
        y_values (1darray): Array of y-values. Each index of this array is related to the same index value of the other two.
        error (1darray): Array of error values. Each index of this array is related to the same index value of the other two.
    Returns:
        x_val (1darray): x_values array sorted in ascending order.
        y_val (1darray): y_values array sorted based on new index order of x_val.
        errors (1darray): error array sorted based on new index order of x_val.
    """
    to_return = []
    # If/Else statement to confirm that the arrays are the same length. Print error message if not.
    if len(x_values) == len(y_values) == len(error):

        for j in range(len(x_values)):
            file_point = (((x_values[j])), (y_values[j]), (error[j]))
            to_return.append(file_point)
            to_return.sort()

        x_val = [x[0] for x in to_return]
        y_val = [y[1] for y in to_return]
        errors = [y[2] for y in to_return]
    else:
        print('the given arrays are of different lengths')
        x_val = False
        y_val = False
        errors = False

    return(x_val, y_val, errors)


def Sort_Points_4(x_values, y_values, xerr, yerr):
    """
    Function which combines three arrays and sorts by ascending values of the first.
    Args:
        x_values (1darray): Array of x-values. This array is the one which sorting will be based. Each index of this array is related to the same index value of the other three.
        y_values (1darray): Array of y-values. Each index of this array is related to the same index value of the other three.
        xerr(1darray): Array of error values which correspond to x_values. Each index of this array is related to the same index value of the other three.

        yerr(1darray): Array of error values which correspond to y_values. Each index of this array is related to the same index value of the other three.
    Returns:
        x_val (1darray): x_values array sorted in ascending order.
        y_val (1darray): y_values array sorted based on new index order of x_val.
        xer (1darray): xerr array sorted based on new index order of x_val.
        yer (1darray): yerr array sorted based on new index order of x_val.
    """
    to_return = []
    # If/Else statement to confirm that the arrays are the same length. Print error message if not.
    if len(x_values) == len(y_values) == len(xerr) == len(yerr):

        for j in range(len(x_values)):
            file_point = ((x_values[j]), (y_values[j]), (xerr[j]), (yerr[j]))
            to_return.append(file_point)
            to_return.sort()

        x_val = [x[0] for x in to_return]
        y_val = [y[1] for y in to_return]
        xer = [y[2] for y in to_return]
        yer = [y[3] for y in to_return]

    else:
        print('the given arrays are of different lengths')
        x_val = False
        y_val = False
        xer = False
        yer = False

    return(x_val, y_val, xer, yer)


def Sort_Points_2(x_values, y_values):
    """
    Function which combines three arrays and sorts by ascending values of the first.
    Args:
        x_values (1darray): Array of x-values. This array is the one which sorting will be based. Each index of this array is related to the same index value of the other.
        y_values (1darray): Array of y-values. Each index of this array is related to the same index value of the other.
    Returns:
        x_val (1darray): x_values array sorted in ascending order.
        y_val (1darray): y_values array sorted based on new index order of x_val.
    """
    to_return = []
    # If/Else statement to confirm that the arrays are the same length. Print error message if not.
    if len(x_values) == len(y_values):

        for j in range(len(x_values)):
            file_point = (((x_values[j])), (y_values[j]))
            to_return.append(file_point)
            to_return.sort()

        x_val = [x[0] for x in to_return]
        y_val = [y[1] for y in to_return]

    else:
        print('the given arrays are of different lengths')
        x_val = False
        y_val = False

    return(x_val, y_val)


def InverseFit(x_vals, y_vals):

    import numpy as np
    import matplotlib.pyplot as plt
    from lmfit.models import PowerLawModel

    fit = PowerLawModel()
    fit_pars = fit.guess(y_vals, x=x_vals)
    fit_out = fit.fit(y_vals, fit_pars, x=x_vals)

    A = fit_out.params['amplitude'].value
    E = fit_out.params['exponent'].value
    print(A, E)
    return(A, E)


def WavelengthToWavenumber(values):
    with np.errstate(divide='ignore'):
        wavenumbers = 1/(values*(10**-8))
    return(wavenumbers)


def WavenumberToWavelength(values):
    with np.errstate(divide='ignore'):
        wavenumbers = 10**8/(values)
    return(wavenumbers)


def LAD_Fit(x_vals, y_vals):

    from scipy.optimize import minimize
    import numpy as np
    import matplotlib.pyplot as plt

    def fit(X, params):
        return X.dot(params)

    def cost_function(params, X, y):
        return np.sum(np.abs(y-fit(X, params)))

    X = np.asarray([np.ones((len(x_vals),)), x_vals]).T
    y = y_vals

    x0 = [0.01, 1]

    output = minimize(cost_function, x0, args=(X, y))
    y_hat = fit(X, output.x)

    return(y_hat, output.x)
