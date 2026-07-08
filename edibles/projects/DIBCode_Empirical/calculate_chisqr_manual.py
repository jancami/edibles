import numpy as np
from edibles.projects.DIBCode_Empirical.composite_model import composite_model
def calculate_chisqr_manual(params, x, data, uncertainties):
    model = composite_model(x, params)
    if uncertainties is not None:
        chisqr = np.sum(((data - model) / uncertainties) ** 2)
    else:
        chisqr = np.sum((data - model) ** 2)
    return chisqr


"""Manually calculates the chi-square value for verification.
  Args:
        params: An lmfit Parameters object containing the model coefficients.
        x: Independent variable data points as a NumPy array.
        data: Observed dependent variable data points.
        uncertainties: Array of uncertainties (sigma) for each data point, or None if unweighted.

           Returns:
        The calculated chi-square value as a float. """
