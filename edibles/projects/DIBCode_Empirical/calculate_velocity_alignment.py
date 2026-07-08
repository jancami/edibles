import numpy as np
from edibles.projects.DIBCode_Empirical.get_shifted_data import get_shifted_data
from edibles.projects.DIBCode_Empirical.calculate_correlation import calculate_correlation
def calculate_velocity_alignment(x_ref, y_obs, model_template, v_grid, mask):
    """
    Calculates correlations across the grid and returns best fit data
    Args:
        x_ref: Initial independent variable data points (wavelengths) as a NumPy array.
        y_obs: Observed flux values from the sightline.
        model_template: Sum of the observed continuum and the initial model fit.
        v_grid: List of possible velocity values to run correlations over.
        mask: Mask of the model fitting region.
    """
    correlations = []
    C_KMS = 299792.458
    for v in v_grid:
        y_shifted = get_shifted_data(-v, x_ref, model_template)

        # Replace out-of-bounds data (prevent extreme extrapolation)

        mask_nan=[]
        for index in range(len(y_shifted)):
            mask_nan.append(not np.isnan(y_shifted[index]))



        # Assumes calculate_correlation is defined in the global scope
        r = calculate_correlation(y_shifted[mask_nan], y_obs[mask_nan])
        correlations.append(r)

    correlations = np.array(correlations)
    best_idx = np.argmax(correlations)
    best_v = v_grid[best_idx]
    max_r = correlations[best_idx]
    best_y_data = get_shifted_data(best_v, x_ref, y_obs)

    return best_v, max_r, correlations, best_y_data
