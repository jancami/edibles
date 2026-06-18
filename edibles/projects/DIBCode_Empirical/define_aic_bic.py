import numpy as np
def calculate_aic(chi_square, n_params, n_data):
    """Calculates the Akaike Information Criterion (AIC).

    A lower AIC value indicates a better-fitting model while penalizing
    unnecessary parameter growth. Note that `n_data` is accepted for signature
    consistency but is not used in the standard linear AIC definition here.

    Args:
        chi_square: The sum of squared residuals from the fit.
        n_params: Number of free parameters in the model.
        n_data: Number of data points used in the fit.

    Returns:
        The calculated AIC value as a float.
    """
    return chi_square + 2 * n_params


def calculate_bic(chi_square, n_params, n_data):
    """Calculates the Bayesian Information Criterion (BIC).

    A lower BIC indicates a better model. BIC penalizes model complexity more
    harshly than AIC, especially for larger data sample sizes.

    Args:
        chi_square: The sum of squared residuals from the fit.
        n_params: Number of free parameters in the model.
        n_data: Number of data points used in the fit.

    Returns:
        The calculated BIC value as a float.
    """
    return chi_square + n_params * np.log(n_data)
