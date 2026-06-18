def calculate_reduced_chi_square(result, n_data):
    """Calculates the reduced chi-square statistic for an lmfit result.

    Divides the total chi-square by the remaining degrees of freedom (number of
    data points minus the number of actively varying parameters).

    Args:
        result: The result object returned after executing an lmfit minimization
          fit.
        n_data: Number of data points processed during the fit.

    Returns:
        The calculated reduced chi-square value as a float.
    """
    return result.chisqr / (n_data - result.nvarys)
