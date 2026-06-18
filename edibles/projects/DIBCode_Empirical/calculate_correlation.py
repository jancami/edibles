import numpy as np
def calculate_correlation(y_data, y_model):
    """
    Calculates correlation coefficient.
    Args:
        y_data: List of the data flux.
        y_model: List of the model flux.
    Returns:
        Value of the correlation coefficient.
    """
    mu_d = np.mean(y_data)
    mu_m = np.mean(y_model)
    num = np.sum((y_data - mu_d) * (y_model - mu_m))
    den = np.sqrt(np.sum((y_data - mu_d) ** 2) * np.sum((y_model - mu_m) ** 2))
    return num / den if den != 0 else 0
