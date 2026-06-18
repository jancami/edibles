import numpy as np
from numpy.polynomial import Chebyshev
def fitting_objective(scale, x_ref, best_y_data, c0, ref_y_profiles, order):
    """
    The objective that the fitting function is trying to optimize
    Args:
        scale: Scaling of the function
        x_ref: Initial independent variable data points (wavelengths) as a NumPy array.
        best_y_data: The best profile for the DIB fit.
        c0: First-order coefficient of the Chebyshev polynomial.
        ref_y_profiles: Reference profiles for each sightline.
        order: Order of the Chebyeshev polynomial.
    Returns:
        Heuristic for the fitting function.
    """
    scaled_total = scale * ref_y_profiles
    poly_target = best_y_data - (c0 + scaled_total)
    poly = Chebyshev.fit(x_ref, poly_target, deg=order)
    y_model = poly(x_ref) + (c0 + scaled_total)
    return np.sum((best_y_data - y_model) ** 2)
