import numpy as np
def composite_model(x, params):
    """Computes a flexible model combining Gaussians, Lorentzians, and Chebyshev continuum.

    The model extracts numbered components from the parameters based on naming conventions:
      - Gaussian i: 'g{i}_amp', 'g{i}_cen', 'g{i}_wid'
      - Lorentzian i: 'l{i}_amp', 'l{i}_cen', 'l{i}_wid'
      - Chebyshev: 'c0', 'c1', 'c2', ...

    Note:
        For absorption features, amplitude parameters should be negative.

    Args:
        x: Independent variable data points as a NumPy array.
        params: An lmfit Parameters object containing the model coefficients.

    Returns:
        The calculated y-values as a NumPy array of the same shape as x.
    """
    y = np.zeros_like(x)

    # Add Gaussian components
    i = 0
    while f'g{i}_amp' in params:
        amp = params[f'g{i}_amp'].value
        cen = params[f'g{i}_cen'].value
        wid = params[f'g{i}_wid'].value
        y += amp * np.exp(-((x - cen) / wid) ** 2)
        # print("These are the Gaussian Model function parameters", amp, cen, wid, "of iteration", i)
        i += 1

    # Add Lorentzian components
    i = 0
    while f'l{i}_amp' in params:
        amp = params[f'l{i}_amp'].value
        cen = params[f'l{i}_cen'].value
        wid = params[f'l{i}_wid'].value
        y += amp * wid ** 2 / ((x - cen) ** 2 + wid ** 2)
        # print("These are the Lorentzian Model function parameters", amp, cen, wid, "of iteration", i)
        i += 1

    # Add Chebyshev polynomial
    cheb_coeffs = []
    i = 0
    while f'c{i}' in params:
        cheb_coeffs.append(params[f'c{i}'].value)
        i += 1

    if cheb_coeffs:
        # Normalize x to [-1, 1] for numerical stability
        x_norm = 2 * (x - x.min()) / (x.max() - x.min()) - 1
        y += chebyshev.chebval(x_norm, cheb_coeffs)
    return y
