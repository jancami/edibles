import numpy as np
from edibles.projects.DIBCode_Empirical.estimate_feature_width import estimate_feature_width
def initialize_lorentzian(x, residuals, prefix='l0_'):
    """Initializes a Lorentzian profile at the location of strongest absorption.

    Args:
        x: Independent variable data points (wavelengths) as a NumPy array.
        residuals: Array of residual values (observed - continuum) used to find
          absorption dips.
        prefix: String identifier prefix to prepend to the lmfit parameter
          dictionary keys.

    Returns:
        A nested dictionary mapping component names to their initial values,
        minimum limits, and maximum limits.
    """
    idx_min = np.argmin(residuals)
    center_guess = x[idx_min]
    amp_guess = residuals[idx_min]  # Should be negative for absorption

    # Estimate width (Lorentzian HWHM)
    width_guess = estimate_feature_width(x, residuals, idx_min)

    print(f"  Initializing {prefix[:-1]}: amp={amp_guess:.2f}, cen={center_guess:.8f}, wid={width_guess:.8f}")

    return {
        f'{prefix}amp': {'value': amp_guess, 'min': -np.inf, 'max': 0},  # ABSORPTION: amp <= 0
        f'{prefix}cen': {'value': center_guess, 'min': x.min(), 'max': x.max()},
        f'{prefix}wid': {'value': width_guess, 'min': np.median(np.diff(x)), 'max': np.ptp(x) * 0.3}
    }
