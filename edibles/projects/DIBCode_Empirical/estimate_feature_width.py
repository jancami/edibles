import numpy as np
def estimate_feature_width(x, residuals, center_idx):
    """Estimates the width of a spectral feature based on its half-maximum.

    Traces the feature outwards from its peak until the absolute amplitude drops
    to half of the peak value (FWHM), converts it to a Gaussian standard
    deviation sigma, and clips the result to physically reasonable bounds.

    Args:
        x: Independent variable data points (wavelengths) as a NumPy array.
        residuals: Array of residual values from which to estimate the feature
          profile.
        center_idx: Index in the arrays corresponding to the feature peak.

    Returns:
        The estimated Gaussian width (sigma) parameter as a float.
    """
    peak_val = abs(residuals[center_idx])
    half_max = peak_val / 2.0

    # Search left
    left_idx = center_idx
    while left_idx > 0 and abs(residuals[left_idx]) > half_max:
        left_idx -= 1

    # Search right
    right_idx = center_idx
    while right_idx < len(residuals) - 1 and abs(residuals[right_idx]) > half_max:
        right_idx += 1

    # FWHM estimate
    fwhm = x[right_idx] - x[left_idx]

    # Convert to sigma (for our Gaussian parameterization)
    width_guess = fwhm / 1.665 if fwhm > 0 else np.ptp(x) * 0.01

    # Sanity check
    x_step = np.median(np.diff(x))
    min_width = 2 * x_step  # At least 2 data points wide
    max_width = np.ptp(x) * 0.2  # At most 20% of range

    width_guess = np.clip(width_guess, min_width, max_width)

    return width_guess
