import numpy as np
def get_shifted_data(v, x_ref, y_obs):
    """
    Interpolates data to a new velocity shift.
    Args:
        v: Velocity of the observed sightline.
        x_ref: Initial independent variable data points (wavelengths) as a NumPy array.
        y_obs: Observed flux values.

    Returns:
        Array of interpolated data, shifted by wavelength x_shifted.
    """
    C_KMS = 299792.458
    x_shifted = x_ref * (1 + v / C_KMS)
    return np.interp(x_ref, x_shifted, y_obs)
