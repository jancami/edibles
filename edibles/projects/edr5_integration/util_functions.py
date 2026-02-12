import numpy as np

def setting_dependent_crop(spec, wave):
    crop_lim_dict = {346: [10, 10], 437: [13, 7], 564: [19, 4], 860: [20, 0]}
    crop_limits = np.array(crop_lim_dict[wave])
    cl_ang = [np.nanmin(spec[0]) + crop_limits[0], np.nanmax(spec[0]) - crop_limits[1]]

    return cl_ang


def crop_spectrum(array_in: np.array, x_min: float, x_max: float) -> np.array:
    """
    Returns a spectrum interval for x_min < wave < x_max.

    Parameters
    ----------
    array_in : np.array([wave, flux, additional_columns])
        Input spectrum.
    x_min : float
        Minimum wave coordinate of slice.
    x_max : float
        Maximum wave coordinate of slice.

    Returns
    -------
    np.array([wave, flux, additional_columns])
        Spectrum slice
    """
    if x_min > x_max:
        raise ValueError('Slice_spectrum error: x_min is larger than x_max!')

    b1 = array_in[0] < x_max  # boolean array of wave values smaller than x_max
    b2 = x_min < array_in[0]  # boolean array of wave values larger than x_min
    bool_array = np.logical_and(b1, b2)  # boolean array of wave values larger than x_min and smaller than x_max

    return array_in[:, bool_array]
