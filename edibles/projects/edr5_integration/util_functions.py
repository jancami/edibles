import numpy as np
from scipy.interpolate import interp1d

def setting_dependent_crop(spec, setting):
    crop_lim_dict = {346: [10, 10], 437: [13, 7], 564: [19, 4], 860: [20, 0]}
    crop_limits = np.array(crop_lim_dict[setting])
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


def normalize_spectrum_linear(spectrum: np.array, cont_1: np.array, cont_2: np.array,
                              additional_normalized_columns: list = None) -> np.array:
    """
    Normalizes an absorption spectrum using a linear continuum model defined by two points con_1 and cont_2.

    Parameters
    ----------
    spectrum : np.array([wave, flux, additional_columns])
        Un-normalized spectrum.
    cont_1 : np.array([x, y])
        First continuum point.
    cont_2 : np.array([x, y])
        Second continuum point.
    additional_normalized_columns : list
        Array of column indices. Specifies additional columns to be normalized, using the same continuum points.
        Mind that the additional columns start with the index 2.
        E.g. [2, 3, 5]

    Returns
    -------
    np.array([wave, flux, additional_columns])
        Normalized spectrum.
    """
    wave = spectrum[0]
    flux = spectrum[1]
    additional_columns = spectrum[2:]
    sx1r = cont_1[0]
    sy1r = cont_1[1]
    sx2r = cont_2[0]
    sy2r = cont_2[1]

    f = interp1d([sx1r, sx2r], [sy1r, sy2r], kind='linear', fill_value="extrapolate")

    flux = flux / f(wave)

    if additional_normalized_columns is not None:
        for i in additional_normalized_columns:
            k = i - 2  # shift index from spectrum numbering to numbering in additional_columns
            additional_columns[k] = additional_columns[k] / f(wave)

    if len(additional_columns) > 0:
        out_spec = np.concatenate((np.array([wave, flux]), additional_columns))
    else:
        out_spec = np.array([wave, flux])

    return out_spec
