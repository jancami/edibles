import numpy as np
def data_cleaning(x, data, minrange, maxrange):
    """Filters data arrays to include only points within a specific range.

    Clips the input data to retain only coordinates where the independent
    variable falls between the provided minimum and maximum bounds.

    Args:
        x: Independent variable data points (wavelengths) as a NumPy array.
        data: Observed dependent variable data points (flux).
        minrange: Lower bound for the wavelength filter.
        maxrange: Upper bound for the wavelength filter.

    Returns:
        A tuple containing:
            - The filtered independent variable array.
            - The corresponding filtered dependent variable array.
    """
    filtered_pairs = [(xi, yi) for xi, yi in zip(x, data) if int(minrange) <= xi <= int(maxrange)]

    if filtered_pairs:
        filtered_x, filtered_data = map(np.array, zip(*filtered_pairs))
    else:
        filtered_x = np.array([])
        filtered_data = np.array([])

    x = filtered_x
    data = filtered_data

    return x, data
