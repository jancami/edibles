import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt


def generate_continuum(data_tuple, delta_v=1.0, n_piece=None):
    '''
    This function fits a continuum to data separated into n sections
    where the x and y-values are the median of each section using a cubic spline

    INPUT:        [Angstroms]
    data_tuple:   [tuple: (ndarray, ndarray)] (wavelength grid, flux values)
    delta_v:      [float]                     desired resolution of continuum (in m/s)
    n_piece:      [int, default=4]            number of sections to split data into

    OUTPUT:  [Angstroms]
    x_cont:  [ndarray] wavelength grid points for fit continuum
    y_cont:  [ndarray] flux value points for fit continuum
    
    Graphic:

    X------|-:    |    :-|-----X
    |      |  :   |   :  |     |
    |      |   :  X  :   |     |
    |      |    :-|-:    |     |
    |______|______|______|_____|

    <----piece---><---piece---->
    <-sec-><-sec-><-sec-><-sec->

    '''
    (x, y) = data_tuple
    
       # check n_piece param
    if n_piece is None: n_piece = 2

    # split x & y arrays into n_piece*2 sections
    x_sections = np.array_split(x, n_piece*2)
    y_sections = np.array_split(y, n_piece*2)

    # initialize list of points to spline fit
    x_points = [np.min(x)]
    y_points = [np.median(y_sections[0])]


    # loop through every other section (1, 3, 5...)
    # make n_piece+1 points to fit a spline through
    # create a spline point on edge of each piece
    for i in range(1, len(x_sections), 2):
        # set x_point 
        x_point = np.max(x_sections[i])

        if i == range(len(x_sections))[-1]:
            span = y_sections[i]
            y_point = np.median(span)
            x_point = np.max(x)

        else:
            span = np.append(y_sections[i], y_sections[i+1])
            y_point = np.median(span)

        x_points.append(x_point)
        y_points.append(y_point)

    spline = CubicSpline(x_points, y_points)
    y_spline = spline(x)
    plt.plot(x_points, y_points, 'kx', markersize='8', label='Points')

    # cont_tuple = (x_spline, y_spline)

    return y_spline
