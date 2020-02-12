import numpy as np
from scipy.interpolate import CubicSpline
import astropy.constants as cst


##### Continuum Fitting Related #####
mode_all = ['spline', 'polynomial']


def polynomial_continuum(data_tuple, n_order=3, errweight = True):
    """
    Fits a n-th order polynominal to input data.

    :param data_tuple: input wavelength grid and flux, in the format (wave, flux)
    :type data_tuple: tuple
    :param n_order: maximal order of the polynomial, default to 3
    :type n_order: int
    :param errweight: if true, Poisson err of flux will be used as weights of fitting, default to True
    :type errweight: bool

    :return: poly, an poly1D object to calculate continuum level
    :rtype: object
    """
    (x, y) = data_tuple
    assert len(x) == len(y), 'x and y must have the same length'
    weights = np.ones_like(x)
    if errweight:
        weights = weights / np.sqrt(y)
    coe_fit = np.polyfit(x, y, n_order, w=weights)
    poly = np.poly1d(coe_fit)

    return poly


def spline_continuum(data_tuple, pars=None, n_piece=2):
    """
    Cuts input data into n_piece sections, and creates a cubic spline continuum from
    median x and y of each sections.
    User can define his/her own anchor points using "pars"

    :param data_tuple: input wavelength grid and flux, in the format (wave, flux)
    :type data_tuple: tuple
    :param pars: user specified anchor points (for y values)
    :type pars: ndarray
    :param n_piece: number of sections, default to 2
    :type n_piece: int

    :return: spline, y_points, cubicspline object to calculate continuum, and anchor points in the fitting
    :rtype: object, ndarray
    """

    (x,y) = data_tuple

    # check n_piece then split x & y arrays into n_piece*2 sections
    if n_piece is None: n_piece = 2
    x_sections = np.array_split(x, n_piece * 2)
    y_sections = np.array_split(y, n_piece * 2)

    # making x_points and y_points used in fitting
    # loop among 1,3,5...,
    x_points = [np.min(x)]
    y_points = [np.median(y_sections[0])]

    for i in range(1,len(x_sections),2):
        if i == range(len(x_sections))[-1]:
            span = y_sections[i]
            y_points.append(np.median(span))
            x_points.append(np.max(x))
        else:
            span = np.append(y_sections[i],y_sections[i+1])
            y_points.append(np.median(span))
            x_points.append(np.max(x_sections[i]))

    # otherwise, import y_points from para
    if pars is not None:
        assert ((n_piece+1) == len(pars)), 'Incorrect number of input y_points parameters'
        y_points = pars

    spline = CubicSpline(x_points, y_points)
    #y_spline = spline(x)

    # plt.plot(x_points, y_points, 'kx', markersize='8', label='Points')

    return spline, y_points


def iterate_continuum(data_tuple,
                       mode='spline', n=3,
                       lsigma=1, usigma=2,
                       iterates=30, min_sigma = 0.2,
                       mask = None):
    """
    Fits continuum of input data in an iteration manner.
    The iteration will break if:
    1) reached the max iteration cycles specified in input
    2) reached the min allowed residual (in terms of n-sigma)

    :param data_tuple: input wavelength grid and flux, in the format (wave, flux)
    :type data_tuple: tuple
    :param mode: fitting method to be used, a member of "mode_all", default to spline
    :type mode: str
    :param n: number of sections in spline mode, and maximal order in polynomial mode, default to 3
    :type n: int
    :param lsigma: allowed n-sigma residual above continuum, default to 2
    :type lsigma: float
    :param usigma: allowed n-sigma residual below continuum, default to 1
    :type usigma: float
    :param iterates: allowed iteration cycles, default to 30
    :type iterates: int
    :param min_sigma: allowed n-sigma residual in both direction,
                      used to break iteration if l/usigma is too small, default to 0.1
    :type min_sigma: float
    :param mask: pre-defined mask to be used in the fitting, default to None
    :type mask: ndarray

    :return: y_cont, idx, cubicspline or poly1D object for continuum and index of continuum region
    :rtype: object, ndarray
    """

    # initializing data
    (x,y) = data_tuple    # original data
    x = np.array(x)
    y = np.array(y)
    (y_res, y_nsigma) = np.zeros_like((y, y))   # residuals and n_sigma
    if mask is None: mask = np.zeros_like(y)
    sigma_old = 9999999
    do_flag = True

    # chosing mode
    assert mode in mode_all, 'Fitting method not available!'
    if mode == 'spline':npiece = n
    if mode == 'polynomial':norder = n

    while do_flag:
        idx = [i for i in range(len(x)) if -usigma <= y_nsigma[i] <= lsigma and mask[i] == 0]
        x_fit = np.copy(x[idx])
        y_fit = np.copy(y[idx])

        if mode == 'spline':
            (y_count, y_points) = spline_continuum((x_fit,y_fit), pars=None, n_piece=npiece)

        if mode == 'polynomial':
            y_count = polynomial_continuum((x_fit,y_fit),norder, errweight= True)

        y_res = y_count(x) - y
        y_sigma = np.std(np.abs(y_res))
        y_nsigma = y_res / y_sigma

        # narrow down the fitting region by more restrict criteria
        if y_sigma >= sigma_old:
            usigma = usigma / 2
            lsigma = lsigma / 2
        sigma_old = y_sigma

        # if iteration should be stopped
        if usigma <= min_sigma or lsigma <= min_sigma:
            do_flag = False
            #print('sigma < 0.1')

        iterates = iterates - 1
        if iterates <= 0:
            do_flag = False
            #print('end of iterates')

        #plt.plot(x, y, color='k')
        #plt.plot(x_fit, y_fit, color='red')
        #plt.plot(x, y_count(x), color='blue')
        #plt.xlabel(str(iterates))
        #plt.show()

    return y_count, idx


##### Wavelength shifting, convert between AA and km/s #####
def velocityShift(wave_array, v_shift):
    """
    A simple function to shift spectrum according to input.
    Used for bary-correction or other on-purpose shifts.

    :param wave_array: wavelength grid
    :type wave_array: ndarray
    :param v_shift: velocity, positive means red-shift
    :type v_shift: float

    :return: wave_out, shifted wavelength grid
    :rtype: ndarray
    """
    wave_out = wave_array + (v_shift / cst.c.to("km's").value) * wave_array
    return wave_out


def waveVelocity(x_in, center=None, mode="wave2velocity"):
    """
    Convert between wavelength and velocity grid

    :param x_in: input x grid
    :type x_in: ndarray
    :param center: zero velocity of wavelength position, wavelength center will be used if None
    :type center: float
    :param mode: wave2velocity or velocity2wave
    :type mode: str

    :return: x_out, center, converted x grid and zero velocity position
    :rtype: ndarray, float
    """
    assert mode in ["wave2velocity", "velocity2wave"], "wave2velocity or velocity2wave"
    wave_in = np.array(x_in)
    ### wave to velocity
    if mode == "wave2velocity":
        if center is None:
            center = np.mean(wave_in)
            print("Zero velocity set to the center of input wavelength")

        x_out = (x_in - center) / center * cst.c.to('km/s').value
    ### velocity to wave
    if mode == "velocity2wave":
        assert center is not None, "Specific center position"
        x_out = (1. + x_in / cst.c.to('km/s').value) * center

    return x_out, center


##### User Interaction and Other Small U#####
def go_message(message):
    """
    Displays a message and collect user's yes-or-no response

    :param message: message to user
    :type message: str

    :return: yes or no decision
    :rtype: bool
    """
    go = ''
    while go.lower() not in ["y", "n"]:
        go = input(message)

    if go.lower() == "y":
        return True
    else:
        return False


def within_boundaries(x, boundary_left, boundary_right):
    """
    Tests if input x is between the left and right boundaries (boundaries NOT included)

    :param x: input x to be examined
    :type x: ndarray
    :param boundary_left: list of left boundaries
    :type boundary_left: ndarray
    :param boundary_right: list of right boundaries, same length as boundary_left
    :type boundary_right: ndarray

    :return: within_all, result for each element in input x, one if within
    :rtype: list
    """
    try:len(x)
    except:x = [x]
    try:len(boundary_left)
    except: boundary_left = [boundary_left]
    try:len(boundary_right)
    except: boundary_right = [boundary_right]

    assert len(boundary_left) == len(boundary_right), "Lengths of boundary lists do not match!"
    assert (boundary_left <= boundary_right).all(0), "Left boundary must be smaller than right boundary"

    within_all = []
    for i in range(len(x)):
        within = 0
        j=0
        while not within and j < len(boundary_left):
            if boundary_left[j] < x[i] < boundary_right[j]: within = 1
            j = j + 1

        within_all.append(within)

    return within_all


def continous_idx(idx, sort=True):
    """
    Split input index array into multiple continuous subsets.

    :param idx: input index to be split
    :type idx: int array
    :param sort: if input index should be sorted
    :type sort: bool

    :return: idx_cont, continuous subset of input index
    :rtype: list
    """

    if sort: idx.sort()

    idx_last = idx[0]
    idx_cont_tmp = [idx[0]]
    idx_cont = []

    for i in range(len(idx)-1):
        if idx[i+1] == (idx_last + 1):
            idx_cont_tmp.append(idx[i+1])
        else:
            idx_cont.append(idx_cont_tmp)
            idx_cont_tmp = [idx[i+1]]

        idx_last = idx[i+1]

    idx_cont.append(idx_cont_tmp)

    return idx_cont


if __name__ == '__main__':
    for i in range(len(mode_all)):
        print(mode_all[i])