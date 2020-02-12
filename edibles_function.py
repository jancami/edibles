import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
from scipy.special import wofz
import astropy.constants as cst



##### Continuum Fitting Related #####
mode_all = ['spline', 'polynomial']

def polynomial_continuum(data_tuple, n_order=3, errweight = True):
    '''
        This function fits a polynomial continuum of given order to the input data
        given the nature of the data, x and y errs were not considered, but the
        sigma (of y) assuming Poisson Error can be used as weights
        draft by H.F. on Jan. 24th, 2020

        INPUT:
        data_tuple:   [tuple: (ndarray, ndarray)] (wavelength grid, flux values)
        errweight:    [bool]                      flag to specific use of Poisson Err as Weight
        n_order:      [int, default=3]            max order of polynomial

        OUTPUT:       [Angstroms or km/s]
        poly:       [PolyCont]                  fitted continuum level
    '''
    (x, y) = data_tuple
    assert len(x) == len(y), 'x and y must have the same length'
    weights = np.ones_like(x)
    if errweight:
        weights = weights / np.sqrt(y)
    coe_fit = np.polyfit(x, y, n_order, w=weights)
    poly = np.poly1d(coe_fit)

    return poly

def spline_continuum(data_tuple, pars=None, n_piece=2):
    '''
        This function fits a continuum to data separated into n sections
        where the x and y-values are the median of each section using a cubic spline
        This function is based on the generate_continuum from continuum_guess

        INPUT:        [Angstroms or km/s]
        data_tuple:   [tuple: (ndarray, ndarray)] (wavelength grid, flux values)
        pars:         [list]                      input y_points to fit spline
        n_piece:      [int, default=2]            number of sections to split data into

        OUTPUT:       [Angstroms or km/s]
        spline:       [CubicSpline]               fitted continuum level
        y_points:     [ndarray]                   y_points used in spline fitting

        Graphic:

        X------|-:    |    :-|-----X
        |      |  :   |   :  |     |
        |      |   :  X  :   |     |
        |      |    :-|-:    |     |
        |______|______|______|_____|

        <----piece---><---piece---->
        <-sec-><-sec-><-sec-><-sec->

        '''
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
    y_spline = spline(x)

    # plt.plot(x_points, y_points, 'kx', markersize='8', label='Points')

    return spline, y_points

def iterate_continuum(data_tuple,
                       mode='spline', n=3,
                       lsigma=1, usigma=2,
                       iterates=30, min_sigma = 0.2,
                       mask = None):
    '''
    This function fit the continuum level in a iterate manner
    Will call for other continuum function in the module

    INPUT:
    data_tuple:   [tuple: (ndarray, ndarray)] (wavelength grid, flux values)
    mode:         [string, default = 'spline'] fitting mode, only spline is available, but more will be added
    n:            [int, default  = 2]          multipurpose, no. of pieces in spline, max order in polynomial
    usigma:       [int, default = 2]           allowed n-sigma outlier in emission
    lsigma:       [int, default = 1]           allowed n-sigma outliers in absorption
    iterates:     [int, default = 30]          allowed max iterate times
    min_sigma:    [float, default = 0.1]       threshold of allowed usigma and lsigma, used to seize the iteration
    mask          [n-d array, default = None]  marker for masked points

    OUTPUT:       [Angstroms]
    y_count:      [CubicSpline/poly1D]         the fitted continuum
    idx:          [ndarray]                    index of continuum region
    '''

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

##### Voigt Calculation Related #####
def voigtMath(x, alpha, gamma):
    """
    Function to return the Voigt line shape centered at cent with Lorentzian component HWHM gamma
    and Gaussian component HWHM alpha.

    Creates a Voigt line profile using the scipy.special.wofz, which returns
    the value of the Faddeeva function.

    WARNING
    scipy.special.wofz is not compaible with np.float128 type parameters.

    Parameters
    ----------
    x : float64
        Dimensionless point/array
    alpha : float64
        Gaussian HWHM component
    gamma : float64
        Lorentzian HWHM component

    Returns
    -------
    ndarray
        Flux array for given input

    """

    sigma = alpha / np.sqrt(2 * np.log(2))

    return np.real(wofz((x + 1j*gamma)/sigma/np.sqrt(2))) / sigma/np.sqrt(2*np.pi)


def voigtOpticalDepth(lam, lam_0, b, d, Nf=1.0):
    """
    Converts parameters to make proper call to voigtMath


    Parameters
    ----------
    lam : float64
        Wavelength grid
    lam_0 : float64
        Central wavelength
    b : float64
        Gaussian standard deviation
    d : float64
        Damping parameter
    Nf : float64
        Scaling parameter, default = 1.0

    Returns
    -------
    ndarray
        Optical depth for given input

    """

    # convert lam & lam_0 to x
    x = lam - lam_0

    # convert b to sigma, then alpha
    sigma = b * lam_0 / cst.c.to('km/s').value
    alpha = sigma * np.sqrt(2. * np.log(2.))

    # convert d to gamma -- [ depends on what units we want to use ]

    # Currently, we are using the Lorentzian HWHM. This can easily be changed...
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gamma = d
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # create y data from voigtMath
    y = voigtMath(x, alpha, gamma)

    # Calculate tau_0
    tau_0 = np.pi * (cst.e.esu.value)**2 * Nf*(1e-8*lam_0)**2 / (cst.m_e.to('g').value*(cst.c.to('cm/s').value)**2)  # cm
    # Convert cm to angstroms
    tau_0 *= 1e8

    # Calculate tau
    tau = tau_0 * y

    # return scaled & shifted data
    return tau


def voigtAbsorptionLine(lam, lam_0, b, d, tau_0=0.1, N=None, f=None):
    """
    Function that takes in physical parameters and returns an absorption line.

    Choose either a tau_0 parameter, or N and f together. Default is tau_0.

    Parameters
    ----------
    lam : float64
        Wavelength grid
    lam_0 : float64
        Central wavelength
    b : float64
        Gaussian standard deviation
    d : float64
        Damping parameter
    N : float64
        Column density
    f : float64
        Oscillator strength
    tau_0 : float64
        Optical depth at center of line

    Returns
    -------
    ndarray
        flux array of light transmission

    """

    # create Nf
    if (N is not None) and (f is not None):
        Nf = N * f
    else:
        Nf = (tau_0 * 1e-8) * cst.m_e.to('g').value * (cst.c.to('cm/s').value)**2 / (np.pi * (cst.e.esu.value)**2 * (1e-8*lam_0)**2)

    tau = voigtOpticalDepth(lam, lam_0, b, d, Nf)

    transmission = np.exp(-tau)

    return transmission


##### User Interaction and Others #####

def go_message(message):
    '''
        A simple function to show a message and collect his/her decision

        INPUT
        message    [string]    message for the usser

        Output
        decision   [boll]      True or False
    '''
    go = ''
    while go not in ['Y','y','N','n']:
        go = input(message)

    if go in ['Y','y']:
        return True
    else:
        return False

def not_within_boundaries (x, boundary_left, boundary_right):
    """
        Test if input x is between the left and right boundaries

        Input:
        x                  float          Input to be examined
        boundary_left      [n-d Array]    List of left boundaries
        boundary_right     [n-d Array]    List of right boundaries; should have same length as boundary_left

        Output:
        within_result      [n-d Array]    Same length as x, each element True or False
    """

    assert len(boundary_left) == len(boundary_right), "Lengths of boundary lists do not match!"
    assert (boundary_left <= boundary_right).all(0), "Left boundary must be smaller than right boundary"

    try:len(boundary_left)
    except: boundary_left = [boundary_left]
    try:len(boundary_right)
    except: boundary_right = [boundary_right]

    not_within = 1
    j = 0
    while not_within and j < len(boundary_left):
        if boundary_left[j] <= x <= boundary_right[j]: not_within = 0
        j = j+1

    return not_within

def continous_idx(idx, sort=True):
    """
    This function split the input index array into multiple continous subsets

    INPUT:
    idx       [n-d int array]        containing the index to split
    sort      [bool, default=True]   if the input index will be sorted

    Output:
    idx_cont  [list of lists]        each element is a set of continous subset of input index
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