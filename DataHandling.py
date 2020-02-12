import numpy as np
from scipy.interpolate import CubicSpline
import astropy.constants as cst


##### Continuum Fitting Related #####
mode_all = ['spline', 'polynomial']

def polynomial_continuum(data_tuple, n_order=3, errweight = True):
    '''
        This function fits a n-th order polynomial continuum to input data
        x and y errs were not considered, but the
        sigma (of y) assuming Poisson Error can be used as weights
        draft by H.F. on Jan. 24th, 2020

        INPUT:
        data_tuple:   [tuple: (ndarray, ndarray)] (wavelength grid, flux values)
        errweight:    [bool]                      flag to specific use of Poisson Err as Weight
        n_order:      [int, default=3]            max order of polynomial

        OUTPUT:
        poly:         [PolyCont]                  fitted continuum level
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

##### Wavelength shifted/ convert between A and km/s #####
# to be added soon