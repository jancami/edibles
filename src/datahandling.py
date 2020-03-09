import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.signal import find_peaks, peak_prominences
import astropy.constants as cst
import edibles.src.math as eMath



##### Continuum Fitting Related #####
mode_all = ["spline", "s", "polynomial", "p"]


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
    if mode[0].lower() == "s":npiece = n
    if mode[0].lower() == "p":norder = n

    while do_flag:
        idx = [i for i in range(len(x)) if -usigma <= y_nsigma[i] <= lsigma and mask[i] == 0]
        x_fit = np.copy(x[idx])
        y_fit = np.copy(y[idx])

        if mode[0].lower() == "s":
            (y_count, y_points) = spline_continuum((x_fit,y_fit), pars=None, n_piece=npiece)

        if mode[0].lower() == "p":
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


def estimateSNR(wave, flux, normalized=True, n=None, silenced=False):
    """
    Estimate the SNR at continuum region of input data from continuum region defined by a a constant continuum
    It's recommended that input data have been normalized. If not, will fit a spline continuum first
    A plot showing the continuum region will be given, but can be silenced
    The approach is different from what was in the repo, which was based on spectrum slicing

    :param wave: input wavelength grid
    :type wave: nparray
    :param flux: input flux, normalized or not normalized, same length as wave
    :type flux: nparray
    :param normalized: if the input flux has been normalized
    :type normalized: bool
    :param silenced: if set true, no figure will be given for visual examine
    :type silenced: bool

    :return: SNR, idx, estimated SNR and index of points identified as continuum region
    :rtype: float, int list
    """

    assert len(wave) == len(flux), "input data must have same length"
    n_points = len(wave)
    # normalized input data first
    if not normalized:
        if n is None:
            print("input data not normalized")
            n=float(input("type the number of sections for your spline continuum"))
        (cont, idx_spline) = iterate_continuum((wave,flux), mode="spline", n=n)

        flux = flux/cont(wave)
        idx_masked = [i for i in range(n_points) if i not in idx_spline]
    else:
        idx_masked = []

    std_res = np.std(np.abs(np.ones(n_points) - flux))
    idx_masked = idx_masked + [i for i in range(n_points) if np.abs(flux[i] - 1.) > 0.5 * std_res]
    idx_masked = np.unique(idx_masked)
    mask = np.zeros(n_points)
    mask[idx_masked] = 1
    (cont, idx_cont) = iterate_continuum((wave,flux), mode='polynomial', n=0, min_sigma=0.1, mask=mask)
    SNR = (1 / np.std(flux[idx_cont]))

    if not silenced:
        idx_cont_split = continous_idx(idx_cont)

        fig, ax = plt.subplots(1,1)
        ax.plot(wave,flux,color="k")
        for i in range(len(idx_cont_split)):
            ax.plot(wave[idx_cont_split[i]], flux[idx_cont_split[i]], marker="x", markersize=3, color="red")
        ax.set_ylabel("Normalized Flux")
        ax.set_xlabel("AA")
        ax.grid()
        plt.show()

    return SNR, idx_cont


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
    wave_out = wave_array + (v_shift / cst.c.to("km/s").value) * wave_array
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




def vac2air(vacw, mode="Ciddor"):
    """
    Convert Vacuum wavelengths to air wavelengths, both in AA, for wavleength > 2000AA
    Two modes are provided:
    Ciddor, based on Ciddor 1996, Applied Optics LP, Vol.35, Issue 9, pp 1556, claimed to be more accurate in IR
    Morton, based on Morton 1991, ApJS, 77, 119

    :param vacw: vacuum wavelengths to be converted
    :type param: nparray
    :param mode: mode to use, Ciddor or Morton, Ciddor by default
    :type mode: str

    :return: airw, converted air wavelength
    :rtype: nparray
    """
    mode_all = ["ciddor", "morton"]
    assert mode.lower() in mode_all, "Available Modes: {mode_all}".format(mode_all = mode_all)

    if mode == "Ciddor":
        airw = eMath.vac2air_ciddor(vacw)

    if mode == "Morton":
        airw = eMath.vac2air_morton(vacw)

    return airw

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

    :return: within_all, idx_all: within_all is a 1/0 array where 1 means the elements is in the boundaries;
                                  idx are the index of all points within the boundaries
    :rtype: list
    """
    try:len(x)
    except:x = [x]
    try:len(boundary_left)
    except: boundary_left = [boundary_left]
    try:len(boundary_right)
    except: boundary_right = [boundary_right]

    assert len(boundary_left) == len(boundary_right), "Lengths of boundary lists do not match!"
    if len(boundary_left) > 1:
        assert (boundary_left <= boundary_right).all(0), "Left boundary must be smaller than right boundary"
    else:
        assert boundary_left <= boundary_right, "Left boundary must be smaller than right boundary"

    within_all = []
    for i in range(len(x)):
        within = 0
        j=0
        while not within and j < len(boundary_left):
            if boundary_left[j] < x[i] < boundary_right[j]: within = 1
            j = j + 1

        within_all.append(within)

    idx_all = [i for i in range(len(within_all)) if within_all[i] == 1]

    return np.array(within_all), idx_all


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


def nearest_point(point_array, spec_array, scale=True):
    """
    Find a point on the spectrum that is closest to each of the input points

    :param point_array: coordinates of input points, a tuple of two lists of the same length
    :type point_array: tuple
    :param spec_array: coordinates of input points, a tuple of two lists of wave and flux
    :type spec_array: tuple
    :param scale: if set true, the coordinates will be scaled to span over 0 and 1, default to true
    :type scal: bool

    :return: point_idx, the index of points in the spec that are most close to the input points
    :rtype: int array
    """
    assert len(point_array) == 2, "Inputs must have 2 dimension"
    assert len(spec_array) == 2, "Inputs must have 2 dimension"

    n_points = len(point_array[0])

    x_min = np.min(spec_array[0])
    x_span = np.max(spec_array[0]) - x_min
    y_min = np.min(spec_array[1])
    y_span = np.max(spec_array[1]) - y_min
    if scale:
        point_array_scal = ((point_array[0] - x_min) / x_span,)
        point_array_scal = point_array_scal + ((point_array[1] - y_min) / y_span,)

        spec_array_scal = ((spec_array[0] - x_min) / x_span,)
        spec_array_scal = spec_array_scal + ((spec_array[1] - y_min) / y_span,)

    else:
        point_array_scal = point_array
        spec_array_scal = spec_array

    d_matrix_all = eMath.all_distance(point_array_scal, spec_array_scal)
    point_idx = []

    for n in range(n_points):
        d_matrix_point = d_matrix_all[n]
        point_idx.append(np.where(d_matrix_point == np.min(d_matrix_point))[0][0])

    return point_idx

def searchPeak(flux, n_peaks=None, normalized=True, SNR=None, prominence=3):
    """
    Searh for peaks in the input flux. Each peak should be significant enough, as defined by SNR and prominence.
    Flux will be normalized if not. SNR will be estimated if not given.

    :param flux: input flux array
    :type flux: nparray
    :param n_peaks: number of peaks to be raised, if set, will return the index of the most significant n peaks
    :type n_peaks: int
    :param normalized: if input flux has been normalized, if not, will be normalized first
    :type normalized: bool
    :param SNR: SNR of input flux, used to define the sigma if prominence
    :type SNR: float
    :param prominence: required n-sigma significance, default to 3
    :type prominence: float

    :return: peak_idx, the index of peaks int he flux array
    :rtype: index array
    """
    if not normalized:
        x = np.arange(len(flux))
        (cont, idx_spline) = iterate_continuum((x, flux), mode="spline", n=7)
        flux = flux / cont(x)

    if SNR is None:
        x = np.arange(len(flux))
        (SNR, idx_cont) = estimateSNR(x, flux, normalized=False, n=7, silenced=True)

    # find the peaks
    flux = -flux # in emission
    peak_idx, _ = find_peaks(flux, prominence = prominence/SNR)
    prominences =peak_prominences(flux,peak_idx)[0]

    # if number of peaks specificed, find the most significant peaks
    if n_peaks is not None and n_peaks < len(peak_idx):
        threshold = np.sort(prominences)[-n_peaks]
        idx = [i for i in range(len(peak_idx)) if prominences[i] >= threshold]
        peak_idx = peak_idx[idx]

    return peak_idx

def parseInput(n,*inputs, checklen=True):
    """
    Sort the inputs and duplicate them into list of length n (iterable).
    If already iterable (list or np.array), check if they have the length n.
    This function is to be used in methods of ediblesSpectrum that might handle multiple spectrum panels

    :param n: required length of the inputs
    :type n: in
    :param inputs: all inputs that need to be duplicated and examined
    :type inputs: tuple
    :param checklen: if set, check if length of all list type inputs equal to n, set to False when convert "panels"
    :type checklen: bool
    :return: inputs_out: duplicated inputs
    :rtype: tuple
    """

    inputs_out = ()
    iterable_types = [type([1,2]), type(np.array([1,2])), (1,2)]

    for input in inputs:
        if type(input) in iterable_types:
            if checklen:
                assert len(input) == n, "Length of each input must be 1 or the same as panel numbers"
            inputs_out = inputs_out + (np.array(input),)
        else:
            if type(input) is list:
                if len(input) == 1:
                    inputs_out = inputs_out + (input * n,)
                else:
                    assert len(input) == n, "Length of each input must be 1 or the same as panel numbers"
            else:
                inputs_out = inputs_out + ([input] * n,)

    if len(inputs_out) == 1: inputs_out = inputs_out[0]
    return inputs_out


if __name__ == '__main__':
    for i in range(len(mode_all)):
        print(mode_all[i])