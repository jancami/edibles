import numpy as np
from scipy.special import wofz
import astropy.constants as cst

##### Voigt Calculation #####
def voigtMath(x, alpha, gamma):
    """
    Returns a NORMALIZED Voigt profile whose peak area is unity.
    The Voigt profile peaks at zero with Lorentzian component HWHM gamma
    and Gaussian component HWHM alpha.
    Voigt profile created using scipy.special.wofz, which returns the value
    of a Faddeeva function.
    WARNING
    scipy.special.wofz is not compaible with np.float128 type parameters.
    Units for x, alpha, and gamma not specified but should be the same

    :param x: wavelength grid centered at zero
    :type x: ndarray of float64
    :param alpha: Gaussian HWHM component, in AA
    :type alpha: float
    :param gamma: Lorentzian HWHM component, in AA
    :type gamma: float

    :return: standard_voigt, y points of a Voigt profile as defined by inputs
    :rtype: ndarray of float
    """
    sigma = alpha / np.sqrt(2 * np.log(2))
    standard_voigt = np.real(wofz((x + 1j*gamma)/sigma/np.sqrt(2))) / sigma/np.sqrt(2*np.pi)

    return standard_voigt


def voigtAbsorption(lam, lam_0, b, d, area = 1.0):
    """
    Returns a Voigt absorption defined by lam_0, b, d, and peak area of "area"
    Can be in flux (e.g. for DIBs) or optical depth unit (e.g. for known ISM line).
    Need to be "exponentialed" if in optical depth unit.
    :param lam: Wavelength grid, in AA
    :type lam: ndarray of float
    :param lam_0: central wavelength, in AA
    :type lam_0: float
    :param b: Gaussian standard deviation, in km/s
    :type b: float
    :param d: Lorentzian damping parameter, in AA
    :type d: float
    :param area: peak area of the profile, in AA, default to 1.0
    :type area: float

    :return: voigt_absorption, y points of voigt absorption as defined by lam_0, b, d, and area
    :rtype: ndarray
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

    voigt_absorption = - y * area
    return voigt_absorption


    # Calculate tau_0
    #tau_0 = np.pi * (cst.e.esu.value) ** 2 * Nf * (1e-8 * lam_0) ** 2 / (
     #           cst.m_e.to('g').value * (cst.c.to('cm/s').value) ** 2)  # cm
    # Convert cm to angstroms
    #tau_0 *= 1e8


def voigtOpticalDepthAbsorption(lam, lam_0, b, d, tau_0=0.1, N=None, f=None):
    """
    Returens an absorption whose OPTICAL DEPTH has Voigt profile, as defined by lam_0, b, d.
    Total Optical depth defined by tau_0 or N the column density and f the oscillator strength.
    N and f, if set, will override tau_0.
    The calculation is tau-based, which is more straight-forward

    :param lam: wavelength grid, in AA
    :type lam: ndarray of float
    :param lam_0: central wavelength, in AA
    :type lam_0: float
    :param b: Gaussian standard deviation, in km/s
    :type b: float
    :param d: Lorentzian damping parameter, in AA
    :type d: float
    :param tau_0: TOTAL optical depth of the absorption, overrode by N and f, default to 1.0
    :type tau_0: float
    :param N: column density, in cm^-2
    :type N: float
    :param f: oscillator strength
    :type f: float

    :return: transmission, ypoints of absorption whose optical depth is Vogit profile defined by inputs
    :rtype: ndarray
    """

    # create tau_0
    if (N is not None) and (f is not None):
        Nf = N * f
        tau_0 = np.pi * (cst.e.esu.value) ** 2 * Nf * (1e-8 * lam_0) ** 2 / (
                    cst.m_e.to('g').value * (cst.c.to('cm/s').value) ** 2)  # cm
        # Convert cm to angstroms
        tau_0 *= 1e8

    voigt_opticaldepth = voigtAbsorption(lam, lam_0, b, d, area = tau_0)

    transmission = np.exp(voigt_opticaldepth)

    return transmission


def euclidean_distance(point1,point2):
    """
    Calculate the Euclidean distance between point A and point B

    :param point1: coordinate array of point A
    :type point1: tuple
    :param point2: coordinate array of point B
    :type point2: tuple of same length as point 1

    :return: distance, the Euclidean distance
    :rtype: float
    """
    assert len(point1) == len(point2), "Inputs points must have the same dimensions"
    dist_ndim = 0
    for i in range(len(point1)):
        dist_ndim = dist_ndim + (point1[i] - point2[i])**2

    dist = np.sqrt(dist_ndim)

    return dist


def all_distance(point_array1, point_array2):
    """
    Calculate the distance matrix of point arraies 1 and 2

    :param point_array1: coordinates for point array A, a tuple containing d (dimension) lists with length of N (points)
    :type point_array1: tuple
    :param point_array2: coordinates for point array A, a tuple containing d (dimension) lists with length of M (points)
    :type point_array2: tuple of same length as point_array1

    :return: d_matrix, the N*M distance matrix
    :rtype: 2d-array
    """
    assert len(point_array1) == len(point_array2), "Inputs points must have the same dimensions"

    d_total = len(point_array1)
    N, M = len(point_array1[0]), len(point_array2[0])

    d_matrix = np.zeros(shape=(N, M))
    for i in range(N):
        point1 = ()
        for d in range(d_total):
            point1 = point1 + (point_array1[d][i],)

        for j in range(M):
            point2 = ()
            for d in range(d_total):
                point2 = point2 + (point_array2[d][j],)

            d_matrix[i][j] = euclidean_distance(point1, point2)

    return d_matrix


def vac2air_morton(vacw):
    """
    Convert vaccum wavelengths to air wavelengths, both in AA
    Uses relations from Morton 1991, ApJS, 77, 119, valid for wavelengths > 2000 AA

    :param vacw: vaccum wavelegnths
    :type vacw: nparray
    :return: airw, converted air wavelength
    :rtype: nparray
    """

    temp = (1e4 / vacw) ** 2
    airw = 1. / (1. + 6.4328e-5 + 2.94981e-2 / (146 - temp) +
                 2.5540e-4 / (41 - temp)) * vacw
    return airw


def vac2air_ciddor(vacw):
    """
    Convert vaccum wavelengths to air wavelengths, both in AA
    users relations from Ciddor 1996, Applied Optics LP, vol. 35, Issue 9, pp1566
    Only for wavelengths > 2000 AA

    :param vacw: vaccum wavelengths
    :type vacw: nparray

    :return: airw
    :rtype: nparray
    """

    k0 = 238.0185
    k1 = 1e-8 * 5792105.
    k2 = 57.362
    k3 = 1e-8 * 167917.
    s2 = (1e4 / vacw) ** 2
    n = 1 + k1 / (k0 - s2) + k3 / (k2 - s2)
    airw = vacw / n

    return airw



if __name__ == '__main__':
    print("Math modular for EDIBLES spectrum")
