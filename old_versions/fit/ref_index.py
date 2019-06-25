#! -*- coding: utf-8 -*-
"""Refractive index of air.

NIST provides an online calculator for calculating refractive index of
air, for light of a certain wave length, under varying atmospheric
conditions. This module implements the equations provided in the
documentation for the online calculator.

In addition to calculating the refractive index, this module also has
functions for converting wave length of light in vacuum to that in air,
and vice-versa.

The documentation for the online calculator is provided at
http://emtoolbox.nist.gov/Wavelength/Documentation.asp, and includes a
link to the online calculator.

The following comments are based on the discussions presented in the
NIST documentation. It is intended as a brief overview. See
http://emtoolbox.nist.gov/Wavelength/Documentation.asp, for detailed
discussions.

Refractive index of air can be caclulated using two different
algorithms: one due to Edlén (updated by Birch and Down), and one due
to Ciddor. The latter has been adopted by the International Association
of Geodesy (IAG) as the reference equation for calculating refractive
index of air. Functions for calculating refractive index using either
of these are defined in this module.

The vacuum to air and air to vacuum wave length conversion functions in
this module use the Ciddor equation, in the form presented in the NIST
documentation.

Uncertainities in refractive index, and hence in wave length
conversions, due to uncertanities in measured values of temperature,
pressure, and humidity exceeds that due to the intrinsic uncertainity
in the equations used.

An uncertainty of 1e-6 in refractive index can result from a
combination of:

  + an error of 1°C (1.8 °F) in air temperature

  + an error of 0.4kPa (3mm of Hg) in air pressure

  + an error of 50% in relative humidity at sufficiently high air
    temperatures (near 35°C)

Valid range for input parameters for the refractive index calculations
are presented below. The online calculator issues a warning if input
parameters are outside a smaller interval within the maximum
range. Functions in this module do not raise a warning by default. But
they accept a keyword ``warn``, which when set to ``True`` will result
in warnings, when the input parameters are outside the accepted range.

  + Wavelength [300nm - 1700nm]

    Warning is issued if value is outside [350nm - 1600nm].

  + Pressure [10kPa - 140kPa]

    Warning is issued if value is outside [60kPa - 120kPa].

  + Temperature [-40∘C - 100∘C].

    Warning is issued if value is outside [0∘C - 40∘C].

  + Humidity [0 - 100]

    Can be given as relative humidity, dew point, frost point or
    partial pressure of water vapour. A warning is given if the mole
    fraction of water vapour exceeds 20% or, equivalently, relative
    humidity exceeds 85%. A warning is issued if relative humidity is
    less than 1%.

  + CO2 concentration [0µmole/mole - 2000µmole/mole]

    The common value to use is 450. Outdoor values are rarely below 300
    and indoor can be as high as 600. A difference of 150 will lead to
    a difference of only ~ 2e-8 in index of refraction.

    A warning is issued if a value other than 450 is used.


In astronomy, the convention is to use the refraction correction for
wave length greater than 200nm, eventhough the equations are not
strictly valid at wave lengths shorter than 300nm. For example, the
popular IDLASTRO IDL code vactoair.pro and airtovac.pro will accept any
wave length greater than 2000Å.

To accomodate this type of usage, instead of limiting the possible
input wave lengths, functions in this module will accept any wave
length value. It is up to the user to decide if a particular wave
length is to be used as an input to the equations.

Comparison with the IDLASTRO vactoair.pro and airtovac.pro algorithms
show that the equivalent functions in this module, vac2air and air2vac,
give results that agree to within 1e-4nm, over a range of wavelengths
from 200nm to 1700nm. This uncertainty translates to a velocity
difference of 150m/s to 17m/s, over the wave length range 1700nm to
200nm.

The IDLASTRO code uses a fixed value of temperature and humidity which
is not documented in the code. The above comparison was carried out at
a temperature of 15∘C and a relative humidity of 0.

The IDL code used for testing was downloaded on 2011/10/07. The
revision history indicates that the IDL code in vactoair.pro and
airtovac.pro were last modified in March 2011.

:author: Prasanth Nair
:contact: prasanthhn@gmail.com
:license: BSD (http://www.opensource.org/licenses/bsd-license.php)
"""
from __future__ import division
from __future__ import print_function
import logging

__version__ = "1.0"

logger = logging.getLogger("ref_index")
logger.setLevel(logging.DEBUG)
h = logging.StreamHandler()
h.setLevel(logging.DEBUG)
f = logging.Formatter("%(name)s: %(levelname)s - %(message)s")
h.setFormatter(f)
logger.addHandler(h)


def f2k(f):
    """Converts Fahrenheit to Kelvin."""
    return (f - 32.0) * (100.0 / 180.0) + 273.15


def k2f(k):
    """Converts Kelvin to Fahrenheit."""
    return (k - 273.15) * (180.0 / 100.0) + 32.0


def c2k(c):
    """Converts Celsius to Kelvin."""
    return c + 273.15


def k2c(k):
    """Converts Kelvin to Celsius."""
    return k - 273.15


def c2f(c):
    """Converts Celsius to Fahrenheit."""
    return c * (180.0 / 100.0) - 32.0


def f2c(f):
    """Converts Fahrenheit to Celsius."""
    return (f - 32.0) * (100.0 / 180.0)


def svp_water(t):
    """Saturation vapour pressure over water at given temperature.

    Parameters
    ----------
    t : float
        Air temperature in degree Celsius.

    Returns
    -------
    p_sv : float
        Saturation vapour pressure over water, at the given
        temperature, in Pascal.

    Notes
    -----
    From section A-I of
    http://emtoolbox.nist.gov/Wavelength/Documentation.asp.

    """
    import math
    K1 = 1.16705214528e+03
    K2 = -7.24213167032e+05
    K3 = -1.70738469401e+01
    K4 = 1.20208247025e+04
    K5 = -3.23255503223e+06
    K6 = 1.49151086135e+01
    K7 = -4.82326573616e+03
    K8 = 4.05113405421e+05
    K9 = -2.38555575678e-01
    K10 = 6.50175348448e+02

    T = t + 273.15
    omega = T + K9 / (T - K10)
    A = omega ** 2 + K1 * omega + K2
    B = K3 * omega ** 2 + K4 * omega + K5
    C = K6 * omega ** 2 + K7 * omega + K8
    X = -B + math.sqrt(B ** 2 - 4 * A * C)

    p_sv = 1.0e6 * ((2.0 * C / X) ** 4)

    return p_sv


def svp_ice(t):
    """Saturation vapour pressure over ice at given temperature.


    Parameters
    ----------
    t : float
        Temperature in degree Celsius.

    Returns
    -------
    p_sv : float
        Saturation vapour pressure over ice, at the given
        temperature, in Pascal.

    Notes
    -----
    From section A-I of
    http://emtoolbox.nist.gov/Wavelength/Documentation.asp.

    """
    import math
    A1 = -13.928169
    A2 = 34.7078238

    t += 273.15
    theta = t / 273.16
    Y = A1 * (1 - theta ** -1.5) + A2 * (1 - theta ** -1.25)

    p_sv = 611.657 * math.exp(Y)

    return p_sv


def dew_point_wvpp(td):
    """Water vapour saturation pressure, given dew point temperature."""
    return svp_water(td)


def frost_point_wvpp(tf):
    """Water vapour saturation pressure, given frost point temperature."""
    return svp_ice(tf)


def rh2wvpp(rh, t):
    """Convert relative humidity to water vapour partial pressure.

    Parameters
    ----------
    rh : float
        Relative humidity as a number between 0 and 100.
    t : float
        Temperature in degree Celsius.

    Returns
    -------
    p_sv : float
        Water vapour partial pressure, in Pascal.

    Notes
    -----
    See section A-II of
    http://emtoolbox.nist.gov/Wavelength/Documentation.asp.

    """
    # t > 0 according to documentation.
    if t >= 0:
        p_sv = svp_water(t)
    elif t < 0:
        p_sv = svp_ice(t)

    return (rh / 100.0) * p_sv


def f_factor(p, t):
    """Enhancement factor for calculating mole fraction.

    Parameters
    ----------
    p : float
        Pressure in Pascal.
    t : float
        Temperature in degree Celsius.

    Returns
    -------
    f : float
        Enhancement factor needed in calculation of mole fraction.

    Notes
    -----
    See section A-II of
    http://emtoolbox.nist.gov/Wavelength/Documentation.asp.

    """
    alpha = 1.00062
    beta = 3.14e-8
    gamma = 5.60e-7

    return alpha + beta * p + gamma * (t ** 2)


def dew_point_mole_fraction(p, t):
    """Water vapour mole fraction for given dew point temperature.
    Parameters
    ----------
    p : float
        Pressure in Pascal.
    t : float
        Temperature in degree Celsius.

    Returns
    -------
    xv : float
        Mole fraction.

    Notes
    -----
    See section A-II of
    http://emtoolbox.nist.gov/Wavelength/Documentation.asp.

    """
    return f_factor(p, t) * dew_point_wvpp(t) / p


def frost_point_mole_fraction(p, t):
    """Water vapour mole fraction for given frost point temperature.
    Parameters
    ----------
    p : float
        Pressure in Pascal.
    t : float
        Temperature in degree Celsius.

    Returns
    -------
    xv : float
        Mole fraction.

    Notes
    -----
    See section A-II of
    http://emtoolbox.nist.gov/Wavelength/Documentation.asp.

    """
    return f_factor(p, t) * frost_point_wvpp(t) / p


def rh2mole_fraction(rh, p, t):
    """Water vapour mole fraction from relative humidity.

    Parameters
    ----------
    rh : float
        Relative humidity as a number between 0 and 100.
    p : float
        Pressure in Pascal.
    t : float
        Temperature in Kelvin.

    Returns
    -------
    xv : float
        Mole fraction.

    Notes
    -----
    See section A-II of
    http://emtoolbox.nist.gov/Wavelength/Documentation.asp.

    """
    return f_factor(p, t) * rh2wvpp(rh, t) / p


def pp2mole_fraction(pv, p, t):
    """Water vapour mole fraction from partial pressure.

    Parameters
    ----------
    rh : float
        Relative humidity as a number between 0 and 100.
    p : float
        Pressure in Pascal.
    t : float
        Temperature in Kelvin.

    Returns
    -------
    xv : float
        Mole fraction.

    Notes
    -----
    See section A-II of
    http://emtoolbox.nist.gov/Wavelength/Documentation.asp.

    """
    return f_factor(p, t) * pv / p


def _check_range(**kwargs):
    """Return True if value is inside accepted range."""
    if not (350 <= kwargs.get('wave', 633) <= 1600):
        logger.warning("Wave length outside [350nm, 1600nm].")
    if not (60000 <= kwargs.get('p', 101325) <= 120000):
        logger.warning("Pressure outside [60000Pa - 120000Pa].")
    if not (0 <= kwargs.get('t', 20) <= 40):
        logger.warning("Temperature outside [0C - 40C].")
    if not (1 < kwargs.get('rh', 50) <= 85):
        logger.warning("Relative humidity outside (1 - 85].")
    if not (kwargs.get('xv', 0.4) >= 0.2):
        logger.warning("Mole fraction less than 0.2.")
    if kwargs.get('co2', 450) != 450:
        logger.warning("CO2 concentration is not 450.")


def ciddor_ri(wave, t, p, xv, co2=450, warn=False):
    """Refractive index of air according to the Ciddor equation.

    Parameters
    ----------
    wave : float or Numpy array of float
        Wavelength in vacuum, in nano-meters. Valid wavelength range is
        300nm - 1700nm.
    t : float
        Temperature in degree Celsius. Valid temperate range is -40 to
        100 degree Celsius.
    p : float
        Pressure in Pascal. Valid range is from 10kPa - 140 kPa.
    xv : float
        Water vapour mole fraction, as a number between 0 and
        1. Default is set to 0.
    co2 : float
        Carbon dioxide concentration in µmole/mole. The default value
        of 450 should be enough for most purposes. Valid range is from
        0 - 2000 µmole/mole.
    warn : bool
        Warning is issued if parameters fall outside accept
        range. Accepted range is smaller than the valid ranges
        mentioned above. See module docstring for accepted ranges.

        The default is False and no warnings are issued.

    Notes
    -----
    See section A-III of
    http://emtoolbox.nist.gov/Wavelength/Documentation.asp.

    See
    """
    if warn:
        _check_range(wave, t, p, xv)

    w0 = 295.235
    w1 = 2.6422
    w2 = -0.03238
    w3 = 0.004028
    k0 = 238.0185
    k1 = 5792105
    k2 = 57.362
    k3 = 167917
    a0 = 1.58123e-6
    a1 = -2.9331e-8
    a2 = 1.1043e-10
    b0 = 5.707e-6
    b1 = -2.051e-8
    c0 = 1.9898e-4
    c1 = -2.376e-6
    d = 1.83e-11
    e = -0.765e-8
    pr1 = 101325
    tr1 = 288.15
    Za = 0.9995922115
    rhovs = 0.00985938
    R = 8.314472
    Mv = 0.018015

    wave = wave * 1.0e-3
    S = 1.0 / wave ** 2

    ras = 1e-8 * ((k1 / (k0 - S)) + (k3 / (k2 - S)))
    rvs = 1.022e-8 * (w0 + w1 * S + w2 * S ** 2 + w3 * S ** 3)

    Ma = 0.0289635 + 1.2011e-8 * (co2 - 400.0)

    raxs = ras * (1 + 5.34e-7 * (co2 - 450.0))

    T = t + 273.15

    Zm = a0 + a1 * t + a2 * t ** 2 + (b0 + b1 * t) * xv + \
        (c0 + c1 * t) * xv ** 2
    Zm *= -(p / T)
    Zm += (p / T ) ** 2 * (d + e * xv ** 2)
    Zm += 1

    rhoaxs = pr1 * Ma / (Za * R * tr1)

    rhov = xv * p * Mv / (Zm * R * T)

    rhoa = (1 - xv) * p * Ma / (Zm * R * T)

    n = 1.0 + (rhoa / rhoaxs) * raxs + (rhov / rhovs) * rvs

    return n


def ciddor(wave, t, p, rh, co2=450, warn=False):
    """Refractive index of air according to the Ciddor equation.

    Accepts relative humidity instead of mole fraction, as done in
    ``ciddor_ri()``.

    Parameters
    ----------
    wave : float or Numpy array of float
        Wavelength in vacuum, in nano-meters. Valid wavelength range is
        300nm - 1700nm.
    t : float
        Temperature in degree Celsius. Valid temperate range is -40 to
        100 degree Celsius.
    p : float
        Pressure in Pascal. Valid range is from 10kPa - 140 kPa.
    rh : float
        Relative humidity [0 - 100].
    co2 : float
        Carbon dioxide concentration in µmole/mole. The default value
        of 450 should be enough for most purposes. Valid range is from
        0 - 2000 µmole/mole.
    warn : bool
        Warning is issued if parameters fall outside accept
        range. Accepted range is smaller than the valid ranges
        mentioned above. See module docstring for accepted ranges.

        The default is False and no warnings are issued.

    Notes
    -----
    See section A-III of
    http://emtoolbox.nist.gov/Wavelength/Documentation.asp.

    """
    if warn:
        _check_range(wave, t, p, rh)
        # turn off warning, so that ciddor_ri doesn't issue duplicate
        # warning.
        warn = False

    xv = rh2mole_fraction(rh=rh, p=p, t=t)
    return ciddor_ri(wave=wave, t=t, p=p, xv=xv, co2=co2, warn=warn)


def edlen_ri(wave, t, p, pv, warn=False):
    """Refractive index of air according to the Edlén equation.

    Parameters
    ----------
    wave : float or Numpy array of float
        Wavelength in vacuum, in nano-meters. Valid wavelength range is
        300nm - 1700nm.
    t : float
        Temperature in degree Celsius. Valid temperate range is -40 to
        100 degree Celsius.
    p : float
        Pressure in Pascal. Valid range is from 10kPa - 140 kPa.
    pv : float
        Water vapour partial pressure, in Pascal.
    warn : bool
        Warning is issued if parameters fall outside accept
        range. Accepted range is smaller than the valid ranges
        mentioned above. See module docstring for accepted ranges.

        The default is False and no warnings are issued.

    Notes
    -----
    See section A-IV of
    http://emtoolbox.nist.gov/Wavelength/Documentation.asp.

    """
    if warn:
        _check_range(wave, t, p)

    A = 8342.54
    B = 2406147
    C = 15998
    D = 96095.43
    E = 0.601
    F = 0.00972
    G = 0.003661

    wave = wave * 1.0e-3
    S = 1.0 / wave ** 2

    ns = 1 + 1e-8 * (A + B / (130.0 - S) + C / (38.9 - S))

    X = (1 + 1e-8 * (E - F * t) * p) / (1 + G * t)

    ntp = 1 + p * (ns - 1) * X / D

    n = ntp - 1e-10 * ((292.75 / (t + 273.15)) * \
                           (3.7345 - 0.0401 * S)) * pv

    return n


def edlen(wave, t, p, rh, warn=False):
    """Refractive index of air according to the Edlén equation.

    Accepts relative humidity instead of water vapour partial pressure,
    as in ``edlen_ri()``.

    Parameters
    ----------
    wave : float or Numpy array of float
        Wavelength in vacuum, in nano-meters. Valid wavelength range is
        300nm - 1700nm.
    t : float
        Temperature in degree Celsius. Valid temperate range is -40 to
        100 degree Celsius.
    p : float
        Pressure in Pascal. Valid range is from 10kPa - 140 kPa.
    rh : float
        Relative humidity in [0 - 100].
    warn : bool
        Warning is issued if parameters fall outside accept
        range. Accepted range is smaller than the valid ranges
        mentioned above. See module docstring for accepted ranges.

        The default is False and no warnings are issued.

    Notes
    -----
    See section A-IV of
    http://emtoolbox.nist.gov/Wavelength/Documentation.asp.

    """
    if warn:
        _check_range(wave, t, p)
        # turn off warning so that edlen_ri() doesn't raise duplicate
        # warning.
        warn = False

    pv = rh2wvpp(rh=rh, t=t)
    return edlen_ri(wave=wave, t=t, p=p, pv=pv, warn=warn)


def vac2air(wave, t=15.0, p=101325, rh=0.0, co2=450, warn=False):
    """Wavelength of light in air, using Ciddor refractive index.

    Parameters
    ----------
    wave : float or Numpy array of float
        Wavelength in nano-meters. Valid range is 300nm - 1700nm.
    t : float
        Temperature in degree Celsius. Valid range is -40 - 100 degree
        Celsius. Default is 15 degree Celsius (288.15 Kelvin).
    p : float
        Pressure in Pascal. Valid range is 10kPa - 140kPa. Default is
        101325 Pa (1 atmosphere).
    rh : float
        Relative humidity as a number between 0 and 100. Default is 0.
    co2 : float
        Carbon dioxide concentration in µmole/mole. The default value
        of 450 is sufficient for most purposes. Valid range is 0 - 2000
        µmole/mole.
    warn : bool
        Warning is issued if parameters fall outside accept
        range. Accepted range is smaller than the valid ranges
        mentioned above. See module docstring for accepted ranges.

        The default is False and no warnings are issued.

    Returns
    -------
    w : float
        Wavelength in air, in nm.

    """
    if warn:
        _check_range(wave, t, p, rh, co2)

    n = ciddor(wave, t, p, rh, co2)
    return wave / n


def air2vac(wave, t=15.0, p=101325, rh=0.0, co2=450, warn=False):
    """Wavelength of light in vacuum, using Ciddor refractive index.

    The refractive index calculation needs wavelength in vacuum. In
    this function, the wavelength in air is used. The errors are on the
    order of 1e-5 nm.

    Parameters
    ----------
    wave : float or Numpy array of float
        Wavelength in nano-meters. Valid range is 300nm - 1700nm.
    t : float
        Temperature in degree Celsius. Valid range is -40 - 100 degree
        Celsius. Default is 15 degree Celsius (288.15 Kelvin).
    p : float
        Pressure in Pascal. Valid range is 10kPa - 140kPa. Default is
        101325 Pa (1 atmosphere).
    rh : float
        Relative humidity as a number between 0 and 100. Default is 0.
    co2 : float
        Carbon dioxide concentration in µmole/mole. The default value
        of 450 is sufficient for most purposes. Valid range is 0 - 2000
        µmole/mole.
    warn : bool
        Warning is issued if parameters fall outside accept
        range. Accepted range is smaller than the valid ranges
        mentioned above. See module docstring for accepted ranges.

        The default is False and no warnings are issued.

    Returns
    -------
    w : float
        Wavelength in vacuum, in nm.

    """
    if warn:
        _check_range(wave=wave, t=t, p=p, rh=rh, co2=co2)

    n = ciddor(wave, t, p, rh, co2)
    return wave * n


def _test_nist_ciddor_1():
    """Compare with NIST output.

    Values from http://emtoolbox.nist.gov/Wavelength/Ciddor.asp.

    fix t at 20, p at 101325, rh at 50
    """
    wave = [321.456, 500, 600.1234, 633.0, 700, 1000.987, 1500.8, 1700.0]
    nist_n = [1.000283543, 1.000273781, 1.000271818, 1.000271373,
              1.000270657, 1.000269038, 1.00026819, 1.000268041]
    nist_w = [321.364879, 499.863147, 599.96032, 632.828268,
              699.810591, 1000.717769, 1500.397608, 1699.544453]

    xv = rh2mole_fraction(50, 101325, 20)

    n = [ciddor_ri(i, 20, 101325, xv) for i in wave]
    wave_n = [vac2air(i, t=20, p=101325, rh=50.0) for i in wave]

    for i, j in zip(n, nist_n):
        assert abs(i - j) < 1e-8

    for i, j in zip(wave_n, nist_w):
        assert abs(i - j) < 1e-6

    n = [ciddor(i, 20, 101325, 50.0) for i in wave]
    wave_n = [vac2air(i, t=20, p=101325, rh=50.0) for i in wave]

    for i, j in zip(n, nist_n):
        assert abs(i - j) < 1e-8

    for i, j in zip(wave_n, nist_w):
        assert abs(i - j) < 1e-6


def _test_nist_ciddor_2():
    """Compare with NIST output.

    Values from http://emtoolbox.nist.gov/Wavelength/Ciddor.asp.

    fix wave at 633.0 p at 101325 rh at 50
    """
    t = [-20.0, 0.0, 20, 26.7982, 40.123, 60.45]
    nist_w = [632.800737, 632.815441, 632.828268, 632.832303, 632.839872,
              632.850953]
    nist_n = [1.00031489, 1.000291647, 1.000271373, 1.000264994, 1.000253031,
              1.000235516]

    xv = [rh2mole_fraction(50, 101325, i) for i in t]
    n = [ciddor_ri(633.0, i, 101325, j) for i, j in zip(t, xv)]

    wave_n = [vac2air(633.0, i, 101325, 50) for i in t]

    for i, j in zip(n, nist_n):
        assert abs(i - j) < 1e-8

    for i, j in zip(wave_n, nist_w):
        assert abs(i - j) < 1e-6


def _test_nist_ciddor_3():
    """Compare with NIST output.

    Values from http://emtoolbox.nist.gov/Wavelength/Ciddor.asp.

    fix wave at 633.0, t at 20, rh at 50. vary p
    """
    p = [1000 * i for i in [10, 50.123, 100.1234, 140.0]]

    nist_n = [1.000026385, 1.000133999, 1.000268148, 1.000375169]
    nist_w = [632.983299, 632.91519, 632.830308, 632.762607]

    xv = [rh2mole_fraction(50, i, 20) for i in p]
    n = [ciddor_ri(633.0, 20, i, j) for i, j in zip(p, xv)]

    wave_n = [vac2air(633.0, 20, i, 50) for i in p]

    for i, j in zip(n, nist_n):
        assert abs(i - j) < 1e-8

    for i, j in zip(wave_n, nist_w):
        assert abs(i - j) < 1e-6


def _test_nist_ciddor_4():
    """Compare with NIST output.

    Values from http://emtoolbox.nist.gov/Wavelength/Ciddor.asp.

    fix wave at 633.0, t at 20, p at 101325, vary rh.
    """
    rh = [0.0, 20.123, 40, 50.9876, 70, 90.7432, 100.0]
    nist_n = [1.0002718, 1.000271627, 1.000271458, 1.000271364,
              1.000271203, 1.000271027, 1.000270949]
    nist_w = [632.827997, 632.828106, 632.828214, 632.828273,
              632.828375, 632.828486, 632.828535]

    xv = [rh2mole_fraction(i, 101325, 20) for i in rh]
    n = [ciddor_ri(633.0, 20, 101325, j) for j in xv]

    wave_n = [vac2air(633.0, 20, 101325, i) for i in rh]

    for i, j in zip(n, nist_n):
        assert abs(i - j) < 1e-8

    for i, j in zip(wave_n, nist_w):
        assert abs(i - j) < 1e-6


def _test_air2vac():
    """Test reversibility with vac2air."""
    wave = [321.456, 500, 600.1234, 633.0, 700, 1000.987, 1500.8, 1700.0]
    wave_o = [air2vac(vac2air(i)) for i in wave]

    for i, j in zip(wave, wave_o):
        assert abs(i - j) < 1e-5


def _test_idlastro():
    # Using IDLASTRO downloaded on 2011/10/07. The vac2air.pro uses a
    # formulation of the Ciddor equation. Previous versions used a
    # different equation.

    # The REVISION HISTORY from the vac2air.pro file is:
    # ; REVISION HISTORY
    # ;       Written W. Landsman                November 1991
    # ;       Use Ciddor (1996) formula for better accuracy in the infrared
    # ;           Added optional output vector, W Landsman Mar 2011
    # ;       Iterate for better precision W.L./D. Schlegel  Mar 2011

    # The REVISION HISTORY from air2vac.pro file is:
    # ; REVISION HISTORY
    # ;	Written, D. Lindler 1982
    # ;	Documentation W. Landsman  Feb. 1989
    # ;       Use Ciddor (1996) formula for better accuracy in the infrared
    # ;           Added optional output vector, W Landsman Mar 2011

    # Velocity errors in m/s for different wave length errors, at
    # different wave lengths.
    # >>> 1e-5/330.0 * 299792458
    # 9.0846199393939404
    # >>> 1e-5/200.0 * 299792458
    # 14.989622900000001
    # >>> 1e-5/1000.0 * 299792458
    # 2.9979245800000003

    # nm
    wave = [200.0, 300.0, 500.0, 800.0, 1200.0, 1600.0, 1700.0]

    # angstrom
    wave_idl_vactoair = [1999.3526550081103323,
                         2999.1255923046301177,
                         4998.6055889614663101, 7997.8003315140686027,
                         11996.7167708424640296, 15995.6298776736693981,
                         16995.3579139663052047]
    wave_vac2air = [vac2air(i, t=15, rh=0) for i in wave]

    # values in wave_idl_vactoair was fed to airtovac idl procedure.
    wave_idl_airtovac = [1999.3526550081103323,
                        3000.0000371189012185,
                        5000.0000183785432455,
                        8000.0000108292333607,
                        12000.0000070745754783,
                        16000.0000052688483265,
                        17000.0000049538284657]
    # Have to convert angstrom to nm.
    wave_air2vac = [air2vac(i / 10.0, t=15, rh=0)
                    for i in wave_idl_vactoair]

    for i, j in zip(wave_vac2air, wave_idl_vactoair):
        # Convert nm to angstrom.
        assert abs(i - j / 10.0) < 1e-4

    # IDL code ignores values under 2000 angstrom.
    for i, j in zip(wave_air2vac[1:], wave_idl_airtovac[1:]):
        # Convert nm to angstrom.
        assert abs(i - j / 10.0) < 1e-4

    #return wave_idl_vactoair, wave_idl_airtovac


def _run_tests():
    import sys
    m = sys.modules[__name__]
    x = dir(m)
    tests = [getattr(m, i) for i in x if i.startswith("_test")]
    for test in tests:
        print(test.__name__)
        test()
