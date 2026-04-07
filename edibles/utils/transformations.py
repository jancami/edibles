"""
Functions for spectrum and magnitude transformations.
"""
import numpy as np
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.time import Time
from astropy import units as u
import pandas as pd

c_light = 299792.458  # speed of light in km/s


def wavenumber_to_rv(wavenumber: np.array, ref_wavenumber: float) -> np.array:
    """
    Transforms np.array of wavenumbers to radial velocities w.r.t. reference wavenumber.

    Parameters
    ----------
    wavenumber : np.array
        Wavenumber array
    ref_wavenumber : float
        Reference wavenumber, e.g. lab wavenumber or wavenumber of line center

    Returns
    -------
    np.array
        Radial velocity array
    """
    wavenumber = np.array(wavenumber)

    z = ref_wavenumber / wavenumber - 1
    beta = (z ** 2 + 2 * z) / (z ** 2 + 2 * z + 2)  # relativistic velocity divided by speed of light

    return np.array(beta * c_light)


def rv_to_wavenumber(rv: np.array, ref_wavenumber: float) -> np.array:
    """
    Transforms np.array of radial velocities to wavenumbers w.r.t. reference wavenumber.

    Parameters
    ----------
    rv : np.array
        Radial velocity array
    ref_wavenumber : float
        Reference wavenumber, e.g. lab wavenumber or wavenumber of line center

    Returns
    -------
    np.array
        Wavenumber array
    """
    rv = np.array(rv)
    beta = rv / c_light  # relativistic velocity divided by speed of light
    z = np.sqrt((1 + beta) / (1 - beta)) - 1  # relativistic doppler shift

    return np.array(ref_wavenumber / (1 + z))


def wavelength_to_rv(wavelength: np.array, ref_wavelength: float) -> np.array:
    """
    Transforms np.array of wavelengths to radial velocities w.r.t. reference wavelength.
    Calculated also for relativistic velocities.

    Parameters
    ----------
    wavelength : np.array
        Wavelength array
    ref_wavelength : float
        Reference wavelength, e.g. lab wavelength or wavelength of line center

    Returns
    -------
    np.array
        Radial velocity array
    """
    wavelength = np.array(wavelength)
    z = wavelength / ref_wavelength - 1  # Doppler shift between observed wavelength and reference wavelength
    beta = (z ** 2 + 2 * z) / (z ** 2 + 2 * z + 2)  # relativistic velocity divided by speed of light

    return np.array(beta * c_light)


def rv_to_wavelength(rv: np.array, ref_wavelength: float) -> np.array:
    """
    Transforms np.array of radial velocities to wavelengths w.r.t. reference wavelength.
    Calculated also for relativistic velocities.

    Parameters
    ----------
    rv : np.array
        Radial velocity array
    ref_wavelength : float
        Reference wavelength, e.g. lab wavelength or wavelength of line center

    Returns
    -------
    np.array
        Wavelength array
    """
    rv = np.array(rv)
    beta = rv / c_light  # relativistic velocity divided by speed of light
    gamma = np.sqrt((1 - beta) / (1 + beta))  # relativistic Lorentz factor
    return np.array(ref_wavelength / gamma)


def bary_corr(wavelength: np.array, star_name: str = None, obs_name: str = None, obs_location: list = None,
              obs_time=None,
              time_format: str = None, return_bc_rv: bool = False, silent=False) -> np.array:
    """
    Calculates barycentric correction of a wavelength array and shifts the wavelength from the observer frame
    to the barycentric rest frame.

    Parameters
    ----------
    wavelength : 1D-np.array
        Wavelengths in observer rest frame.

    star_name : str
        Simbad name of the target.

    obs_name : str
        Name of the observatory.

    obs_location : list(float)
        A list containing [lon, lat, heigth] of the observatory as floats. Only needed if no obs_name is provided.

    obs_time : float or string
        Observation date.

    time_format : str
        Format of the observation time.
        See available formats for **astropy.time.Time** object.
        https://docs.astropy.org/en/stable/time/index.html

    return_bc_rv : bool
        If set to True, the functions returns additionally the barycentric radial velocity correction in km/s.

    silent : bool
        If set to True, the function will not print the radial velocity correction. (Default: False)

    Returns
    -------
    np.array(barycentric_wavelength)
        1D-Wavelength array in barycentric rest frame.
    """
    wavelength = np.array(wavelength)
    # check if kwargs are missing, i.e. still None
    if star_name is None:
        raise ValueError('Missing star_name.')
    if obs_time is None:
        raise ValueError('Missing obs_time.')
    if time_format is None:
        raise ValueError('Missing time_format.')

    coord = SkyCoord.from_name(star_name)  # get coordinates of target
    obs_time = Time(obs_time, format=time_format)  # get the observation time

    if obs_name is None:
        if obs_location is None:
            raise ValueError('Missing obs_location.')
        else:
            obs_loc = EarthLocation.from_geodetic(*obs_location)
    else:
        obs_loc = EarthLocation.of_site(obs_name)  # get location of observation site

    # calculate barycentric correction velocity
    bary_corr_unit_less = coord.radial_velocity_correction(location=obs_loc, obstime=obs_time)
    z = bary_corr_unit_less.to(u.cds.c) / u.cds.c  # transform radial velocity into units of speed of light
    wavelength = wavelength * (1 + z)  # shift the spectrum

    if not silent:
        print('Radial velocity correction: ', z * c_light, ' km/s')  # print RV correction

    if return_bc_rv:
        return wavelength, z * c_light
    else:
        return wavelength


def lsr_corr(wavelength: np.array, star_name: str = None, ref='schoenrich2010', return_lsr_corr_rv: bool = False) \
        -> np.array:
    """
    Calculates LSR correction of a wavelength array and shifts the wavelength from the barycentric rest frame
    to the LSR.

    Parameters
    ----------
    wavelength : 1D-np.array
        Wavelengths in barycentric rest frame.

    star_name : str
        Simbad name of the target.

    ref : str
        Reference for the solar peculiar velocity.

    return_lsr_corr_rv : bool
        If set to True, the functions returns additionally the LSR radial velocity correction in km/s.

    Returns
    -------
    np.array(LSR_wavelength)
        1D-Wavelength array in LSR.
    """

    wavelength = np.array(wavelength)
    coord = SkyCoord.from_name(star_name)
    l_gal, b_gal = coord.galactic.l.radian, coord.galactic.b.radian

    if ref == "zbinden":
        u_lsr, v_lsr, w_lsr = 11.9, 11.9, 6.4  # solar peculiar motion in km/s (Zbinden)
    elif ref == "schoenrich2010":
        u_lsr, v_lsr, w_lsr = 11.1, 12.24, 7.25  # solar peculiar motion in km/s (Schönrich 2010)
    else:
        print("No peculiar velocity specified!")

    lsr_vector = np.array([u_lsr, v_lsr, w_lsr])  # solar peculiar motion as vector

    # View direction in the Galactic XYZ-coordinates as vector
    view_direction = np.array([np.cos(l_gal) * np.cos(b_gal), np.sin(l_gal) * np.cos(b_gal), np.sin(b_gal)])

    # Calculate the LSR correction along the view direction
    lsr_corr_rv = np.dot(lsr_vector, view_direction)
    z = lsr_corr_rv / c_light

    wavelength = wavelength * (1 + z)  # non-relativistic doppler shift

    print('LSR correction: ', lsr_corr_rv, ' km/s')  # print RV correction

    if return_lsr_corr_rv:
        return wavelength, lsr_corr_rv
    else:
        return wavelength


def angstrom_vac_to_air(angstrom: np.array) -> np.array:
    """
    For the conversion of λ_vac to λ_air the IAU standard formula provided by Donald Morton [RD27] λ_air = λ_vac/n
    is used.

    Parameters
    ----------
    angstrom : np.array
        Vacuum wavelength 1D-array in Ångströms.

    Returns
    -------
    np.array
        Air wavelength 1D-array in Ångströms.
    """
    angstrom = np.array(angstrom)
    s = 10 ** 4 / angstrom
    n = 1 + 0.0000834254 + 0.02406147 / (130 - s ** 2) + 0.00015998 / (38.9 - s ** 2)

    angstrom = angstrom / n

    return angstrom


def angstrom_air_to_vac(angstrom: np.array) -> np.array:
    """
    The reverse transform λ_air to λ_vac was derived by N.Piskunov.
    (see: https://www.astro.uu.se/valdwiki/Air-to-vacuum%20conversion)

    Parameters
    ----------
    angstrom : np.array
        Air wavelength 1D-array in Ångströms.

    Returns
    -------
    np.array
        Vacuum wavelength 1D-array in Ångströms.
    """
    angstrom = np.array(angstrom)
    s = 10 ** 4 / angstrom
    n = 1 + 0.000083366242 + 0.0240892687 / (130.106592452 - s ** 2) + 0.00015997409 / (38.925687933 - s ** 2)

    angstrom = angstrom * n

    return angstrom


def angstrom_to_wavenumber(angstrom: np.array, air_to_vac=True) -> np.array:
    """
    Transforms np.array of wave lengths in Ångströms to array of wavenumbers.

    Parameters
    ----------
    angstrom : np.array
        Wavelength array in Ångströms
    air_to_vac : bool
        If True, wavelengths are assumed in air and get transform into vacuum before wavenumber transformation.
        Default: True

    Returns
    -------
    np.array
        Wavenumber array
    """
    angstrom = np.array(angstrom)
    if air_to_vac:
        angstrom = angstrom_air_to_vac(angstrom)
    return np.array(10 ** 8 / angstrom)


def wavenumber_to_angstrom(wavenumber: np.array, vac_to_air=True) -> np.array:
    """
    Transforms np.array of wavenumbers to array of wave lengths in angstrom.

    Parameters
    ----------
    wavenumber : np.array
        Wavenumber array
    vac_to_air : bool
        If True, wavelengths are transformed to air.
        Default: True

    Returns
    -------
    np.array
        Ångströms array
    """
    wavenumber = np.array(wavenumber)
    angstrom = 10 ** 8 / wavenumber
    if vac_to_air:
        angstrom = angstrom_vac_to_air(angstrom)
    return angstrom


def doppler_shift_wl(wavelength: np.array, rv: float) -> np.array:
    """
    Applies Doppler shift to wavelength array.

    Parameters
    ----------
    wavelength : np.array
        Wavelength array
    rv : float
        Radial velocity in km/s

    Returns
    -------
    np.array
        Wavelength array
    """
    wavelength = np.array(wavelength)
    beta = rv / c_light  # relativistic velocity divided by speed of light
    gamma = np.sqrt((1 - beta) / (1 + beta))  # relativistic Lorentz factor
    return np.array(wavelength / gamma)


def doppler_shift_wn(wavenumber: np.array, rv: float) -> np.array:
    """
    Applies Doppler shift to wavelength array.

    Parameters
    ----------
    wavenumber : np.array
        Wavenumber array
    rv : float
        Radial velocity in km/s

    Returns
    -------
    np.array
        Wavenumber array
    """
    wavenumber = np.array(wavenumber)
    beta = rv / c_light  # relativistic velocity divided by speed of light
    gamma = np.sqrt((1 - beta) / (1 + beta))  # relativistic Lorentz factor
    return np.array(wavenumber * gamma)
