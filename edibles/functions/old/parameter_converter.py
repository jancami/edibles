import astropy.units as u
import astropy.constants as cst
import numpy as np


def param_convert(params):
    """
    Function to convert voigt parameteres from astronomical
    to mathematical version.

    Parameters
    ----------
    cent : float
        central wavelength
    b_eff : float
        velocity broadening
    Gamma : float
        Damping constant
    scaling : float
        tau_0 ???

    Returns
    -------
    list
        converted parameters 

    """

    cent, b_eff, Gamma, scaling = params

    cent = cent * u.AA
    b_eff = b_eff * u.km / u.s
    Gamma = Gamma * 1 / u.s

    # cent
    nu0 = cent.to(u.Hz, equivalencies=u.spectral())
    # b_eff
    delta_nu_D = nu0 * b_eff / cst.c.to("km/s")  # freq * km/s / km/s
    sigma = delta_nu_D / np.sqrt(2)
    alpha_Hz = sigma * np.sqrt(2 * np.log(2))
    alpha = alpha_Hz * cst.c / (nu0 ** 2)  # Hz * m/s  / Hz^2
    alpha = alpha.decompose().to(u.AA)
    # Gamma
    gam = Gamma / (4 * np.pi) * 1e-10

    converted_params = (cent.value, alpha.value, gam.value, scaling)
    return converted_params


if __name__ == "__main__":

    # Proper parameters
    alpha = 0.0576265588185308
    gamma = 0.00048255778745462673

    # Given parameters
    cent = 5980
    b_eff = 3.47
    Gamma = 6.064e7
    scaling = 1.0

    params = (cent, b_eff, Gamma, scaling)

    new_params = param_convert(params)
    print("Parameters should be:")
    print("cent|alpha             |gamma                 |scaling")
    print(cent, alpha, gamma, scaling)
    print()
    print("Parameters are:")
    print(new_params)
