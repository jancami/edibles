import numpy as np
from scipy.special import wofz


def voigt_math(x, cent, alpha, gamma):
    """
    Return the Voigt line shape centered at cent with Lorentzian component HWHM gamma
    and Gaussian component HWHM alpha.

    Input:  [Angstroms]

    x:         [ndarray]  Data grid
    cent:      [float]    Peak of the Voigt profile
    alpha:     [float]    Gaussian HWHM component
    gamma:     [float]    Lorentzian HWHM component

    OUTPUT:

    y:         [ndarray]  Flux data for given inputs

    """

    sigma = alpha / np.sqrt(2 * np.log(2))

    y = np.real(wofz(((x - cent) + 1j*gamma)/sigma/np.sqrt(2))) / sigma/np.sqrt(2*np.pi)

    return y






if __name__ == "__main__":

    from edibles.fit.make_grid import make_grid
    from astropy import constants as cst
    import matplotlib.pyplot as plt

        # set params
    alpha = 0.0576265588185308
    gamma = 0.00048255778745462673
    delta_v = 1000
    x_min = 5978
    x_max = 5982
    cent = [5980]

    b_eff=[3.47]
    Gamma=[6.064e7]


    # generate wavelength grid with resolving power delta_v (R = c/delta_v)
    R = cst.c.value / delta_v
    wave = make_grid(x_min, x_max, resolution=R)


    flux_norm = voigt_math(wave, cent, alpha, gamma)

    plt.plot(wave, flux_norm, markersize='1', label='Data')
    plt.show()