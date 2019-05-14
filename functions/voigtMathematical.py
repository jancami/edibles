import numpy as np
from scipy.special import wofz


def voigt_math(x, cent, alpha, gamma, amplitude=1.0):
    """
    Return the Voigt line shape centered at cent with Lorentzian component HWHM gamma
    and Gaussian component HWHM alpha.

    Input:  [Angstroms]

    x:         [ndarray]  Data grid
    cent:      [float]    Peak of the Voigt profile
    amplitude: [float]    Height of the peak
    alpha:     [float]    Gaussian HWHM component
    gamma:     [float]    Lorentzian HWHM component
    norm:      [bool]     area under peak equal to 1

    OUTPUT:

    y:         [ndarray]  Flux data for given inputs

    """

    sigma = alpha / np.sqrt(2 * np.log(2))



    z = (x - cent + 1j*gamma)/sigma/np.sqrt(2)
    y = amplitude * wofz(z).real / (sigma*np.sqrt(2*np.pi))
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
    amplitude = 1.

    b_eff=[3.47]
    Gamma=[6.064e7]


    # generate wavelength grid with resolving power delta_v (R = c/delta_v)
    R = cst.c.value / delta_v
    wave = make_grid(x_min, x_max, resolution=R)


    flux_norm = voigt_math(wave, cent, alpha, gamma, amplitude)

    plt.plot(wave, flux_norm, markersize='1', label='Data')
    plt.show()