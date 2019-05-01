from __future__ import print_function
from edibles.functions.astro_wrapper import voigt_astro
import numpy as np

def multi_voigt(x, cents, b_effs, Gammas, n=1):
    '''
    Function that takes in wavelength grid and parameters for each absorption line,
    and outputs a spectrum

    INPUT:
    x        [ndarray]   input wavelength grid
    cents    [list]      guesses for each central peak
    b_effs   [list]      guesses for each b_eff
    Gammas   [list]      guesses for each Gamma
    n        [int]       number of lines

    OUTPUT:
    y        [ndarray]   spectrum from given parameters
    '''

    # Input checking
    assert (n == len(cents)), 'Incorrect number of central wavelengths!'
    assert (n == len(b_effs)), 'Incorrect number of b_effs!'
    assert (n == len(Gammas)), 'Incorrect number of Gammas!'

    # Create initial continuum flux values
    cont = np.zeros_like(x)
    y = cont

    for i in range(n):
        cent = cents[i]
        b_eff = b_effs[i]
        Gamma = Gammas[i]

        abs_line = voigt_astro(x, cent, b_eff, Gamma)
        y += abs_line

    return y


if __name__ == "__main__":
    from edibles.fit.make_grid import make_grid
    from astropy import constants as cst
    import matplotlib.pyplot as plt

    # set params

    delta_v = 1000
    x_min = 5978
    x_max = 5982
    cent = [5980, 5981]

    b_eff= [3.47, 10]
    Gamma= [6.064e7, 6.064e7]

    # generate wavelength grid with resolving power delta_v (R = c/delta_v)
    R = cst.c.value / delta_v
    wave = make_grid(x_min, x_max, resolution=R)

    flux = multi_voigt(wave, cents=cent, b_effs=b_eff, Gammas=Gamma, n=2)
    plt.plot(wave,flux)
    plt.show()
