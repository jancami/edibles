from edibles.functions.voigtMathematical import voigt_math
from edibles.functions.parameter_converter import param_convert


def voigt_astro(x, cent, b_eff, Gamma, scaling=1.0):
    '''
    INPUT:
    x:         [ndarray]    Wavelength grid
    cent:      [float]      Central wavelength
    b_eff:     [float]      Velocity width [km/s]
    Gamma:     [float]      Lorentzian HWHM component * 4pi
    scaling: [float]      Height of the line

    OUTPUT:
    y:     [ndarray]    Voigt profile
    '''


    params = (cent, b_eff, Gamma, scaling)
    new_params = param_convert(params)


    cent, alpha, gam, scaling = new_params

    y = voigt_math(x, cent, alpha, gam, scaling)

    return y
