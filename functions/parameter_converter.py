import astropy.units as u
import astropy.constants as cst
import numpy as np
from edibles.functions.make_grid import make_grid



def param_convert(params):
    '''
    DESCRIPTION:
        Function to convert voigt parameteres between 
        mathematical version and astronomical version

    INPUT:
        params:
            cent     [float]     angstroms   central wavelength
            b_eff    [float]     km/s        velocity broadening
            Gamma    [float]     1/s         Damping constant
            scaling  [float]     ?

        method:




    OUTPUT:
        converted_params: 

    '''

    cent, b_eff, Gamma, scaling = params

    cent = cent * u.AA
    b_eff = b_eff * u.km/ u.s

    Gamma = Gamma * 1 / u.s



    # data

    # data_Hz = data.to(u.Hz, equivalencies=u.spectral()).value


    # cent

    nu0 = cent.to(u.Hz, equivalencies=u.spectral())



    # b_eff


        # b_eff -> delta_nu_D -> sigma -> alpha

    delta_nu_D = nu0 * b_eff / cst.c.to('km/s')  # freq * km/s / km/s
    sigma = delta_nu_D / np.sqrt(2)
    alpha_Hz = sigma * np.sqrt(2 * np.log(2))
    alpha = alpha_Hz * cst.c / (nu0**2)    # Hz * m/s  / Hz^2
    alpha = alpha.decompose().to(u.AA)

    # Gamma

    gam = Gamma / (4 * np.pi) * 1e-10





    converted_params = (cent.value, alpha.value, gam.value, scaling)

    return converted_params



if __name__ == "__main__":


    alpha = 0.0576265588185308
    gamma = 0.00048255778745462673
    delta_v = 1000
    x_min = 5977
    x_max = 5983
    cent = 5980
    scaling = 1.0
    n_piece = 3

    b_eff=3.47
    Gamma=6.064e7

    R = cst.c.value / delta_v
    x_nonbroad = make_grid(x_min, x_max, resolution=R)
    wave = np.array(x_nonbroad)

    params = (cent, b_eff, Gamma, scaling)

    new_params = param_convert(wave, params, method=0)
    print(new_params)