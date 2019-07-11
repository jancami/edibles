import numpy as np
from scipy.special import wofz
import astropy.constants as cst


def voigtNorm(x, alpha, gamma):
    """
    Return the Voigt line shape centered at cent with Lorentzian component HWHM gamma
    and Gaussian component HWHM alpha.

    INPUT:

    x:        [float64]    dimensionless point/array
    alpha:    [float64]    Gaussian HWHM component
    gamma:    [float64]    Lorentzian HWHM component

    OUTPUT:

    y:        [float64]  Flux for given input

    """

    sigma = alpha / np.sqrt(2. * np.log(2.))

    z = (x + 1j*gamma)/sigma/np.sqrt(2.)

    y = wofz(z).real / (sigma*np.sqrt(2.*np.pi))

    return y


def voigtOpticalDepth(lam, lam_0, b, d, Nf=1.0):
    """
    Converts parameters to make proper call to voigtNorm

    INPUT:

    lam:      [float64]  (Angstroms)  Wavelength grid
    lam_0:    [float64]  (Angstroms)  Central wavelength
    b:        [float64]  (km/s)       Gaussian standard deviation
    d:        [float64]  (units)      Damping parameter
    Nf:       [float64]  (units)      Scaling parameter, default = 1.0

    OUTPUT:

    tau:      [float64]  (units)  Optical depth for given input

    """

    # convert lam & lam_0 to x

    x = lam - lam_0

    # convert b to alpha    -->    b = sigma / lam_0 * cst.c.to('km/s')

    sigma = b * lam_0 / cst.c.to('km/s').value
    alpha = sigma * np.sqrt(2. * np.log(2.))

    # convert d to gamma -- ??? is d: [ gamma OR Gamma OR a ] ???

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gamma = d  
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # create y data from voigtNorm

    y = voigtNorm(x, alpha, gamma)

    # shift data to desired wavelength location -- I don't think this is a real step

    # multiply by tau_0

    tau_0 = np.pi * (cst.e.esu.value)**2 * Nf*1e8*lam_0 / (cst.m_e.to('g').value*(cst.c.to('cm/s').value)**2)

    tau = tau_0 * y

    # return scaled & shifted data
    return tau


def voigtAbsorptionLine(lam, lam_0, b, d, tau_0, N=None, f=None):
    """
    Converts parameters to make proper call to voigtOpticalDepth

    INPUT:
    lam:          [float64]  (Angstroms)  Wavelength grid
    lam_0:        [float64]  (Angstroms)  Central wavelength
    b:            [float64]  (km/s)       Gaussian standard deviation
    d:            [float64]  (units)      Damping parameter
        ======================= EITHER ========================
        N:        [float64]  (units)      Column density
        f:        [float64]  (units)      Oscillator strength
        ========================  OR  ========================
        tau_0:    [float64]  (units)      Scaling parameter, default = 0.1


    OUTPUT:
    transmission:    [float64]  (units)  

    """

    # create Nf
    if (N is not None) and (f is not None):
        Nf = N * f
    else:
        Nf = tau_0 * cst.m_e.to('g').value * (cst.c.to('cm/s').value)**2 / (np.pi * (cst.e.esu.value)**2 * 1e8*lam_0)

    tau = voigtOpticalDepth(lam, lam_0, b, d, Nf)

    transmission = np.exp(-tau)

    return transmission


if __name__ == "__main__":

    from edibles.fit.make_grid import make_grid
    import matplotlib.pyplot as plt

        # set params
    alpha = 0.0576265588185308
    gamma = 0.00048255778745462673
    delta_v = 1000.
    x_min = 5888.
    x_max = 5892.

    R = cst.c.value / delta_v
    lam = make_grid(x_min, x_max, resolution=R)

    lam_0 = 5890.
    b = 3
    d = 0.00048255778745462673
    Nf = 28747080.71319038
    tau_0 = 0.05

    x = lam - lam_0
    one = voigtNorm(x, alpha, gamma)
    plt.plot(x, one)
    plt.show()

    two = voigtOpticalDepth(lam=lam, lam_0=lam_0, b=b, d=d, Nf=Nf)
    plt.plot(lam, two)
    plt.show()

    three = voigtAbsorptionLine(lam=lam, lam_0=lam_0, b=b, d=d, tau_0=tau_0)
    plt.plot(lam, three)

    plt.show()

