import numpy as np
from scipy.special import wofz
import astropy.constants as cst


class Voigt:
    """
    An object with different vesions of the Voigt function.


    Attributes
    ----------
    type : str
        type of voigt function

    """

    def __init__(self):
        self


    def voigtMath(self, x, alpha, gamma):
        """
        Function to return the Voigt line shape centered at cent with Lorentzian component HWHM gamma
        and Gaussian component HWHM alpha.

        Creates a Voigt line profile using the scipy.special.wofz, which returns 
        the value of the Faddeeva function. 

        WARNING
        scipy.special.wofz is not compaible with np.float128 type parameters. 

        Input:
        ----------
        x : float64
            Dimensionless point/array
        alpha : float64
            Gaussian HWHM component
        gamma : float64
            Lorentzian HWHM component

        Output:
        -------
        ndarray
            Flux array for given input

        """

        sigma = alpha / np.sqrt(2 * np.log(2))

        return (
            np.real(wofz((x + 1j * gamma) / sigma / np.sqrt(2)))
            / sigma
            / np.sqrt(2 * np.pi)
        )


    def voigtOpticalDepth(self, lam, lam_0, b, d, Nf=1.0):
        """
        Converts parameters to make proper call to voigtMath


        Input:
        ----------
        lam : float64
            Wavelength grid
        lam_0 : float64
            Central wavelength
        b : float64
            Gaussian standard deviation
        d : float64
            Damping parameter
        Nf : float64
            Scaling parameter, default = 1.0

        Output:
        -------
        ndarray
            Optical depth for given input

        """

        # convert lam & lam_0 to x
        x = lam - lam_0

        # convert b to sigma, then alpha
        sigma = b * lam_0 / cst.c.to("km/s").value
        alpha = sigma * np.sqrt(2.0 * np.log(2.0))

        # convert d to gamma -- [ depends on what units we want to use ]

        # Currently, we are using the Lorentzian HWHM. This can easily be changed...
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        gamma = d
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        # create y data from voigtMath
        y = self.voigtMath(x, alpha, gamma)

        # Calculate tau_0
        tau_0 = (
            np.pi
            * (cst.e.esu.value) ** 2
            * Nf
            * (1e-8 * lam_0) ** 2
            / (cst.m_e.to("g").value * (cst.c.to("cm/s").value) ** 2)
        )  # cm
        # Convert cm to angstroms
        tau_0 *= 1e8

        # Calculate tau
        tau = tau_0 * y

        # return scaled & shifted data
        return tau


    def voigtAbsorptionLine(self, lam, lam_0, b, d, tau_0=0.1, N=None, f=None):
        """
        Function that takes in physical parameters and returns an absorption line. 

        Choose either a tau_0 parameter, or N and f together. Default is tau_0.

        Input:
        -----
        lam : float64
            Wavelength grid
        lam_0 : float64
            Central wavelength
        b : float64
            Gaussian standard deviation
        d : float64
            Damping parameter
        N : float64
            Column density
        f : float64
            Oscillator strength
        tau_0 : float64
            Optical depth at center of line

        Output:
        -------
        ndarray
            flux array of light transmission

        """

        # create Nf
        if (N is not None) and (f is not None):
            Nf = N * f
        else:
            Nf = (
                (tau_0 * 1e-8)
                * cst.m_e.to("g").value
                * (cst.c.to("cm/s").value) ** 2
                / (np.pi * (cst.e.esu.value) ** 2 * (1e-8 * lam_0) ** 2)
            )

        tau = self.voigtOpticalDepth(lam, lam_0, b, d, Nf)

        transmission = np.exp(-tau)

        return transmission


if __name__ == "__main__":

    from edibles.edibles.functions.make_grid import make_grid
    import matplotlib.pyplot as plt

    # set params
    alpha = 0.0576265588185308
    gamma = 0.00048255778745462673
    delta_v = 1000.0
    x_min = 5888.0
    x_max = 5892.0

    R = cst.c.value / delta_v
    lam = make_grid(x_min, x_max, resolution=R)

    lam_0 = 5890.0
    b = 3
    d = 0.00048255778745462673
    Nf = 28747080.71319038
    tau_0 = 0.05

    V = Voigt()

    x = lam - lam_0
    one = V.voigtMath(x=x, alpha=alpha, gamma=gamma)
    plt.plot(x, one)
    plt.show()

    two = V.voigtOpticalDepth(lam=lam, lam_0=lam_0, b=b, d=d, Nf=Nf)
    plt.plot(lam, two)
    plt.show()

    three = V.voigtAbsorptionLine(lam=lam, lam_0=lam_0, b=b, d=d, tau_0=tau_0)
    plt.plot(lam, three)

    plt.show()
