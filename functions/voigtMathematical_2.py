class Voigt1D(Fittable1DModel):
    """
    One dimensional model for the Voigt profile.

    Parameters
    ----------
    x_0 : float
        Position of the peak
    amplitude_L : float
        The Lorentzian amplitude
    fwhm_L : float
        The Lorentzian full width at half maximum
    fwhm_G : float
        The Gaussian full width at half maximum

    See Also
    --------
    Gaussian1D, Lorentz1D

    Notes
    -----
    Algorithm for the computation taken from
    McLean, A. B., Mitchell, C. E. J. & Swanston, D. M. Implementation of an
    efficient analytical approximation to the Voigt function for photoemission
    lineshape analysis. Journal of Electron Spectroscopy and Related Phenomena
    69, 125-132 (1994)

    Examples
    --------
    .. plot::
        :include-source:

        import numpy as np
        from astropy.modeling.models import Voigt1D
        import matplotlib.pyplot as plt

        plt.figure()
        x = np.arange(0, 10, 0.01)
        v1 = Voigt1D(x_0=5, amplitude_L=10, fwhm_L=0.5, fwhm_G=0.9)
        plt.plot(x, v1(x))
        plt.show()
    """

    x_0 = Parameter(default=0)
    amplitude_L = Parameter(default=1)
    fwhm_L = Parameter(default=2/np.pi)
    fwhm_G = Parameter(default=np.log(2))

    _abcd = np.array([
        [-1.2150, -1.3509, -1.2150, -1.3509],  # A
        [1.2359, 0.3786, -1.2359, -0.3786],    # B
        [-0.3085, 0.5906, -0.3085, 0.5906],    # C
        [0.0210, -1.1858, -0.0210, 1.1858]])   # D

# [docs]    @classmethod
    def evaluate(cls, x, x_0, amplitude_L, fwhm_L, fwhm_G):

        A, B, C, D = cls._abcd
        sqrt_ln2 = np.sqrt(np.log(2))
        X = (x - x_0) * 2 * sqrt_ln2 / fwhm_G
        X = np.atleast_1d(X)[..., np.newaxis]
        Y = fwhm_L * sqrt_ln2 / fwhm_G
        Y = np.atleast_1d(Y)[..., np.newaxis]

        V = np.sum((C * (Y - A) + D * (X - B))/(((Y - A) ** 2 + (X - B) ** 2)), axis=-1)

        return (fwhm_L * amplitude_L * np.sqrt(np.pi) * sqrt_ln2 / fwhm_G) * V


# [docs]    @classmethod
    def fit_deriv(cls, x, x_0, amplitude_L, fwhm_L, fwhm_G):

        A, B, C, D = cls._abcd
        sqrt_ln2 = np.sqrt(np.log(2))
        X = (x - x_0) * 2 * sqrt_ln2 / fwhm_G
        X = np.atleast_1d(X)[:, np.newaxis]
        Y = fwhm_L * sqrt_ln2 / fwhm_G
        Y = np.atleast_1d(Y)[:, np.newaxis]
        constant = fwhm_L * amplitude_L * np.sqrt(np.pi) * sqrt_ln2 / fwhm_G

        alpha = C * (Y - A) + D * (X - B)
        beta = (Y - A) ** 2 + (X - B) ** 2
        V = np.sum((alpha / beta), axis=-1)
        dVdx = np.sum((D/beta - 2 * (X - B) * alpha / np.square(beta)), axis=-1)
        dVdy = np.sum((C/beta - 2 * (Y - A) * alpha / np.square(beta)), axis=-1)

        dyda = [-constant * dVdx * 2 * sqrt_ln2 / fwhm_G,
                constant * V / amplitude_L,
                constant * (V / fwhm_L + dVdy * sqrt_ln2 / fwhm_G),
                -constant * (V + (sqrt_ln2 / fwhm_G) * (2 * (x - x_0) * dVdx + fwhm_L * dVdy)) / fwhm_G]
        return dyda


    @property
    def input_units(self):
        if self.x_0.unit is None:
            return None
        else:
            return {'x': self.x_0.unit}

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        return OrderedDict([('x_0', inputs_unit['x']),
                            ('fwhm_L', inputs_unit['x']),
                            ('fwhm_G', inputs_unit['x']),
                            ('amplitude_L', outputs_unit['y'])])