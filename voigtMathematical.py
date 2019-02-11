import numpy as np
from scipy.special import wofz
# import matplotlib.pyplot as plt


def voigt_math(x, alpha, gamma, x_pts=None, y_pts=None, cont=False):
    """
    Return the Voigt line shape at x with Lorentzian component HWHM gamma
    and Gaussian component HWHM alpha.

    Input:
    x:      [ndarray]  Data grid
    alpha:  [float]    Gaussian HWHM component
    gamma:  [float]    Lorentzian HWHM component
    x_pts:  [ndarray]  x-components of initial continuum guess
    y_pts:  [ndarray]  y-components of initial continuum guess
    cont:   [bool]     Use a continuum term? Default is False

    """
    sigma = alpha / np.sqrt(2 * np.log(2))


    if cont is True:
        # call continuum function here?
        print("Do something")

    return np.real(wofz((x + 1j*gamma)/sigma/np.sqrt(2))) / sigma/np.sqrt(2*np.pi)



# alpha, gamma = 0.1, 0.1
# x = np.linspace(-0.8,0.8,1000)
# y = voigt_math(x, alpha, gamma)

# plt.plot(x, y, label='Voigt')
# # plt.xlim(-0.8,0.8)
# plt.legend()
# plt.show()