import numpy as np
import matplotlib.pyplot as plt
from edibles.projects.DIBCode_Empirical.composite_model import composite_model
from numpy.polynomial import chebyshev
def plot_final_with_components(x, data, params, figpath):
    """
    Creates the final plot of the components.
    Args:
        x: Independent variable data points (wavelengths) as a NumPy array.
        data: List of coadded flux as a function of wavelength for each sightline
        params: Parameters of the model components
    Returns:
        Final plot of the components.
    """
    total_model = composite_model(x, params)
    cheb_coeffs = []
    i = 0
    while f'c{i}' in params:
        cheb_coeffs.append(params[f'c{i}'].value)
        i += 1

    x_norm = 2 * (x - x.min()) / (x.max() - x.min()) - 1
    continuum = chebyshev.chebval(x_norm, cheb_coeffs)

    plt.figure(figsize=(12, 6))
    plt.plot(x, data, "k.", alpha=0.4, label="Data")
    plt.plot(x, total_model, "r-", lw=2, label="Total Model")
    plt.plot(x, continuum, "g-", lw=1.5, label="Continuum")

    # --- Plot Gaussians on continuum ---
    i = 0
    while f'g{i}_amp' in params:
        amp = params[f'g{i}_amp'].value
        cen = params[f'g{i}_cen'].value
        wid = params[f'g{i}_wid'].value

        gauss = amp * np.exp(-((x - cen) / wid) ** 2)
        plt.plot(x, continuum + gauss, "--", lw=1.5, label=f"Gaussian {i}")
        i += 1

        # --- Plot Lorentzians on continuum ---
    i = 0
    while f'l{i}_amp' in params:
        amp = params[f'l{i}_amp'].value
        cen = params[f'l{i}_cen'].value
        wid = params[f'l{i}_wid'].value

        lor = amp * wid ** 2 / ((x - cen) ** 2 + wid ** 2)
        plt.plot(x, continuum + lor, "--", lw=1.5, label=f"Lorentzian {i}")
        i += 1

    plt.legend(fontsize=8, ncol=2)
    plt.title("Final Model with All Components (Continuum + Each Line)")
    plt.xlabel("Wavelength")
    plt.ylabel("Flux")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{figpath}/Final_Model_With_Components.png")
    plt.show()

