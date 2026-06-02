import os
from edibles.utils.edibles_spectrum import EdiblesSpectrum
import numpy as np
import matplotlib.pyplot as plt
from voigt_profile import voigt_profile
from scipy.optimize import curve_fit
import csv


data_folder = "/Users/astridsaldana/Desktop/MITACS/EDIBLES/EDR5/orders/"
plot_folder = "CN_plots/"
output_csv = "CN_results.csv"

os.makedirs(plot_folder, exist_ok=True)

# línea principal CN
R0_rest = 3874.608

results = []


def voigt_model(wave, center, sigma, depth, continuum):
    profile = voigt_profile(wave - center, sigma, 0.004)
    profile /= np.max(profile)
    return continuum * np.exp(-depth * profile)


for file in os.listdir(data_folder):

    if file.endswith(".fits") and "437nm_blue" in file:

        try:
            print(f"Procesando: {file}")

            sp = EdiblesSpectrum("orders/" + file)

            wave = sp.wave
            flux = sp.flux

            flux_norm = flux / np.median(flux)

            mask = (wave > 3873.0) & (wave < 3876.0)
            wave = wave[mask]
            flux_norm = flux_norm[mask]

            if len(wave) == 0:
                continue

            continuum = np.median(flux_norm)

            # detectar mínimo real
            idx = np.argmin(flux_norm)
            center_guess = wave[idx]
            depth_guess = -np.log(np.min(flux_norm) / continuum)
            sigma_guess = 0.01

            # filtro: sin CN
            if depth_guess < 0.02:
                print(f"Sin CN: {file}")
                continue

            # ✅ FIT REAL
            popt, _ = curve_fit(
                voigt_model,
                wave,
                flux_norm,
                p0=[center_guess, sigma_guess, depth_guess, continuum],
                maxfev=5000
            )

            center_fit, sigma_fit, depth_fit, cont_fit = popt

            fit_flux = voigt_model(wave, *popt)

            results.append([file, center_fit, sigma_fit, depth_fit])

            # plot
            plt.figure(figsize=(8,5))
            plt.plot(wave, flux_norm, label="Observed", color="blue")
            plt.plot(wave, fit_flux, label="Voigt Fit", color="red")
            plt.legend()
            plt.xlim(3873, 3876)

            plt.savefig(f"{plot_folder}/{file}.png")
            plt.close()

        except Exception as e:
            print(f"Error con {file}: {e}")


# guardar CSV
with open(output_csv, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["file", "center", "sigma", "depth"])
    writer.writerows(results)


