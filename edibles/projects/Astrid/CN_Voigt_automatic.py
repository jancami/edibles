import os
from edibles.utils.edibles_spectrum import EdiblesSpectrum
import numpy as np
import matplotlib.pyplot as plt
from voigt_profile import voigt_profile
import csv


# =========
# 📂 Carpeta datos
# =========

data_folder = "/Users/astridsaldana/Desktop/MITACS/EDIBLES/EDR5/orders/"


# =========
# Outputs
# =========

output_csv = "CN_results.csv"
plot_folder = "CN_plots/"

os.makedirs(plot_folder, exist_ok=True)


# =========
# Líneas CN
# =========

R0 = 3874.608
R1 = 3873.998
P1 = 3875.763


# =========
# Parámetros
# =========

sigma = 0.014
gamma = 0.004


results = []


# =========
# LOOP
# =========

for file in os.listdir(data_folder):

    # ✅ SOLO espectros útiles (CLAVE)
    if file.endswith(".fits") and "437nm_blue" in file:

        try:
            print(f"Procesando: {file}")

            # ✅ CORRECTO para EDIBLES
            sp = EdiblesSpectrum("orders/" + file)

            wave = sp.bary_wave
            flux = sp.flux

            flux_norm = flux / np.median(flux)

            mask = (wave > 3873.0) & (wave < 3876.0)
            wave = wave[mask]
            flux_norm = flux_norm[mask]

            # ✅ EXTRA SEGURIDAD
            if len(wave) == 0:
                print(f"Saltado (sin CN): {file}")
                continue

            continuum = np.median(flux_norm)

            tau0_main = -np.log(np.min(flux_norm) / continuum)

            profile_R0 = voigt_profile(wave - R0, sigma, gamma)
            profile_R1 = voigt_profile(wave - R1, sigma, gamma)
            profile_P1 = voigt_profile(wave - P1, sigma, gamma)

            profile_R0 /= np.max(profile_R0)
            profile_R1 /= np.max(profile_R1)
            profile_P1 /= np.max(profile_P1)

            tau_R0 = 0.9 * tau0_main
            tau_R1 = 0.08 * tau0_main
            tau_P1 = 0.06 * tau0_main

            tau_total = (
                tau_R0 * profile_R0 +
                tau_R1 * profile_R1 +
                tau_P1 * profile_P1
            )

            fit_flux = continuum * np.exp(-tau_total)

            results.append([file, tau_R0, tau_R1, tau_P1, sigma, gamma])

            plt.figure(figsize=(8,5))
            plt.plot(wave, flux_norm, label="Observed", color="blue")
            plt.plot(wave, fit_flux, label="Voigt CN", color="red")

            plt.xlim(3873, 3876)
            plt.legend()

            plot_name = os.path.join(plot_folder, file.replace(".fits", ".png"))
            plt.savefig(plot_name)
            plt.close()

        except Exception as e:
            print(f"Error con {file}: {e}")


# =========
# GUARDAR CSV
# =========

with open(output_csv, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["file", "tau_R0", "tau_R1", "tau_P1", "sigma", "gamma"])
    writer.writerows(results)


print("✅ Todo terminado!")
