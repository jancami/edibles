from edibles.utils.edibles_spectrum import EdiblesSpectrum
import numpy as np
import matplotlib.pyplot as plt
from voigt_profile import voigt_profile


# =========
# 1. Cargar espectro
# =========

filename = "orders/HD149757_2015-07-20T02:54:27_437nm_blue_O5.fits"
sp = EdiblesSpectrum(filename)

wave = sp.bary_wave
flux = sp.flux


# =========
# 2. Normalizar
# =========

flux_norm = flux / np.median(flux)


# =========
# 3. Seleccionar región CN
# =========

mask = (wave > 3873.0) & (wave < 3876.0)

wave = wave[mask]
flux_norm = flux_norm[mask]


# =========
# 4. Líneas CN (valores físicos)
# =========

R0 = 3874.608
R1 = 3873.998
P1 = 3875.763


# =========
# 5. Continuo y profundidad base
# =========

continuum = np.median(flux_norm)
tau0_main = -np.log(np.min(flux_norm) / continuum)


# =========
# 6. Parámetros físicos (de R0 ya afinados)
# =========

sigma = 0.014
gamma = 0.004


# =========
# 7. Perfiles Voigt individuales
# =========

profile_R0 = voigt_profile(wave - R0, sigma, gamma)
profile_R1 = voigt_profile(wave - R1, sigma, gamma)
profile_P1 = voigt_profile(wave - P1, sigma, gamma)

# Normalizar cada perfil
profile_R0 = profile_R0 / np.max(profile_R0)
profile_R1 = profile_R1 / np.max(profile_R1)
profile_P1 = profile_P1 / np.max(profile_P1)


# =========
# 8. Profundidades relativas (IMPORTANTE)
# =========

# 🔬 Estas son aproximadas — puedes ajustarlas
tau_R0 = tau0_main
tau_R1 = 0.09 * tau0_main
tau_P1 = 0.07 * tau0_main




# =========
# 9. Modelo multilínea (suma de τ)
# =========

tau_total = (
    tau_R0 * profile_R0
    + tau_R1 * profile_R1
    + tau_P1 * profile_P1
)

fit_flux = continuum * np.exp(-tau_total)


# =========
# 10. Plot final
# =========

plt.figure(figsize=(10,6))

plt.plot(wave, flux_norm, label="Observed Spectrum", color="blue")
plt.plot(wave, fit_flux, label="CN Voigt Model (R0+R1+P1)", color="red")

plt.xlim(3873.0, 3876.0)

plt.xlabel("Wavelength (Å)")
plt.ylabel("Normalized Flux")

plt.legend()
plt.show()