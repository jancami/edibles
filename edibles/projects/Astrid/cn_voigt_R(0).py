from edibles.utils.edibles_spectrum import EdiblesSpectrum
import numpy as np
import matplotlib.pyplot as plt
from voigt_profile import voigt_profile


# =========
# 1. Cargar espectro
# =========

filename = "orders/HD149757_2015-07-20T02:54:27_437nm_blue_O5.fits"
sp = EdiblesSpectrum(filename)

wave = sp.wave
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
# 4. Encontrar línea principal
# =========

idx_min = np.argmin(flux_norm)
center = wave[idx_min]


# =========
# 5. Parámetros físicos
# =========

continuum = np.median(flux_norm)

# profundidad correcta en tau
tau0 = -np.log(np.min(flux_norm) / continuum)

# ancho (ajustable ligeramente)
sigma = 0.016
gamma = 0.005


# =========
# 6. Perfil Voigt (CORRECTO)
# =========

profile = voigt_profile(wave - center, sigma, gamma)
profile = profile / np.max(profile)


# =========
# 7. Modelo físico FINAL
# =========

# factor para evitar que se vuelva plano
tau = 0.95 * tau0 * profile

fit_flux = continuum * np.exp(-tau)


# =========
# 8. Plot final
# =========

plt.figure(figsize=(10,6))

plt.plot(wave, flux_norm, label="Observed Spectrum", color="blue")
plt.plot(wave, fit_flux, label="Voigt Model (Final)", color="red")

plt.xlim(3873.0, 3876.0)

plt.xlabel("Wavelength (Å)")
plt.ylabel("Normalized Flux")

plt.legend()
plt.show()