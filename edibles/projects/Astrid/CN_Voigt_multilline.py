from edibles.utils.edibles_spectrum import EdiblesSpectrum
import numpy as np
import matplotlib.pyplot as plt
from voigt_profile import voigt_optical_depth


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

R0 = 3874.6059
R1 = 3873.9939
P1 = 3875.759


# =========
# 5. Continuo y profundidad base
# =========

continuum = np.median(flux_norm)


# =========
# 6. Parámetros físicos (de R0 ya afinados)
# =========

sigma = 0.014
gamma = 0.004

n = 1e23
v_rad=-14

f0 = 0.342 #Para R(0) = 3874.6059
f1= 0.0228 #Para R(1) = 3873.9939
f2= 0.0114 #Para P(1) = 3875.759


# =========
# 7. Perfiles Voigt individuales
# =========

profile_R0 = voigt_optical_depth(wave, lambda0=R0, b=sigma, gamma=gamma, N=n, v_rad=v_rad, f=f0)
profile_R1 = voigt_optical_depth(wave, lambda0=R1, b=sigma, gamma=gamma, N=n, v_rad=v_rad, f=f1)
profile_P1 = voigt_optical_depth(wave, lambda0=P1, b=sigma, gamma=gamma, N=n, v_rad=v_rad, f=f2)

print(profile_R0, profile_R1, profile_P1)

# =========
# 8. Profundidades relativas (IMPORTANTE)
# =========




# =========
# 9. Modelo multilínea (suma de τ)
# =========

tau_total = profile_R0 + profile_R1 + profile_P1


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


