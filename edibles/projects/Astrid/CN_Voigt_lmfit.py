from edibles.utils.edibles_spectrum import EdiblesSpectrum
import numpy as np
import matplotlib.pyplot as plt
from edibles.utils.voigt_profile import voigt_optical_depth
from lmfit import Model


# =========
# 1. Cargar espectro
# =========


filename = "orders/HD149757_2015-07-20T02:54:27_437nm_blue_O5.fits"
sp = EdiblesSpectrum(filename)


wave = sp.bary_wave
flux = sp.flux



# =========
# 3. Seleccionar región CN
# =========


mask = (wave > 3873.0) & (wave < 3876.0)
wave = wave[mask]
flux = flux[mask]

# =========
# 2. Normalizar
# =========


flux_norm = flux / np.median(flux)


# =========
# 4. Líneas CN (valores de las tablas)
# =========


R0 = 3874.6059
R1 = 3873.9939
P1 = 3875.759


# =========
# 5. Parámetros físicos (de R0 ya afinados)
# =========


f0 = 0.342  # Para R(0) = 3874.6059
f1 = 0.0228 # Para R(1) = 3873.9939
f2 = 0.0114 # Para P(1) = 3875.759


# =========
# 6. Definición del Modelo Físico
# =========


def cn_model(wave, sigma, N0, N1, v_rad):
  gamma = 1e6
  profile_R0 = f0 * voigt_optical_depth(wave, lambda0=R0, b=sigma, gamma=gamma, N=N0, v_rad=v_rad)
  profile_R1 = f1 * voigt_optical_depth(wave, lambda0=R1, b=sigma, gamma=gamma, N=N1, v_rad=v_rad)
  profile_P1 = f2 * voigt_optical_depth(wave, lambda0=P1, b=sigma, gamma=gamma, N=N1, v_rad=v_rad)

  print("profile_R0:", profile_R0)

  tau_total = profile_R0 + profile_R1 + profile_P1

  print("tau min/max:", np.min(tau_total), np.max(tau_total))

  return np.exp(-tau_total)


# =========
# 7. Ajuste con lmfit
# =========


# Creamos el objeto Model vinculándolo a tu función cn_model

model = Model(cn_model, independent_vars=['wave'])


params = model.make_params()

params['N0'].value = 1e30
params['N0'].min = 1e12
params['N0'].max = 1e30

params['N1'].value = 1e13
params['N1'].min = 1e12
params['N1'].max = 1e15

params['sigma'].value = 1
params['sigma'].min = 0.1
params['sigma'].max = 10

params[('v_rad')].value = 0
params['v_rad'].min = -100
params['v_rad'].max = +100


guess=model.eval(params=params,wave=wave)
print("guess:", guess)


# Ejecución del ajuste espectral
result = model.fit(flux_norm, params, wave=wave)
print(result.fit_report())


# =========
# 8. Plot final
# =========


plt.figure(figsize=(10,6))


plt.plot(wave, flux_norm, label="Observed", color="blue")
plt.plot(wave, result.best_fit, label="Voigt fit (lmfit)", color="red")

plt.plot(wave, guess, label="guess", color="violet")

plt.legend()
plt.xlim(3873, 3876)


plt.xlabel("Wavelength (Å)")
plt.ylabel("Normalized Flux")


plt.show()
