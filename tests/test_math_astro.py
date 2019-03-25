from __future__ import print_function

from astropy import constants as cst
import edibles.fit.avoigt as fit
import edibles.fit.avoigt_fit as ft
import edibles.fit.make_grid as mg
import numpy as np
import matplotlib.pyplot as plt
from edibles.voigtMathematical import voigt_math
from edibles.astro_wrapper import voigt_astro



delta_v = 1000


R = cst.c.value / delta_v
grid = mg.make_grid(5886, 5894, resolution=R)
yy = 1 + fit.voigt(grid, lambda_peak=5890, b_eff=3.47, log_N=12.843, gamma=6.064e+07, osc_freq=0.631, resolving_power=R)
flux = yy
noise = np.random.normal(0, 0.02, len(flux))
flux = flux + noise


# interpolate into fixed grid
wave = np.linspace(5887, 5893, 1000)
flux = np.interp(wave, grid, flux)

obj = ft.fitSpectrum()
obj.afit(wave, flux, [5890], lines=['NaI_5891'], cheb_order=1, resolving_power=R)


# compute scaled uncertainties
DOF = len(wave) - len(obj.fit_parm)
PCERROR = obj.fit_err * np.sqrt(obj.fit_norm/DOF)

print('')
print('                    *** the fitting results ***')
print('               real parameter           fitted parameter')
print('  lambda_peak :     5890                       {:7.2f} +- {:1.2f}'.format(obj.fit_parm[4], PCERROR[4]))
print('  b_eff       :     3.47                       {:7.2f} +- {:1.2f}'.format(obj.fit_parm[6], PCERROR[6]))
print('  Gamma       :     6.064e7                    {:7.3f} +- {:1.2f}'.format(obj.fit_parm[5], PCERROR[5]))

print(obj.fit_parm)

print('')



# plotting
plt.figure()
plt.gcf().canvas.set_window_title('avoigt')
plt.plot(wave, flux, 'gray', marker='.')
plt.plot(wave, obj.yfit, 'red')
plt.xlim(5888.5, 5891.5)
plt.xlabel('Wavelength ($\AA$)')
plt.ylabel('Flux')
plt.legend()




cent = obj.fit_parm[4]
Gamma = obj.fit_parm[5]
b_eff = obj.fit_parm[6]

x2, y2 = voigt_astro(grid, cent, b_eff, Gamma)
plt.figure()
plt.plot(grid, y2, 'red', linestyle=':', label='astro')


# b_eff=3.47, Gamma=6.064e+07

alpha = 0.05760020918460115
gamma = 0.00048255778745462673

y3 = voigt_math(grid, cent, alpha, gamma)

plt.plot(grid, y3, 'blue', label='math')
plt.gcf().canvas.set_window_title('voigt_math and voigt_astro')
plt.legend()


difference = y3 - y2
plt.figure()
plt.plot(grid, difference)
plt.gcf().canvas.set_window_title('difference between voigt_math and voigt_astro')


plt.show()
