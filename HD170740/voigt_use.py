from __future__ import print_function
import sys
import os

path = os.getcwd()
print(sys.path)
os.chdir('..')
os.chdir('..')
sys.path.append(os.getcwd())
os.chdir(path)
print(sys.path)
from edibles import fit, more
import edibles.fit.avoigt as fit
import edibles.fit.avoigt_fit as ft
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits


sys.dont_write_bytecode = True

# ===================== Test Data ==========================
#   Three components, made to look like the real data from
#   the file 'txt/HD170740_20160613_reduced.txt'
# ==========================================================

# creat multi component dataset - with noise
wave = np.linspace(7664, 7667, 150)
cont_const = 0.05 * (wave-7665)+1
# cont_const = (0.5*(wave-7665)**(2)) + 10
# cont_const = 10.0
yy1 = 1*fit.voigt(wave, lambda_peak=7665.05, b_eff=2.85, log_N=11.983, gamma=6.064e+07, osc_freq=0.631)
yy2 = 1*fit.voigt(wave, lambda_peak=7666.15, b_eff=3.12, log_N=12.047, gamma=6.064e+07, osc_freq=0.631)
yy3 = 1*fit.voigt(wave, lambda_peak=7664.55, b_eff=2.5, log_N=11.7, gamma=6.064e+07, osc_freq=0.631)
flux = cont_const + yy1 + yy2 + yy3
noise = np.random.normal(0, 0.01, len(flux))
flux = flux + noise

obj = ft.fitSpectrum()
obj.afit(wave, flux, [[7665], [7666]], lines=['KI_7667'], cheb_order=1)
# obj.afit(wave, flux, [7665], lines=['KI_7667'], cheb_order=1)

DOF = len(wave) - len(obj.fit_parm)
PCERROR = obj.fit_err * np.sqrt(obj.fit_norm/DOF)


# --------------------------
# Fitting a line to the data
# --------------------------
mean = np.mean(flux)
weights = flux / mean
# weights = 1 - abs(flux - mean)

# # quad
# fit = np.polyfit(wave, flux, 2)
# better_fit = fit[0]*wave**2 + fit[1]*wave + fit[2]

# linear
fit = np.polyfit(wave, flux, 1, w=weights)
better_fit = fit[0]*wave + fit[1]

# --------------------------
# Cleaning the spectrum
# --------------------------
new_flux = flux - obj.yfit


print('')
print('  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
print('                              Test Data')
print('')
print('                    *** the fitting results ***')
print('               real parameter           fitted parameter')
print('  continuum :     line                      {:7.2f} +- {:1.2f}'.format(obj.fit_parm[0], PCERROR[0]))
print('')
print('  lambda0_1 :     7665.05                   {:7.2f} +- {:1.2f}'.format(obj.fit_parm[4], PCERROR[4]))
print('  b_eff_1   :     2.85                      {:7.2f} +- {:1.2f}'.format(obj.fit_parm[6], PCERROR[6]))
print('  log_N_1   :     11.983                    {:7.3f} +- {:1.2f}'.format(obj.fit_parm[7], PCERROR[7]))
print('')
print('  lambda0_2 :     7666.15                   {:7.2f} +- {:1.2f}'.format(obj.fit_parm[8], PCERROR[8]))
print('  b_eff_2   :     3.12                      {:7.2f} +- {:1.2f}'.format(obj.fit_parm[11], PCERROR[11]))
print('  log_N_2   :     12.047                    {:7.3f} +- {:1.2f}'.format(obj.fit_parm[12], PCERROR[12]))
print('')

# plotting

fig, ax1 = plt.subplots()
plt.gcf().canvas.set_window_title('test data')
ax1.plot(wave, flux, 'gray', marker='.', label='data')
ax1.plot(wave, obj.yfit, 'magenta', label='fit')
cont_fit = np.ones_like(wave) * obj.fit_parm[0]
ax1.plot(wave, cont_fit, 'b--', label='average')
ax1.plot(wave, better_fit, 'r--', label='poly')


# ax2 = ax1.twinx()

ax1.plot(wave, new_flux, 'green', marker='.', label='cleaned')
plt.xlabel('Wavelength ($\AA$)')
ax1.set_ylabel('Flux')

plt.legend()
# --------------------------
# data / red line
# --------------------------














# ===================== Raw Data ==========================
#             HD170740 from June 13, 2016
# ==========================================================
file = 'txt/HD170740_20160613_reduced.txt'
wave, flux = np.loadtxt(file, unpack=True)

obj = ft.fitSpectrum()
obj.afit(wave, flux, [[7665.05], [7666.15], [7664.5]], lines=['KI_7667'], cheb_order=1)

# print(''
# print('  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% '
# print('                              20160613'
# print('  Known values                        Estimated values'
# print('--------------------------------------------------------------------'
# print('  lambda0   :     ?                  {:7.2f}'.format(obj.fit_parm[4])
# print('  b_eff     :     ?                  {:7.2f}'.format(obj.fit_parm[6])
# print('  log_N     :     ?                  {:7.3f}'.format(obj.fit_parm[7])

DOF = len(wave) - len(obj.fit_parm)
PCERROR = obj.fit_err * np.sqrt(obj.fit_norm/DOF)


# --------------------------
# Fitting a line to the data
# --------------------------
mean = np.mean(flux)
weights = (flux / mean)

# weights = 1 - abs(flux - mean)

# # quad
# fit = np.polyfit(wave, flux, 2)
# better_fit = fit[0]*wave**2 + fit[1]*wave + fit[2]

# linear
fit = np.polyfit(wave, flux, 1, w=weights)
better_fit = fit[0]*wave + fit[1]

print('')
print('  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
print('                              Raw data')
print('')
print('                    *** the fitting results ***')
print('                  fitted parameter')
print('  continuum :     {:7.2f} +- {:1.2f}'.format(obj.fit_parm[0], PCERROR[0]))
print('  lambda0_1 :     {:7.2f} +- {:1.2f}'.format(obj.fit_parm[4], PCERROR[4]))
print('  b_eff_1   :     {:7.2f} +- {:1.2f}'.format(obj.fit_parm[6], PCERROR[6]))
print('  log_N_1   :     {:7.3f} +- {:1.2f}'.format(obj.fit_parm[7], PCERROR[7]))
print('')
print('  lambda0_2 :     {:7.2f} +- {:1.2f}'.format(obj.fit_parm[8], PCERROR[8]))
print('  b_eff_2   :     {:7.2f} +- {:1.2f}'.format(obj.fit_parm[11], PCERROR[11]))
print('  log_N_2   :     {:7.3f} +- {:1.2f}'.format(obj.fit_parm[12], PCERROR[12]))
print('')


fig, ax1 = plt.subplots()
plt.gcf().canvas.set_window_title(file)
ax1.plot(wave, flux, 'gray', marker='.', label='HD170740')
ax1.plot(wave, obj.yfit, 'magenta', label='fit')
cont_fit = np.ones_like(wave) * obj.fit_parm[0]
ax1.plot(wave, cont_fit, 'blue', label='average')
ax1.plot(wave, better_fit, 'red', label='poly')

# ax2 = ax1.twinx()
# ax2.plot(wave, weights, 'green')
plt.xlabel('Wavelength ($\AA$)')
ax1.set_ylabel('Flux')
plt.legend()


# maxflux = np.max(flux)
# flux = flux / maxflux

# obj = ft.fitSpectrum()
# obj.afit(wave, flux, [[7665.05], [7666.15]], lines=['KI_7667'], cheb_order=1)

# print('')
# print('  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
# print('                              Normalized data')
# print('  In this test, we added the following component to test 3')
# print('  to check the multi-component Voigt fitting of edibles package')
# print('')
# print('                    *** the fitting results ***')
# print('                  fitted parameter')
# print('  continuum :     {:7.2f} +- {:1.2f}'.format(obj.fit_parm[0], PCERROR[0]))
# print('  lambda0_1 :     {:7.2f} +- {:1.2f}'.format(obj.fit_parm[4], PCERROR[4]))
# print('  b_eff_1   :     {:7.2f} +- {:1.2f}'.format(obj.fit_parm[6], PCERROR[6]))
# print('  log_N_1   :     {:7.3f} +- {:1.2f}'.format(obj.fit_parm[7], PCERROR[7]))
# print('')
# print('  lambda0_2 :     {:7.2f} +- {:1.2f}'.format(obj.fit_parm[8], PCERROR[8]))
# print('  b_eff_2   :     {:7.2f} +- {:1.2f}'.format(obj.fit_parm[11], PCERROR[11]))
# print('  log_N_2   :     {:7.3f} +- {:1.2f}'.format(obj.fit_parm[12], PCERROR[12]))
# print('')


# plt.figure()
# plt.gcf().canvas.set_window_title('normalized')
# plt.plot(wave, flux, 'gray', marker='.')
# plt.plot(wave, obj.yfit, 'magenta')
# cont_fit = np.ones_like(wave) * obj.fit_parm[0]
# plt.plot(wave, cont_fit, 'blue')
# plt.xlabel('Wavelength ($\AA$)')
# plt.ylabel('Flux')


plt.show()
