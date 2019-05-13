from __future__ import print_function
import edibles.fit.avoigt as fit
import edibles.new_fit.fit2 as ft
import edibles.fit.make_grid as mg
import random
import numpy as np
import matplotlib.pyplot as plt


# ========================== Test 1 ===============================
#   creat dataset - without noise
#         lambda_c = 5890, R = 80,000, column density = 12.34,
#         Doppler parameter = 3.47 , continuum = 1.0, error = 0
# =================================================================
grid = mg.make_grid(5886, 5894, resolution=80000)
cont_const = 0.0
yy = 1 + fit.voigt(grid, lambda_peak=5890, b_eff=3.47, log_N=12.843, gamma=6.064e+07, osc_freq=0.631, resolving_power=80000)
flux = cont_const + yy

# interpolate into fixed grid
wave = np.linspace(5887, 5893, 1000)
flux = np.interp(wave, grid, flux)

# fit the line
obj = ft.fitSpectrum()
obj.afit(wave, flux, [5890], cheb_order=2, resolving_power=80000)


# compute scaled uncertainties
DOF = len(wave) - len(obj.fit_parm)
PCERROR = obj.fit_err * np.sqrt(obj.fit_norm/DOF)

print('')
print('  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
print('                              Test 1:')
print('  The resolution for this test is R=80,000 and the error on data')
print('  is zero further the continuum is the constant value 1.')
print('')
print('                    *** the fitting results ***')
print('               real parameter           fitted parameter')
print('  lambda_peak   :     5890                       {:7.2f} +- {:1.2f}'.format(obj.fit_parm[4], PCERROR[4]))
print('  b_eff         :     3.47                       {:7.2f} +- {:1.2f}'.format(obj.fit_parm[6], PCERROR[6]))
print('  log_N         :     12.843                     {:7.3f} +- {:1.2f}'.format(obj.fit_parm[7], PCERROR[7]))
print('')



# plotting
plt.figure()
plt.gcf().canvas.set_window_title('Test 1 - no error')
plt.plot(wave, flux, 'gray', marker='.')
plt.plot(wave, obj.yfit, 'red')
plt.xlim(5888.5,5891.5)
plt.xlabel('Wavelength ($\AA$)')
plt.ylabel('Flux')


plt.show()