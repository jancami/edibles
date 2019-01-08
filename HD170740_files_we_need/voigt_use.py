import sys
sys.path.append('/export/home/klay/github/')

import edibles.fit.avoigt as fit
import edibles.fit.avoigt_fit as ft
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np


# ========================== Test 7 ===============================
#   In this test, we added the following component to test 3
#   to check the multi-component Voigt fitting of edibles package
# =================================================================

# creat multi component dataset - with noise
wave = np.linspace(7664, 7667, 150)
# cont_const = 0.05 * (wave-7665) + 1000
cont_const = 0.05*(wave-7665)**(2)
yy1 = fit.voigt(wave, lambda_peak=7665.05, b_eff=3.47, log_N=12.34, gamma=6.064e+07, osc_freq=0.631)
yy2 = fit.voigt(wave, lambda_peak=7666.15, b_eff=4.13, log_N=12.34, gamma=6.064e+07, osc_freq=0.631)
yy3 = fit.voigt(wave, lambda_peak=7664.55, b_eff=2.5, log_N=12.0, gamma=6.064e+07, osc_freq=0.631)
flux = cont_const + yy1 + yy2 + yy3
noise = np.random.normal(0, 0.01, len(flux))
flux = flux + noise

obj = ft.fitSpectrum()
obj.afit(wave, flux, [[7665], [7666]], lines=['KI_7667'], cheb_order=1)



print ''
print '  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% '
print '                              Test 7'
print '  In this test, we added the following component to test 3'
print '  to check the multi-component Voigt fitting of edibles package'
print ''
print '                    *** the fitting results ***'
print '               real parameter           fitted parameter'
print '  continuum :     12.3 + 0.05 * x            {:7.2f}'.format(obj.fit_parm[0])
print '  lambda0_1 :     5890                       {:7.2f}'.format(obj.fit_parm[4])
print '  b_eff_1   :     3.47                       {:7.2f}'.format(obj.fit_parm[6])
print '  log_N_1   :     12.843                     {:7.3f}'.format(obj.fit_parm[7])
print ''
print '  lambda0_2 :     5890.12                    {:7.2f}'.format(obj.fit_parm[8])
print '  b_eff_2   :     4.13                       {:7.2f}'.format(obj.fit_parm[11])
print '  log_N_2   :     11.17                      {:7.3f}'.format(obj.fit_parm[12])
print ''

# plotting
plt.figure()
plt.gcf().canvas.set_window_title('Test 7 - two component')
plt.plot(wave, flux, 'gray', marker='.')
plt.plot(wave, obj.yfit, 'magenta')
plt.xlim(7663.5, 7667.5)
plt.xlabel('Wavelength ($\AA$)')
plt.ylabel('Flux')




# # ========================== Test 8 ===============================
# #   In this test, we added the following component to test 3
# #   to check the multi-component Voigt fitting of edibles package
# # =================================================================

# # creat multi transition dataset - with noise
# wave = np.arange(5887, 5898, 0.02)
# cont_const = 0.05 * (wave-5893) + 12.3
# yy1 = 1 + fit.voigt(wave, lambda0=5890, b_eff=3.47, log_N=12.843, gamma=6.064e+07, osc_freq=0.631)
# yy2 = 1 + fit.voigt(wave, lambda0=5896.1, b_eff=3.47, log_N=12.843, gamma=6.098e+07, osc_freq=0.318)
# flux = cont_const + yy1 + yy2
# noise = np.random.normal(0, 0.02, len(flux))
# flux = flux + noise

# obj = ft.fitSpectrum()
# obj.afit(wave, flux, [5890, 5896.1], lines=['NaI_5891', 'NaI_5897'], cheb_order=1)

# print ''
# print '  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% '
# print '                              Test 8'
# print '  In this test, we added the following component to test 3'
# print '  to check the multi-component NaI line'
# print ''
# print '                    *** the fitting results ***'
# print '               real parameter           fitted parameter'
# print '  continuum :     12.3 + 0.05 * x            {:7.2f}'.format(obj.fit_parm[0])
# print '  lambda0_1 :     5890                       {:7.2f}'.format(obj.fit_parm[3])
# print '  lambda0_2 :     5896.1                     {:7.2f}'.format(obj.fit_parm[5])
# print '  b_eff_1   :     3.47                       {:7.2f}'.format(obj.fit_parm[8])
# print '  log_N_1   :     12.843                     {:7.3f}'.format(obj.fit_parm[9])
# print ''

# # plotting
# plt.figure()
# plt.gcf().canvas.set_window_title('Test 8 - two component')
# plt.plot(wave, flux, 'gray', marker='.')
# plt.plot(wave, obj.yfit, 'magenta')
# plt.xlim(5888.5,5898)
# plt.xlabel('Wavelength ($\AA$)')
# plt.ylabel('Flux')


# # ========================== Test 9 ===============================
# #   This is CaII-K line from Welsh et al. 2010
# # =================================================================

# wave, flux = np.loadtxt('../data/caii.ascii', unpack=True)
# # fit the line
# obj = ft.fitSpectrum()
# obj.afit(wave, flux, [3933.8], lines=['CaII_3934'], cheb_order=1)


# print ''
# print '  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% '
# print '                              Test 9'
# print '  In this test, we added the following component to test 3'
# print '  to check the multi-component NaI line'
# print '  lambda0   :     3933.66                   {:7.2f}'.format(obj.fit_parm[4])
# print '  b_eff     :     2.5                       {:7.2f}'.format(obj.fit_parm[6])
# print '  log_N     :     11.67                     {:7.3f}'.format(obj.fit_parm[7])

# plt.figure()
# plt.gcf().canvas.set_window_title('Welsh et al. 2010 CaII-K')
# plt.plot(wave, flux, 'gray', marker='.')
# plt.plot(wave, obj.yfit, 'magenta')
# plt.xlabel('Wavelength ($\AA$)')
# plt.ylabel('Flux')
# plt.show()



# ========================== Use 10 ===============================
#   This is a real line somewhere
# =================================================================

# hdulist = fits.open('')
# hdu = hdulist[0]

# start = hdu.header['CRVAL1']
# step = hdu.header['CDELT1']
# flux = hdu.data
# length = len(flux)
# stop = start + step * length
# wave = np.linspace(start, stop, length)


file = 'txt/HD170740_20160613_reduced.txt'
wave, flux = np.loadtxt(file, unpack=True)

obj = ft.fitSpectrum()
obj.afit(wave, flux, [[7665.05], [7666.15]], lines=['KI_7667'], cheb_order=1)

# print ''
# print '  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% '
# print '                              20160613'
# print '  Known values                        Estimated values'
# print '--------------------------------------------------------------------'
# print '  lambda0   :     ?                  {:7.2f}'.format(obj.fit_parm[4])
# print '  b_eff     :     ?                  {:7.2f}'.format(obj.fit_parm[6])
# print '  log_N     :     ?                  {:7.3f}'.format(obj.fit_parm[7])


plt.figure()
plt.gcf().canvas.set_window_title(file)
plt.plot(wave, flux, 'gray', marker='.')
plt.plot(wave, obj.yfit, 'magenta')
plt.xlabel('Wavelength ($\AA$)')
plt.ylabel('Flux')



maxflux = np.max(flux)
flux = flux / maxflux

obj = ft.fitSpectrum()
obj.afit(wave, flux, [[7665.05], [7666.15]], lines=['KI_7667'], cheb_order=1)

plt.figure()
plt.gcf().canvas.set_window_title('normalized')
plt.plot(wave, flux, 'gray', marker='.')
plt.plot(wave, obj.yfit, 'magenta')
plt.xlabel('Wavelength ($\AA$)')
plt.ylabel('Flux')


cleaned = flux - obj.yfit


# obj = ft.fitSpectrum()
# obj.afit(wave, cleaned, [[7665.05], [7666.15]], lines=['KI_7667'], cheb_order=1)

plt.figure()
plt.plot(wave, cleaned)

# plt.gcf().canvas.set_window_title('cleaned')
# plt.plot(wave, flux, 'gray', marker='.')
# plt.plot(wave, obj.yfit, 'magenta')
# plt.xlabel('Wavelength ($\AA$)')
# plt.ylabel('Flux')


# file2 = 'txt/HD170740_20140915.txt'
# wave, flux = np.loadtxt(file2, unpack=True)

# obj = ft.fitSpectrum()
# obj.afit(wave, flux, [7664.1, 7665.2], lines=['KI_7667'], cheb_order=1)
# # obj.afit(wave, flux, [[7664.1], [7665.2]], lines=['KI_7667'], cheb_order=1)

# print ''
# print '  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% '
# print '                              20140915'
# print '  Known values                        Estimated values'
# print '--------------------------------------------------------------------'
# print '  lambda0   :     ?                  {:7.2f}'.format(obj.fit_parm[4])
# print '  b_eff     :     ?                  {:7.2f}'.format(obj.fit_parm[6])
# print '  log_N     :     ?                  {:7.3f}'.format(obj.fit_parm[7])


# plt.figure()
# plt.gcf().canvas.set_window_title(file2)
# plt.plot(wave, flux, 'gray', marker='.')
# plt.plot(wave, obj.yfit, 'magenta')
# plt.xlabel('Wavelength ($\AA$)')
# plt.ylabel('Flux')




# file3 = 'txt/HD170740_20140916.txt'
# wave, flux = np.loadtxt(file3, unpack=True)

# obj = ft.fitSpectrum()
# obj.afit(wave, flux, [7664.1, 7665.2], lines=['KI_7667'], cheb_order=1)
# # obj.afit(wave, flux, [[7664.1], [7665.2]], lines=['KI_7667'], cheb_order=1)

# print ''
# print '  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% '
# print '                              20140916'
# print '  Known values                        Estimated values'
# print '--------------------------------------------------------------------'
# print '  lambda0   :     ?                  {:7.2f}'.format(obj.fit_parm[4])
# print '  b_eff     :     ?                  {:7.2f}'.format(obj.fit_parm[6])
# print '  log_N     :     ?                  {:7.3f}'.format(obj.fit_parm[7])


# plt.figure()
# plt.gcf().canvas.set_window_title(file3)
# plt.plot(wave, flux, 'gray', marker='.')
# plt.plot(wave, obj.yfit, 'magenta')
# plt.xlabel('Wavelength ($\AA$)')
# plt.ylabel('Flux')


# file4 = 'txt/HD170740_20150626.txt'
# wave, flux = np.loadtxt(file4, unpack=True)

# obj = ft.fitSpectrum()
# obj.afit(wave, flux, [7664.1, 7665.2], lines=['KI_7667'], cheb_order=1)
# # obj.afit(wave, flux, [[7664.1], [7665.2]], lines=['KI_7667'], cheb_order=1)

# print ''
# print '  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% '
# print '                              20150626'
# print '  Known values                        Estimated values'
# print '--------------------------------------------------------------------'
# print '  lambda0   :     ?                  {:7.2f}'.format(obj.fit_parm[4])
# print '  b_eff     :     ?                  {:7.2f}'.format(obj.fit_parm[6])
# print '  log_N     :     ?                  {:7.3f}'.format(obj.fit_parm[7])


# plt.figure()
# plt.gcf().canvas.set_window_title(file4)
# plt.plot(wave, flux, 'gray', marker='.')
# plt.plot(wave, obj.yfit, 'magenta')
# plt.xlabel('Wavelength ($\AA$)')
# plt.ylabel('Flux')


# file5 = 'txt/HD170740_20170705.txt'
# wave, flux = np.loadtxt(file5, unpack=True)

# obj = ft.fitSpectrum()
# obj.afit(wave, flux, [7664.1, 7665.2], lines=['KI_7667'], cheb_order=1)
# # obj.afit(wave, flux, [[7664.1], [7665.2]], lines=['KI_7667'], cheb_order=1)

# print ''
# print '  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% '
# print '                              20170705'
# print '  Known values                        Estimated values'
# print '--------------------------------------------------------------------'
# print '  lambda0   :     ?                  {:7.2f}'.format(obj.fit_parm[4])
# print '  b_eff     :     ?                  {:7.2f}'.format(obj.fit_parm[6])
# print '  log_N     :     ?                  {:7.3f}'.format(obj.fit_parm[7])


# plt.figure()
# plt.gcf().canvas.set_window_title(file5)
# plt.plot(wave, flux, 'gray', marker='.')
# plt.plot(wave, obj.yfit, 'magenta')
# plt.xlabel('Wavelength ($\AA$)')
# plt.ylabel('Flux')


plt.show()
