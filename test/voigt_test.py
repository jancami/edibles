import edibles.fit.avoigt as fit
import edibles.fit.avoigt_fit as ft
import pickle
import numpy as np
import matplotlib.pyplot as plt


# ========================== Test 1 ===============================
#   creat dataset - without noise
#         lambda_c = 5890, R = 80,000, column density = 12.34,
#         Doppler parameter = 3.47 , continuum = 1.0, error = 0
# =================================================================

# create the dataset
wave = np.arange(5887, 5893, 0.07)
cont_const = 0.0
yy = 1 + fit.voigt(wave, lambda0=5890, b_eff=3.47, log_N=12.843, gamma=6.064e+07, osc_freq=0.631)
flux = cont_const + yy

# fit the line
obj = ft.fitSpectrum()
obj.afit(wave, flux, [5890], lines=['NaI_5891'], cheb_order=2)

print ''
print '  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% '
print '                              Test 1:'
print '  The resolution for this test is R=80,000 and the error on data'
print '  is zero further the continuum is the constant value 1.'
print ''
print '                    *** the fitting results ***'
print '               real parameter           fitted parameter'
print '  continuum :     1.0                        {:7.2f}'.format(obj.fit_parm[0])
print '  lambda0   :     5890                       {:7.2f}'.format(obj.fit_parm[4])
print '  b_eff     :     3.47                       {:7.2f}'.format(obj.fit_parm[6])
print '  log_N     :     12.843                     {:7.3f}'.format(obj.fit_parm[7])
print ''

# plotting
plt.figure()
plt.gcf().canvas.set_window_title('Test 1 - no error')
plt.plot(wave, flux, 'gray', marker='.')
plt.plot(wave, obj.yfit, 'magenta')
plt.xlim(5888.5,5891.5)
plt.xlabel('Wavelength ($\AA$)')
plt.ylabel('Flux')



# ========================== Test 2 ===============================
#   creat dataset - with noise
#         lambda_c = 5890, R = 80,000, column density = 12.34,
#         Doppler parameter = 3.47 , continuum = 1.0
# =================================================================

# generate dataset - with noise
wave = np.arange(5887, 5893, 0.07)
cont_const = 0.0
yy = 1 + fit.voigt(wave, lambda0=5890, b_eff=3.47, log_N=12.843, gamma=6.064e+07, osc_freq=0.631)
flux = cont_const + yy
np.random.seed(1)
noise = np.random.normal(0, 0.02, len(flux))
flux = flux + noise

obj = ft.fitSpectrum()
obj.afit(wave, flux, [5890], lines=['NaI_5891'], cheb_order=1)

print ''
print '  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% '
print '                              Test 2'
print '  This test exactly takes the test 1 parameters, however, we add'
print '  the error on the intensity with a Gaussian distribution of'
print '  width 0.05 of normalized intensity.'
print ''
print '                    *** the fitting results ***'
print '               real parameter           fitted parameter'
print '  continuum :     1.0                        {:7.2f}'.format(obj.fit_parm[0])
print '  lambda0   :     5890                       {:7.2f}'.format(obj.fit_parm[4])
print '  b_eff     :     3.47                       {:7.2f}'.format(obj.fit_parm[6])
print '  log_N     :     12.843                     {:7.3f}'.format(obj.fit_parm[7])
print ''

# plotting
plt.figure()
plt.gcf().canvas.set_window_title('Test 2 - with error')
plt.plot(wave, flux, 'gray', marker='.')
plt.plot(wave, obj.yfit, 'magenta')
plt.xlim(5888.5,5891.5)
plt.xlabel('Wavelength ($\AA$)')
plt.ylabel('Flux')


# ========================== Test 3 ===============================
#   creat dataset - with noise
#         lambda_c = 5890, R = 80,000, column density = 12.843,
#         Doppler parameter = 4.36 , continuum = 0.05 * (wave-5890) + 12.3
# =================================================================

wave = np.arange(5887, 5893, 0.07)
cont_const = 0.05 * (wave-5890) + 12.3
yy = 1 + fit.voigt(wave, lambda0=5890, b_eff=3.47, log_N=12.843, gamma=6.064e+07, osc_freq=0.631)
flux = cont_const + yy
np.random.seed(1)
noise = np.random.normal(0, 0.02, len(flux))
flux = flux + noise

obj = ft.fitSpectrum()
obj.afit(wave, flux, [5890], lines=['NaI_5891'], cheb_order=1)

print ''
print '  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% '
print '                              Test 3'
print '  This test exactly takes the test 2 parameters, however, we add'
print '  the continuum line of cnt = 0.005 * wave + 0.002'
print ''
print '                    *** the fitting results ***'
print '               real parameter           fitted parameter'
print '  continuum :     12.3 + 0.05 * x            {:7.2f}'.format(obj.fit_parm[0])
print '  lambda0   :     5890                       {:7.2f}'.format(obj.fit_parm[4])
print '  b_eff     :     3.47                       {:7.2f}'.format(obj.fit_parm[6])
print '  log_N     :     12.843                     {:7.3f}'.format(obj.fit_parm[7])
print ''

# plotting
plt.figure()
plt.gcf().canvas.set_window_title('Test 3 - with continuum')
plt.plot(wave, flux, 'gray', marker='.')
plt.plot(wave, obj.yfit, 'magenta')
plt.xlim(5888.5,5891.5)
plt.xlabel('Wavelength ($\AA$)')
plt.ylabel('Flux')



# ========================== Test 5 ===============================
#   This test is similar to test 3, with delat_lambnda = 0.02
# =================================================================

wave = np.arange(5887, 5893, 0.02)
cont_const = 0.05 * (wave-5890) + 12.3
yy = 1 + fit.voigt(wave, lambda0=5890, b_eff=3.47, log_N=12.843, gamma=6.064e+07, osc_freq=0.631)
flux = cont_const + yy
np.random.seed(1)
noise = np.random.normal(0, 0.02, len(flux))
flux = flux + noise

obj = ft.fitSpectrum()
obj.afit(wave, flux, [5890], lines=['NaI_5891'], cheb_order=1)

print ''
print '  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% '
print '                              Test 5'
print '  This test exactly takes the test 3 parameters, with fixed delta_lambda = 0.02'
print ''
print '                    *** the fitting results ***'
print '               real parameter           fitted parameter'
print '  continuum :     12.3 + 0.05 * x            {:7.2f}'.format(obj.fit_parm[0])
print '  lambda0   :     5890                       {:7.2f}'.format(obj.fit_parm[4])
print '  b_eff     :     3.47                       {:7.2f}'.format(obj.fit_parm[6])
print '  log_N     :     12.843                     {:7.3f}'.format(obj.fit_parm[7])
print ''

# plotting
plt.figure()
plt.gcf().canvas.set_window_title('Test 5 - delta_lambda=0.02')
plt.plot(wave, flux, 'gray', marker='.')
plt.plot(wave, obj.yfit, 'magenta')
plt.xlim(5888.5,5891.5)
plt.xlabel('Wavelength ($\AA$)')
plt.ylabel('Flux')





# ========================== Test 6 ===============================
#   This test is similar to test 3, with delat_lambnda = 0.02
# =================================================================

wave = 5887 + (5893-5887)*np.random.rand(1,200)
wave = wave[0]
wave.sort()
cont_const = 0.05 * (wave-5890) + 12.3
yy = 1 + fit.voigt(wave, lambda0=5890, b_eff=3.47, log_N=12.843, gamma=6.064e+07, osc_freq=0.631)
flux = cont_const + yy
np.random.seed(1)
noise = np.random.normal(0, 0.02, len(flux))
flux = flux + noise

obj = ft.fitSpectrum()
obj.afit(wave, flux, [5890], lines=['NaI_5891'], cheb_order=1)

print ''
print '  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% '
print '                              Test 6'
print '  This test is similar to test 3, with a nonlinear wavelength grid'
print ''
print '                    *** the fitting results ***'
print '               real parameter           fitted parameter'
print '  continuum :     12.3 + 0.05 * x            {:7.2f}'.format(obj.fit_parm[0])
print '  lambda0   :     5890                       {:7.2f}'.format(obj.fit_parm[4])
print '  b_eff     :     3.47                       {:7.2f}'.format(obj.fit_parm[6])
print '  log_N     :     12.843                     {:7.3f}'.format(obj.fit_parm[7])
print ''

# plotting
plt.figure()
plt.gcf().canvas.set_window_title('Test 6 - delta_lambda=0.02')
plt.plot(wave, flux, 'gray', marker='.')
plt.plot(wave, obj.yfit, 'magenta')
plt.xlim(5888.5,5891.5)
plt.xlabel('Wavelength ($\AA$)')
plt.ylabel('Flux')




# ========================== Test 7 ===============================
#   In this test, we added the following component to test 3
#   to check the multi-component Voigt fitting of edibles package
# =================================================================

# creat multi component dataset - with noise
wave = np.linspace(5887, 5893, 1000)
cont_const = 0.05 * (wave-5890) + 12.3
yy1 = 1 + fit.voigt(wave, lambda0=5890, b_eff=3.47, log_N=12.34, gamma=6.064e+07, osc_freq=0.631)
yy2 = 1 + fit.voigt(wave, lambda0=5890.12, b_eff=4.13, log_N=11.17, gamma=6.064e+07, osc_freq=0.631)
flux = cont_const + yy1 + yy2
noise = np.random.normal(0, 0.02, len(flux))
flux = flux + noise

obj = ft.fitSpectrum()
obj.afit(wave, flux, [[5890], [5890.1]], lines=['NaI_5891'], cheb_order=1)


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
plt.xlim(5888.5,5891.5)
plt.xlabel('Wavelength ($\AA$)')
plt.ylabel('Flux')



# ========================== Test 8 ===============================
#   In this test, we added the following component to test 3
#   to check the multi-component Voigt fitting of edibles package
# =================================================================

# creat multi transition dataset - with noise
wave = np.arange(5887, 5898, 0.02)
cont_const = 0.05 * (wave-5893) + 12.3
yy1 = 1 + fit.voigt(wave, lambda0=5890, b_eff=3.47, log_N=12.843, gamma=6.064e+07, osc_freq=0.631)
yy2 = 1 + fit.voigt(wave, lambda0=5896.1, b_eff=3.47, log_N=12.843, gamma=6.098e+07, osc_freq=0.318)
flux = cont_const + yy1 + yy2
noise = np.random.normal(0, 0.02, len(flux))
flux = flux + noise

obj = ft.fitSpectrum()
obj.afit(wave, flux, [5890, 5896.1], lines=['NaI_5891', 'NaI_5897'], cheb_order=1)

print ''
print '  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% '
print '                              Test 8'
print '  In this test, we added the following component to test 3'
print '  to check the multi-component NaI line'
print ''
print '                    *** the fitting results ***'
print '               real parameter           fitted parameter'
print '  continuum :     12.3 + 0.05 * x            {:7.2f}'.format(obj.fit_parm[0])
print '  lambda0_1 :     5890                       {:7.2f}'.format(obj.fit_parm[3])
print '  lambda0_2 :     5896.1                     {:7.2f}'.format(obj.fit_parm[5])
print '  b_eff_1   :     3.47                       {:7.2f}'.format(obj.fit_parm[8])
print '  log_N_1   :     12.843                     {:7.3f}'.format(obj.fit_parm[9])
print ''

# plotting
plt.figure()
plt.gcf().canvas.set_window_title('Test 8 - two component')
plt.plot(wave, flux, 'gray', marker='.')
plt.plot(wave, obj.yfit, 'magenta')
plt.xlim(5888.5,5898)
plt.xlabel('Wavelength ($\AA$)')
plt.ylabel('Flux')




# ========================== Test 9 ===============================
#   This is CaII-K line from Welsh et al. 2010
# =================================================================

wave, flux = np.loadtxt('../data/caii.ascii', unpack=True)
# fit the line
obj = ft.fitSpectrum()
obj.afit(wave, flux, [3933.8], lines=['CaII_3934'], cheb_order=1)


print ''
print '  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% '
print '                              Test 9'
print '  In this test, we added the following component to test 3'
print '  to check the multi-component NaI line'
print '  lambda0   :     3933.66                   {:7.2f}'.format(obj.fit_parm[4])
print '  b_eff     :     2.5                       {:7.2f}'.format(obj.fit_parm[6])
print '  log_N     :     11.67                     {:7.3f}'.format(obj.fit_parm[7])

plt.figure()
plt.gcf().canvas.set_window_title('Welsh et al. 2010 CaII-K')
plt.plot(wave, flux, 'gray', marker='.')
plt.plot(wave, obj.yfit, 'magenta')
plt.xlabel('Wavelength ($\AA$)')
plt.ylabel('Flux')
plt.show()
