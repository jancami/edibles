import edibles.fit.avoigt as fit
import edibles.fit.avoigt_fit as ft
import edibles.fit.make_grid as mg
import random
import numpy as np
import matplotlib.pyplot as plt
from edibles.fit.output import print_results

def broad(x, y, res):
    pxs = np.diff(x)[0] / x[0] * 299792.458
    fwhm_instrumental = res
    sigma_instrumental = fwhm_instrumental / 2.35482 / pxs
    LSF = gaussian(len(x)/2, sigma_instrumental)
    LSF = LSF / LSF.sum()
    y_profile = fftconvolve(y, LSF, 'same')

    return y_profile



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
obj.afit(wave, flux, [5890], lines=['NaI_5891'], cheb_order=2, resolving_power=80000)


# compute scaled uncertainties
DOF = len(wave) - len(obj.fit_parm)
PCERROR = obj.fit_err * np.sqrt(obj.fit_norm/DOF)

print ''
print '  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% '
print '                              Test 1:'
print '  The resolution for this test is R=80,000 and the error on data'
print '  is zero further the continuum is the constant value 1.'
print ''
print '                    *** the fitting results ***'
print '               real parameter           fitted parameter'
print '  lambda_peak   :     5890                       {:7.2f} +- {:1.2f}'.format(obj.fit_parm[4], PCERROR[4])
print '  b_eff         :     3.47                       {:7.2f} +- {:1.2f}'.format(obj.fit_parm[6], PCERROR[6])
print '  log_N         :     12.843                     {:7.3f} +- {:1.2f}'.format(obj.fit_parm[7], PCERROR[7])
print ''



# plotting
plt.figure()
plt.gcf().canvas.set_window_title('Test 1 - no error')
plt.plot(wave, flux, 'gray', marker='.')
plt.plot(wave, obj.yfit, 'red')
plt.xlim(5888.5,5891.5)
plt.xlabel('Wavelength ($\AA$)')
plt.ylabel('Flux')





# ========================== Test 2 ===============================
#   creat dataset - with noise
#         lambda_c = 5890, R = 80,000, column density = 12.34,
#         Doppler parameter = 3.47 , continuum = 1.0
# =================================================================

# generate dataset - with noise
grid = mg.make_grid(5886, 5894, resolution=80000)
cont_const = 0.0
yy = 1 + fit.voigt(grid, lambda_peak=5890, b_eff=3.47, log_N=12.843, gamma=6.064e+07, osc_freq=0.631, resolving_power=80000)
flux = cont_const + yy
np.random.seed(2)
noise = np.random.normal(0, 0.02, len(flux))
flux = flux + noise

# interpolate into fixed grid
wave = np.linspace(5887, 5893, 1000)
flux = np.interp(wave, grid, flux)

obj = ft.fitSpectrum()
obj.afit(wave, flux, [5890], lines=['NaI_5891'], cheb_order=1, resolving_power=80000)

# compute scaled uncertainties
DOF = len(wave) - len(obj.fit_parm)
PCERROR = obj.fit_err * np.sqrt(obj.fit_norm/DOF)

print ''
print '  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% '
print '                              Test 2'
print '  This test exactly takes the test 1 parameters, however, we add'
print '  the error on the intensity with a Gaussian distribution of'
print '  width 0.05 of normalized intensity.'
print ''
print '                    *** the fitting results ***'
print '               real parameter           fitted parameter'
print '  lambda_peak   :     5890                       {:7.2f} +- {:1.2f}'.format(obj.fit_parm[4], PCERROR[4])
print '  b_eff         :     3.47                       {:7.2f} +- {:1.2f}'.format(obj.fit_parm[6], PCERROR[6])
print '  log_N         :     12.843                     {:7.3f} +- {:1.2f}'.format(obj.fit_parm[7], PCERROR[7])
print ''




# plotting
plt.figure()
plt.gcf().canvas.set_window_title('Test 2 - with error')
plt.plot(wave, flux, 'gray', marker='.')
plt.plot(wave, obj.yfit, 'red')
plt.xlim(5888.5,5891.5)
plt.xlabel('Wavelength ($\AA$)')
plt.ylabel('Flux')


# ========================== Test 3 ===============================
#   creat dataset - with noise
#         lambda_c = 5890, R = 80,000, column density = 12.843,
#         Doppler parameter = 4.36 , continuum = 0.05 * (wave-5890) + 12.3
# =================================================================

grid = mg.make_grid(5886, 5894, resolution=80000)
cont_const = 0.05 * (grid-5890) + 12.3
yy = 1 + fit.voigt(grid, lambda_peak=5890, b_eff=3.47, log_N=12.843, gamma=6.064e+07, osc_freq=0.631, resolving_power=80000)
flux = cont_const + yy
noise = np.random.normal(0, 0.02, len(flux))
flux = flux + noise


# interpolate into fixed grid
wave = np.linspace(5887, 5893, 1000)
flux = np.interp(wave, grid, flux)

obj = ft.fitSpectrum()
obj.afit(wave, flux, [5890], lines=['NaI_5891'], cheb_order=1, resolving_power=80000)


# compute scaled uncertainties
DOF = len(wave) - len(obj.fit_parm)
PCERROR = obj.fit_err * np.sqrt(obj.fit_norm/DOF)

print ''
print '  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% '
print '                              Test 3'
print '  This test exactly takes the test 2 parameters, however, we add'
print '  the continuum line of cnt = 0.005 * wave + 0.002'
print ''
print '                    *** the fitting results ***'
print '               real parameter           fitted parameter'
print '  lambda_peak :     5890                       {:7.2f} +- {:1.2f}'.format(obj.fit_parm[4], PCERROR[4])
print '  b_eff       :     3.47                       {:7.2f} +- {:1.2f}'.format(obj.fit_parm[6], PCERROR[6])
print '  log_N       :     12.843                     {:7.3f} +- {:1.2f}'.format(obj.fit_parm[7], PCERROR[7])
print ''



# plotting
plt.figure()
plt.gcf().canvas.set_window_title('Test 3 - with continuum')
plt.plot(wave, flux, 'gray', marker='.')
plt.plot(wave, obj.yfit, 'red')
plt.xlim(5888.5,5891.5)
plt.xlabel('Wavelength ($\AA$)')
plt.ylabel('Flux')







# ========================== Test 4 ===============================
#   creat dataset - with noise
#         lambda_c = 5890, R = 10,000, column density = 12.843,
#         Doppler parameter = 4.36 , continuum = 0.05 * (wave-5890) + 12.3
# =================================================================

grid = mg.make_grid(5886, 5894, resolution=10000)
cont_const = 0.05 * (grid-5890) + 12.3
yy = 1 + fit.voigt(grid, lambda_peak=5890, b_eff=2.47, log_N=12.843, gamma=6.064e+07, osc_freq=0.631, resolving_power=10000)
flux = cont_const + yy
noise = np.random.normal(0, 0.02, len(flux))
flux = flux + noise


# interpolate into fixed grid
wave = np.linspace(5887, 5893, 1000)
flux = np.interp(wave, grid, flux)

obj = ft.fitSpectrum()
obj.afit(wave, flux, [5890], lines=['NaI_5891'], cheb_order=1, resolving_power=10000)


# compute scaled uncertainties
DOF = len(wave) - len(obj.fit_parm)
PCERROR = obj.fit_err * np.sqrt(obj.fit_norm/DOF)


print ''
print '  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% '
print '                              Test 4'
print '  Similar to test 3 with R=10,000'
print ''
print '                    *** the fitting results ***'
print '               real parameter           fitted parameter'
print '  lambda_peak   :     5890                       {:7.2f} +- {:1.2f}'.format(obj.fit_parm[4], PCERROR[4])
print '  b_eff         :     3.47                       {:7.2f} +- {:1.2f}'.format(obj.fit_parm[6], PCERROR[6])
print '  log_N         :     12.843                     {:7.3f} +- {:1.2f}'.format(obj.fit_parm[7], PCERROR[7])
print ''

# plotting
plt.figure()
plt.gcf().canvas.set_window_title('Test 4 - R=10,000')
plt.plot(wave, flux, 'gray', marker='.')
plt.plot(wave, obj.yfit, 'red')
plt.xlim(5888.5,5891.5)
plt.xlabel('Wavelength ($\AA$)')
plt.ylabel('Flux')





# ========================== Test 5 ===============================
#   This test is similar to test 3, with delat_lambnda = 0.02
# =================================================================

wave = np.arange(5887, 5893, 0.02)
cont_const = 0.05 * (wave-5890) + 12.3
yy = 1 + fit.voigt(wave, lambda_peak=5890, b_eff=3.47, log_N=12.843, gamma=6.064e+07, osc_freq=0.631, resolving_power=80000)
flux = cont_const + yy
noise = np.random.normal(0, 0.02, len(flux))
flux = flux + noise

obj = ft.fitSpectrum()
obj.afit(wave, flux, [5890], lines=['NaI_5891'], cheb_order=1, resolving_power=80000)


# compute scaled uncertainties
DOF = len(wave) - len(obj.fit_parm)
PCERROR = obj.fit_err * np.sqrt(obj.fit_norm/DOF)


print ''
print '  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% '
print '                              Test 5'
print '  This test exactly takes the test 3 parameters, with fixed delta_lambda = 0.02'
print ''
print '                    *** the fitting results ***'
print '               real parameter           fitted parameter'
print '  lambda_peak   :     5890                       {:7.2f} +- {:1.2f}'.format(obj.fit_parm[4], PCERROR[4])
print '  b_eff         :     3.47                       {:7.2f} +- {:1.2f}'.format(obj.fit_parm[6], PCERROR[6])
print '  log_N         :     12.843                     {:7.3f} +- {:1.2f}'.format(obj.fit_parm[7], PCERROR[7])
print ''

# plotting
plt.figure()
plt.gcf().canvas.set_window_title('Test 5 - delta_lambda=0.02')
plt.plot(wave, flux, 'gray', marker='.')
plt.plot(wave, obj.yfit, 'red')
plt.xlim(5888.5,5891.5)
plt.xlabel('Wavelength ($\AA$)')
plt.ylabel('Flux')





# ========================== Test 6 ===============================
#   This test is similar to test 3, with non-linear grid
# =================================================================
grid = mg.make_grid(5886, 5894, resolution=80000)
cont_const = 0.05 * (grid-5890) + 12.3
yy = 1 + fit.voigt(grid, lambda_peak=5890, b_eff=3.47, log_N=12.843, gamma=6.064e+07, osc_freq=0.631, resolving_power=80000)
flux = cont_const + yy
noise = np.random.normal(0, 0.02, len(flux))
flux = flux + noise


# interpolate into fixed grid
wave = np.linspace(5887, 5893, 1000).tolist()
index = random.randrange(len(wave))
for it in range(120):
    index = random.randrange(len(wave))
    del wave[index]
flux = np.interp(wave, grid, flux)


obj = ft.fitSpectrum()
obj.afit(wave, flux, [5890], lines=['NaI_5891'], cheb_order=1, resolving_power=80000)

# compute scaled uncertainties
DOF = len(wave) - len(obj.fit_parm)
PCERROR = obj.fit_err * np.sqrt(obj.fit_norm/DOF)


print ''
print '  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% '
print '                              Test 6'
print '  This test is similar to test 3, with a nonlinear wavelength grid'
print ''
print '                    *** the fitting results ***'
print '               real parameter           fitted parameter'
print '  lambda_peak   :     5890                       {:7.2f} +- {:1.2f}'.format(obj.fit_parm[4], PCERROR[4])
print '  b_eff         :     3.47                       {:7.2f} +- {:1.2f}'.format(obj.fit_parm[6], PCERROR[6])
print '  log_N         :     12.843                     {:7.3f} +- {:1.2f}'.format(obj.fit_parm[7], PCERROR[7])
print ''

# plotting
plt.figure()
plt.gcf().canvas.set_window_title('Test 6 - delta_lambda=0.02')
plt.plot(wave, flux, 'gray', marker='.')
plt.plot(wave, obj.yfit, 'red')
plt.xlim(5888.5,5891.5)
plt.xlabel('Wavelength ($\AA$)')
plt.ylabel('Flux')




# ========================== Test 7 ===============================
#   In this test, we added the following component to test 3
#   to check the multi-component Voigt fitting of edibles package
# =================================================================

# create multi-component dataset + with noise
grid = mg.make_grid(5886, 5894, resolution=80000)
cont_const = 0.05 * (grid-5890) + 12.3
yy1 = 1 + fit.voigt(grid, lambda_peak=5890, b_eff=3.47, log_N=12.34, gamma=6.064e+07, osc_freq=0.631, resolving_power=80000)
yy2 = 1 + fit.voigt(grid, lambda_peak=5890.12, b_eff=4.13, log_N=11.17, gamma=6.064e+07, osc_freq=0.631, resolving_power=80000)
flux = cont_const + yy1 + yy2
noise = np.random.normal(0, 0.02, len(flux))
flux = flux + noise


# interpolate into fixed grid
wave = np.linspace(5887, 5893, 1000)
flux = np.interp(wave, grid, flux)


obj = ft.fitSpectrum()
obj.afit(wave, flux, [[5890], [5890.1]], lines=['NaI_5891'], cheb_order=1, resolving_power=80000)


# compute scaled uncertainties
DOF = len(wave) - len(obj.fit_parm)
PCERROR = obj.fit_err * np.sqrt(obj.fit_norm/DOF)

print ''
print '  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% '
print '                              Test 7'
print '  In this test, we added the following component to test 3'
print '  to check the multi-component Voigt fitting of edibles package'
print ''
print '                    *** the fitting results ***'
print '               real parameter           fitted parameter'
print '  lambda_peak_1 :     5890                       {:7.2f} +- {:1.2f}'.format(obj.fit_parm[3], PCERROR[3])
print '  b_eff_1       :     3.47                       {:7.2f} +- {:1.2f}'.format(obj.fit_parm[6], PCERROR[6])
print '  log_N_1       :     12.843                     {:7.3f} +- {:1.2f}'.format(obj.fit_parm[7], PCERROR[7])
print ''
print '  lambda_peak_2 :     5890.12                    {:7.2f} +- {:1.2f}'.format(obj.fit_parm[8], PCERROR[8])
print '  b_eff_2       :     4.13                       {:7.2f} +- {:1.2f}'.format(obj.fit_parm[11], PCERROR[11])
print '  log_N_2       :     11.17                      {:7.3f} +- {:1.2f}'.format(obj.fit_parm[12], PCERROR[12])
print ''



# plotting
plt.figure()
plt.gcf().canvas.set_window_title('Test 7 - two component')
plt.plot(wave, flux, 'gray', marker='.')
plt.plot(wave, obj.yfit, 'red')
plt.xlim(5888.5,5891.5)
plt.xlabel('Wavelength ($\AA$)')
plt.ylabel('Flux')



# ========================== Test 8 ===============================
#   In this test, we added the following component to test 3
#   to check the multi-component Voigt fitting of edibles package
# =================================================================

# create multi transition dataset - with noise
grid = mg.make_grid(5886, 5899, resolution=80000)
cont_const = 0.05 * (grid-5893) + 12.3
yy1 = 1 + fit.voigt(grid, lambda_peak=5890, b_eff=3.47, log_N=12.843, gamma=6.064e+07, osc_freq=0.631, resolving_power=80000)
yy2 = 1 + fit.voigt(grid, lambda_peak=5896.1, b_eff=3.47, log_N=12.843, gamma=6.098e+07, osc_freq=0.318, resolving_power=80000)
flux = cont_const + yy1 + yy2
noise = np.random.normal(0, 0.02, len(flux))
flux = flux + noise


# interpolate into fixed grid
wave = np.linspace(5887, 5898, 1000)
flux = np.interp(wave, grid, flux)



obj = ft.fitSpectrum()
obj.afit(wave, flux, [5890, 5896.1], lines=['NaI_5891', 'NaI_5897'], cheb_order=1, resolving_power=80000)


# compute scaled uncertainties
DOF = len(wave) - len(obj.fit_parm)
PCERROR = obj.fit_err * np.sqrt(obj.fit_norm/DOF)


print ''
print '  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% '
print '                              Test 8'
print '  In this test, we added the following component to test 3'
print '  to check the multi-component NaI line'
print ''
print '                    *** the fitting results ***'
print '               real parameter           fitted parameter'
print '  lambda_peak_1 :     5890                       {:7.2f} +- {:1.2f}'.format(obj.fit_parm[3], PCERROR[3])
print '  lambda_peak_2 :     5896.1                     {:7.2f} +- {:1.2f}'.format(obj.fit_parm[5], PCERROR[5])
print '  b_eff         :     3.47                       {:7.2f} +- {:1.2f}'.format(obj.fit_parm[8], PCERROR[8])
print '  log_N         :     12.843                     {:7.3f} +- {:1.2f}'.format(obj.fit_parm[9], PCERROR[9])
print ''

# plotting
plt.figure()
plt.gcf().canvas.set_window_title('Test 8 - two component')
plt.plot(wave, flux, 'gray', marker='.')
plt.plot(wave, obj.yfit, 'red')
plt.xlim(5888.5,5898)
plt.xlabel('Wavelength ($\AA$)')
plt.ylabel('Flux')






# ========================== Test 9 ===============================
#   In this test, we added the following component to test 3
#   to check the multi-component Voigt fitting of edibles package
# =================================================================

# create multi transition dataset - with noise
grid = mg.make_grid(5886, 5899, resolution=80000)
cont_const = 0.05 * (grid-5893) + 12.3
yy1 = 1 + fit.voigt(grid, lambda_peak=5890, b_eff=3.47, log_N=12.843, gamma=6.064e+07, osc_freq=0.631, resolving_power=80000)
yy2 = 1 + fit.voigt(grid, lambda_peak=5896, b_eff=3.47, log_N=12.843, gamma=6.098e+07, osc_freq=0.318, resolving_power=80000)
yy3 = 1 + fit.voigt(grid, lambda_peak=5890.3, b_eff=4.13, log_N=11.17, gamma=6.064e+07, osc_freq=0.631, resolving_power=80000)
yy4 = 1 + fit.voigt(grid, lambda_peak=5896.3, b_eff=4.13, log_N=11.17, gamma=6.098e+07, osc_freq=0.318, resolving_power=80000)
flux = cont_const + yy1 + yy2 + yy3 + yy4
noise = np.random.normal(0, 0.02, len(flux))
flux = flux + noise


# interpolate into fixed grid
wave = np.linspace(5887, 5898, 1000)
flux = np.interp(wave, grid, flux)


obj = ft.fitSpectrum()
obj.afit(wave, flux, [ [5890, 5896], [5890.3, 5896.3] ], lines=['NaI_5891', 'NaI_5897'], cheb_order=1, resolving_power=80000)


# compute scaled uncertainties
DOF = len(wave) - len(obj.fit_parm)
PCERROR = obj.fit_err * np.sqrt(obj.fit_norm/DOF)


print ''
print '  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% '
print '                              Test 9'
print '  In this test, we added the following component to test 3'
print '  to check the multi-component NaI line'
print ''
print '                    *** the fitting results ***'
print '               real parameter           fitted parameter'
print '  lambda1_1 :     5890                       {:7.2f} +- {:1.2f}'.format(obj.fit_parm[3], PCERROR[3])
print '  lambda1_2 :     5896                       {:7.2f} +- {:1.2f}'.format(obj.fit_parm[5], PCERROR[5])
print '  b_eff_1   :     3.47                       {:7.2f} +- {:1.2f}'.format(obj.fit_parm[8], PCERROR[8])
print '  log_N_1   :     12.843                     {:7.3f} +- {:1.2f}'.format(obj.fit_parm[9], PCERROR[9])
print ''
print '  lambda2_1 :     5890.3                     {:7.2f} +- {:1.2f}'.format(obj.fit_parm[10], PCERROR[10])
print '  lambda2_2 :     5896.3                     {:7.2f} +- {:1.2f}'.format(obj.fit_parm[12], PCERROR[12])
print '  b_eff_2   :     4.13                       {:7.2f} +- {:1.2f}'.format(obj.fit_parm[15], PCERROR[15])
print '  log_N_2   :     11.17                      {:7.3f} +- {:1.2f}'.format(obj.fit_parm[16], PCERROR[16])
print ''


# plotting
plt.figure()
plt.gcf().canvas.set_window_title('Test 9 - two components')
plt.plot(wave, flux, 'gray', marker='.')
plt.plot(wave, obj.yfit, 'red')
plt.xlim(5888.5,5898)
plt.xlabel('Wavelength ($\AA$)')
plt.ylabel('Flux')




# ========================== Test 10 ===============================
#   This is CaII-K line from Welsh et al. 2010
# =================================================================

wave, flux = np.loadtxt('../data/caii.ascii', unpack=True)
# fit the line
obj = ft.fitSpectrum()
obj.afit(wave, flux, [3933.8], lines=['CaII_3934'], cheb_order=1, resolving_power=80000)

# compute scaled uncertainties
DOF = len(wave) - len(obj.fit_parm)
PCERROR = obj.fit_err * np.sqrt(obj.fit_norm/DOF)


print ''
print '  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% '
print '                              Test 10'
print '  In this test, we added the following component to test 3'
print '  to check the multi-component NaI line'
print '  lambda_peak   :     3933.66                   {:7.2f} +- {:1.2f}'.format(obj.fit_parm[4], PCERROR[4])
print '  b_eff         :     2.5                       {:7.2f} +- {:1.2f}'.format(obj.fit_parm[6], PCERROR[6])
print '  log_N         :     11.67                     {:7.3f} +- {:1.2f}'.format(obj.fit_parm[7], PCERROR[7])

plt.figure()
plt.gcf().canvas.set_window_title('Welsh et al. 2010 CaII-K')
plt.plot(wave, flux, 'gray', marker='.')
plt.plot(wave, obj.yfit, 'red')
plt.xlabel('Wavelength ($\AA$)')
plt.ylabel('Flux')
plt.show()
