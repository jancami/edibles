import numpy as np
import matplotlib.pyplot as plt
from edibles.utils import voigt_profile as vp
from edibles.utils.benchmark_voigt.parameter_modelling import file_reader
import astropy.constants as cst
import scipy.signal as ss

wavelengths, normflux = np.asarray(file_reader('omiper.m95.7698.txt'))

b= np.array([1.00])
N = np.array([1e12])
# need good initial guess for v_rad
# calc flux weighted wavelength and convert to v_rad
wave_weight = sum((1 - normflux) * wavelengths / sum(1-normflux))
print(wave_weight)

v_rad_guess = (wave_weight - 7698.974)/7698.974 * cst.c.to("km/s").value
print(v_rad_guess)
v_rad = np.array([v_rad_guess])

fit = vp.fit_multi_voigt_absorptionlines(wavegrid=wavelengths,
                                    ydata=normflux,
                                    restwave=7698.974,
                                    f=3.393e-1,
                                    gamma=3.8e7,
                                    b=b,
                                    N=N,
                                    v_rad=v_rad,
                                    v_resolution= 0.56,
                                    n_step=25)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(wavelengths, normflux, c = 'k', label = 'data')
#ax.plot(wavelengths, fit.best_fit, 'b', label = '1 component')
#ax.plot(wavelengths, normflux/fit.best_fit, c='r', label = 'flux/best fit')

peak_position = np.argwhere(normflux/fit.best_fit == np.max(normflux/fit.best_fit))
temp_array = normflux/fit.best_fit

local_min_1 = wavelengths[np.argmin(temp_array[:peak_position[0,0]])]
print('min 1', local_min_1)

local_min_2 = wavelengths[peak_position[0,0] + np.argmin(temp_array[peak_position[0,0]:])]
print('min 2', local_min_2)

wave_weight_1 = sum((1 - normflux[:peak_position[0,0]]) * wavelengths[:peak_position[0,0]] / sum(1-normflux[:peak_position[0,0]]))
print(wave_weight_1)
wave_weight_2 = sum((1 - normflux[peak_position[0,0]:]) * wavelengths[peak_position[0,0]:] / sum(1-normflux[peak_position[0,0]:]))
print(wave_weight_2)
v_rad[0] = (wave_weight_1 - 7698.974)/7698.974 * cst.c.to("km/s").value

v_rad = np.append(v_rad,(wave_weight_2 - 7698.974)/7698.974 * cst.c.to("km/s").value)


b =np.array([1,1])
N = np.array([1e12,1e12])
print('v_rad', v_rad)
print('b',b)
print('N',N)



fit_2 = vp.fit_multi_voigt_absorptionlines(wavegrid=wavelengths,
                                        ydata=normflux,
                                        restwave=7698.974,
                                        f=3.393e-1,
                                        gamma=3.8e7,
                                        b=b,
                                        N=N,
                                        v_rad=v_rad,
                                        v_resolution= 0.56,
                                        n_step=25)

temp_array = normflux/fit_2.best_fit
print(temp_array[np.argwhere((wavelengths > 7699.22) & (wavelengths < 7699.26))])
print('argrelextrma', ss.argrelextrema(temp_array[np.argwhere((wavelengths > 7699.23) & (wavelengths < 7699.4) & (temp_array<1))], np.less))
#print('find_peaks', ss.find_peaks(temp_array))
dydx = np.zeros(len(normflux)-1)
dydx_prime = np.zeros(len(normflux)-2)
dx = np.zeros(len(wavelengths)-1)
for i in range(len(normflux)-1):
    dy = temp_array[i] - temp_array[i-1]
    dx[i] = wavelengths[i] - wavelengths[i-1]
    dydx[i] = dy/dx[i]
for i in range(len(dydx)-1):
    dy_prime = dydx[i] - dydx[i-1]
    dx_prime = dx[i] - dx[i-1]
    dydx_prime[i] = dy_prime/dx_prime

ax.plot(wavelengths, fit_2.best_fit, c = 'orange', label = '2 components')
ax.plot(wavelengths, normflux/fit_2.best_fit, label = 'comparison 2')
#ax.plot(wavelengths[1:], (dydx/60) +1, label = 'dy/dx')
#ax.plot(wavelengths[2:], (dydx_prime/ 2e13) + 1, label = 'dy/dx prime')
plt.legend()

fit_2.params.pretty_print()
print('')
fit.params.pretty_print()

plt.show()

