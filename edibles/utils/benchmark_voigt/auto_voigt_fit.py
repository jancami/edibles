import numpy as np
import matplotlib.pyplot as plt
from edibles import PYTHONDIR
from pathlib import Path
import pandas as pd
from edibles.utils import voigt_profile as vp
import os.path
from edibles.utils.benchmark_voigt.parameter_modelling import file_reader
import astropy.constants as cst

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
ax.plot(wavelengths, fit.best_fit, 'b', label = '1 component')
ax.plot(wavelengths, normflux/fit.best_fit, c='r', label = 'flux/best fit')

peak_position = np.argwhere(normflux/fit.best_fit == np.max(normflux/fit.best_fit))
temp_array = normflux/fit.best_fit

local_min_1 = wavelengths[np.argmin(temp_array[:peak_position[0,0]])]
print('min 1', local_min_1)

local_min_2 = wavelengths[peak_position[0,0] + np.argmin(temp_array[peak_position[0,0]:])]
print('min 2', local_min_2)

v_rad[0] = (wave_weight - local_min_1)/local_min_1 * cst.c.to("km/s").value

v_rad = np.append(v_rad,(wave_weight - local_min_2)/local_min_2 * cst.c.to("km/s").value)


b[0] = fit.params['b0'].value
new_b = fit.params['b0'].stderr * 2 + b[0]
multiplier = 2
#print('calculating new_b')
# issue with std being too small leading to values which are too close together breaking code so this loop
# will increase the number of std's away until the new component is suitably different from the first one
while (new_b - b[0]) < 0.1:
    multiplier += 1
    new_b = fit.params['b0'].stderr * multiplier + b[0]
b = np.append(b,new_b)
#print('b multiplier is:', multiplier)

N[0] = fit.params['N0'].value
new_N = fit.params['N0'].stderr * 2 + N[0]
N = np.append(N,new_N)

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


ax.plot(wavelengths, fit_2.best_fit, c = 'orange', label = '2 components')
#ax.plot(wavelengths, fit_2_comp_1.best_fit, label = 'comp1')
#ax.plot(wavelengths, fit_2_comp_2.best_fit, label = 'comp2')
plt.legend()

fit_2.params.pretty_print()
print('')
fit.params.pretty_print()

plt.show()

