from sherpa.data import Data1D
import numpy as np
import edibles.fit.avoigt as fit
import edibles.fit.make_grid as mg
import matplotlib.pyplot as plt
from sherpa.plot import DataPlot


# Parameters to change
delta_v = 1000
n_piece = 4
n_points = n_piece + 1


dataset = 2

# ===================================================================


if dataset == 1:
    # ============
    # Fake dataset
    # ============

    grid = mg.make_grid(5886, 5899, resolution=80000)
    cont_const = 0.05 * (grid) + 12.3

    yy1 = 1 + fit.voigt(grid, lambda_peak=5890, b_eff=3.47, log_N=12.843, 
                        gamma=6.064e+07, osc_freq=0.631, resolving_power=80000)
    yy2 = 1 + fit.voigt(grid, lambda_peak=5896, b_eff=3.47, log_N=12.843, 
                        gamma=6.098e+07, osc_freq=0.318, resolving_power=80000)

    flux =  cont_const + yy1 #+ yy2
    noise = np.random.normal(0, 0.005, len(flux))
    flux = flux + noise

    # interpolate into fixed grid
    spec_wave = np.linspace(5887, 5898, 1000)
    spec_flux = np.interp(spec_wave, grid, flux)

    wave_subset = spec_wave
    flux_subset = spec_flux


# ====================================================================
if dataset == 2:
    # =======================
    # Real data from HD170740
    # =======================


    # load data from file
    from astropy.io import fits
    hdu = fits.open('/data/DR3_fits/HD170740/RED_860/HD170740_w860_redl_20160613_O12.fits')
    # hdu = fits.open('/data/DR3_fits/HD170740/RED_860/HD170740_w860_redl_20140915_O12.fits')
    spec_flux = hdu[0].data
    crval1 = hdu[0].header["CRVAL1"]
    cdelt1 = hdu[0].header["CDELT1"]
    nwave = len(spec_flux)
    wave = np.arange(0, nwave, 1)
    spec_wave = (wave) * cdelt1 + crval1

    # create data subset
    x_min = 7662.
    x_max = 7670.
    min_idx = (np.abs(spec_wave - x_min)).argmin()
    max_idx = (np.abs(spec_wave - x_max)).argmin()
    wave_subset = spec_wave[min_idx:max_idx]
    flux_subset = spec_flux[min_idx:max_idx]

    # error of data for different stats options
    err = (hdu[0].header["CRDER1"] + hdu[0].header["CSYER1"]) * np.ones_like(flux_subset)

    # normalize data
    flux_subset = flux_subset / (np.max(flux_subset) - np.min(flux_subset)) 
    err = err / (np.max(flux_subset) - np.min(flux_subset))


# =========================================================================

# find peaks of spectrum
from scipy.signal import find_peaks

prominence = (np.max(flux_subset) - np.min(flux_subset)) * 0.1

peaks, _ = find_peaks(-flux_subset, prominence=prominence)

plt.plot(wave_subset, flux_subset)
plt.plot(wave_subset[peaks], flux_subset[peaks], 'x')
plt.show()


# =========================================================================

# CREATE MODEL
from sherpa.astro.optical import AbsorptionVoigt
from cont_model import Cont1D



# create initial continuum guess spline points
from edibles.functions.continuum_guess import generate_continuum
y_spline, y_points= generate_continuum((wave_subset, flux_subset), 
                                        delta_v=delta_v, n_piece=n_piece)

print(np.sum(np.abs(flux_subset - y_spline)))

cont = Cont1D()

# always at least 2 points / 1 piece
if n_points >= 1:
    cont.y1            = y_points[0]
    cont.y1.frozen     = False
if n_points >= 2:
    cont.y2            = y_points[1]
    cont.y2.frozen     = False
if n_points >= 3:
    cont.y3            = y_points[2]
    cont.y3.frozen     = False
if n_points >= 4:
    cont.y4            = y_points[3]
    cont.y4.frozen     = False
if n_points >= 5:
    cont.y5            = y_points[4]
    cont.y5.frozen     = False
if n_points >= 6:
    cont.y6            = y_points[5]
    cont.y6.frozen     = False

cont.n_piece           = n_piece
print(cont)

# ==========
model = cont
# ==========



# voigt_models = [v1, v2, v3, v4, v5, v6, v7, v8]
# voigt_models = []
for i in range(len(peaks)):

    temp = AbsorptionVoigt()
    temp.center        = wave_subset[peaks[i]]  #  7664.87
    temp.center.frozen = False
    temp.ew            = .3
    temp.fwhm          = 6.
    temp.lg            = 2.
    print(temp)

    model += temp
    temp=0



# # telluric on left 
# v1 = AbsorptionVoigt()
# v1.center        = wave_subset[peaks[1]]  #  7664.87
# v1.center.frozen = False
# v1.ew            = .343141
# v1.fwhm          = 7.21035
# v1.lg            = 1.81324
# print(v1)

# 7665.94

# # telluric on right
# v2 = AbsorptionVoigt()
# v2.center        = wave_subset[peaks[2]]  #  7665.94      # v1.center + (7665.94-7664.87)
# v2.center.frozen = False
# v2.ew            = 0.37439
# v2.fwhm          = 6.78904
# v2.lg            = 2.3472
# print(v2)


# # small IS on left
# v3 = AbsorptionVoigt()
# v3.center        = wave_subset[peaks[0]]  #  7664.43
# v3.center.frozen = False
# v3.ew            = 0.191724
# v3.fwhm          = 5.34348
# v3.lg            = 5.91221
# print(v3)

# model = cont + v1 + v2 + v3



d = Data1D('ex', wave_subset, flux_subset)

dplot = DataPlot()
dplot.prepare(d)
dplot.plot()

from sherpa.plot import ModelPlot
mplot = ModelPlot()
mplot.prepare(d, model)
dplot.plot()
mplot.overplot()
plt.show()



from sherpa.stats import LeastSq
stat = LeastSq()
staterr = stat.calc_staterror


from sherpa.optmethods import LevMar
opt = LevMar()
print(opt)


from sherpa.fit import Fit
vfit = Fit(d, model, stat=stat, method=opt)
print(vfit)

vres = vfit.fit()
print(vres.succeeded)

print(vres.format())


from sherpa.plot import FitPlot
fplot = FitPlot()
mplot.prepare(d, model)
fplot.prepare(dplot, mplot)
fplot.plot()


print(np.sum(np.abs(flux_subset - model(d.x))))

title = 'number of spline points: ' + str(n_piece+1)
plt.title(title)


plt.plot(wave_subset, flux_subset-model(wave_subset))
plt.show()








