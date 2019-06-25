import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

from edibles.functions.continuum_guess import generate_continuum


from edibles.models import Cont1D, VoigtAbsorptionLine
from edibles.functions.load_fits_range import load_fits_range

from sherpa.data import Data1D
from sherpa.stats import LeastSq
from sherpa.optmethods import LevMar, MonCar, NelderMead
from sherpa.fit import Fit
from sherpa.plot import DataPlot, ModelPlot, FitPlot

import time

start = time.time()

# ===========
# Main script
# ===========

# file params
star_name = 'HD170740'
file = '/data/DR3_fits/HD170740/RED_860/HD170740_w860_n20_20140915_L.fits'
xmin = 7661.
xmax = 7670.

# spline params
n_points = 4
n_piece = n_points - 1

# line params

peak_cutoff = 0.035

b_1       = 2.6
d_1       = .005
tau_0_2   = 0.08

b_2       = 2.2
d_2       = .006
tau_0_3   = 0.08

tau_0_4 = 0.08

tau_0_0 = 0.01
tau_0_1 = 0.01
tau_0_5 = 0.01
tau_0_6 = 0.01

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# MODEL 1 BEGIN


# Get data
# --------


data = load_fits_range(file, xmin, xmax)

wave, flux = data

    # Auto-generate continuum guess parameters - no need to do it manually



y_spline1, y_points1= generate_continuum((wave, flux), delta_v=1000, n_piece=n_piece)



# create cont object and define parameters
# ----------------------------------------

cont = Cont1D()

    # always at least 2 points / 1 piece
if n_points >= 1:
    cont.y1            = y_points1[0]
    cont.y1.frozen     = False
if n_points >= 2:
    cont.y2            = y_points1[1]
    cont.y2.frozen     = False
if n_points >= 3:
    cont.y3            = y_points1[2]
    cont.y3.frozen     = False
if n_points >= 4:
    cont.y4            = y_points1[3]
    cont.y4.frozen     = False
if n_points >= 5:
    cont.y5            = y_points1[4]
    cont.y5.frozen     = False
if n_points >= 6:
    cont.y6            = y_points1[5]
    cont.y6.frozen     = False
if n_points >= 7:
    cont.y7            = y_points1[6]
    cont.y7.frozen     = False
if n_points >= 8:
    cont.y8            = y_points1[7]
    cont.y8.frozen     = False
# print(cont)


    # add cont to model
model = cont



# create voigt lines
# ------------------

# find peaks
prominence = (np.max(flux) - np.min(flux)) * peak_cutoff
peaks, _ = find_peaks(-flux, prominence=prominence)
# print(peaks)

# plt.plot(wave, flux)
# plt.plot(wave[peaks], flux[peaks], 'x')

# plt.show()

line1 = VoigtAbsorptionLine()
line1.lam_0          = wave[peaks[2]]
line1.b              = b_1
line1.d              = d_1
line1.tau_0          = tau_0_2
# print(line1)

line2 = VoigtAbsorptionLine()
line2.lam_0          = wave[peaks[3]]
line2.b              = b_2
line2.d              = d_2
line2.tau_0          = tau_0_3
# print(line2)

line3 = VoigtAbsorptionLine()
line3.lam_0          = wave[peaks[4]]
line3.b              = b_2
line3.d              = d_2
line3.tau_0          = tau_0_4
# print(line3)



line4 = VoigtAbsorptionLine()
line4.lam_0          = wave[peaks[0]]
line4.b              = 1.5
line4.d              = d_2
line4.tau_0          = tau_0_0
# print(line4)
line5 = VoigtAbsorptionLine()
line5.lam_0          = wave[peaks[1]]
line5.b              = 1.5
line5.d              = d_2
line5.tau_0          = tau_0_1
# print(line5)
line6 = VoigtAbsorptionLine()
line6.lam_0          = wave[peaks[5]]
line6.b              = 1.5
line6.d              = d_2
line6.tau_0          = tau_0_5
# print(line6)
line7 = VoigtAbsorptionLine()
line7.lam_0          = wave[peaks[6]]
line7.b              = 1.5
line7.d              = d_2
line7.tau_0          = tau_0_6
# print(line7)


# multiply lines by model
model *= line1
model *= line2
model *= line3
# model *= line4
# model *= line5
# model *= line6
# model *= line7


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# setup to fit / plot

d = Data1D('2014-09-15', wave, flux)


# ==========================================
# Initial guesses

    # Dataset 1
dplot = DataPlot()
dplot.prepare(d)
# dplot.plot()

mplot = ModelPlot()
mplot.prepare(d, model)
# dplot.plot()
# mplot.overplot()
# plt.show()


# =========================================
# Fitting happens here - don't break please
stat = LeastSq()

opt = NelderMead()
# opt = LevMar()

print(opt)

vfit = Fit(d, model, stat=stat, method=opt)


print(vfit)
vres = vfit.fit()

print()
print()
print('Did the fit succeed? [bool]')
print(vres.succeeded)
print()
print()
print(vres.format())

# =========================================
# Plotting after fit

    # Dataset 1
fplot = FitPlot()
mplot.prepare(d, model)
fplot.prepare(dplot, mplot)
fplot.plot()

    # residual
title = '2014-09-15'
plt.title(title)
plt.plot(wave, flux-model(wave))

# plt.xaxis(fontsize = )
plt.xlabel('Wavelength (AA)', fontsize=12)
plt.ylabel('Flux', fontsize=12)
plt.tick_params(axis='both', labelsize=12)


duration = time.time() - start
print('Time taken: ' + str(duration))


plt.show()
