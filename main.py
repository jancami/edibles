import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from astropy.io import fits

from edibles.functions.continuum_guess import generate_continuum
# from edibles.new_fit.cont_model import Cont1D
# from edibles.new_fit.astro_v_model import AstroVoigt1D
from edibles.new_fit.models import Cont1D, AstroVoigt1D

from sherpa.data import Data1D, DataSimulFit
from sherpa.stats import LeastSq
from sherpa.optmethods import LevMar, MonCar, NelderMead
from sherpa.fit import Fit, SimulFitModel
from sherpa.plot import DataPlot, ModelPlot, FitPlot, SplitPlot

# ===========
# Main script
# ===========

star_name = 'HD170740'
number_of_lines = 4


# set model to zero
model1 = 0
model2 = 0

# set SOME initial params
n_points = 4
n_piece = n_points - 1

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# MODEL 1 BEGIN

file1 = '/data/DR3_fits/HD170740/BLUE_346/HD170740_w346_n6_20160612_B.fits'
xmin1 = 3300.
xmax1 = 3305.

hdu = fits.open(file1)
spec_flux = hdu[0].data
crval1 = hdu[0].header["CRVAL1"]
cdelt1 = hdu[0].header["CDELT1"]
nwave = len(spec_flux)
wave = np.arange(0, nwave, 1)
spec_wave = (wave) * cdelt1 + crval1

# create data subset
min_idx = (np.abs(spec_wave - xmin1)).argmin()
max_idx = (np.abs(spec_wave - xmax1)).argmin()
x1 = spec_wave[min_idx:max_idx]
y1 = spec_flux[min_idx:max_idx]

# Auto-generate continuum guess parameters - no need to do it manually
y_spline1, y_points1= generate_continuum((x1, y1), delta_v=1000, n_piece=n_piece)

# create cont object and define parameters
cont1 = Cont1D()

# always at least 2 points / 1 piece
if n_points >= 1:
    cont1.y1            = y_points1[0]
    cont1.y1.frozen     = False
if n_points >= 2:
    cont1.y2            = y_points1[1]
    cont1.y2.frozen     = False
if n_points >= 3:
    cont1.y3            = y_points1[2]
    cont1.y3.frozen     = False
if n_points >= 4:
    cont1.y4            = y_points1[3]
    cont1.y4.frozen     = False
if n_points >= 5:
    cont1.y5            = y_points1[4]
    cont1.y5.frozen     = False
if n_points >= 6:
    cont1.y6            = y_points1[5]
    cont1.y6.frozen     = False
if n_points >= 7:
    cont1.y7            = y_points1[6]
    cont1.y7.frozen     = False
if n_points >= 8:
    cont1.y8            = y_points1[7]
    cont1.y8.frozen     = False
cont1.n_piece           = n_piece
print(cont1)

# add cont to model
model1 += cont1

# find peaks
peak_cutoff = 0.3
prominence = (np.max(y1) - np.min(y1)) * peak_cutoff
peaks, _ = find_peaks(-y1, prominence=prominence)


# create voigt lines

# params:
b_eff1 = 3.4
Gamma1 = 1.1e+09
scaling1 = 60

b_eff2 = 3.5
Gamma2 = 2.05e+09
scaling2 = 40

# central wavelengths from NIST
c1 = 3302.369
c2 = 3302.979
diff_bw_c1_c2 = c2 - c1

obj1 = AstroVoigt1D()
obj1.cent          = x1[peaks[0]]
obj1.cent.frozen   = False
obj1.b_eff         = b_eff1
# obj1.b_eff.frozen  = True
obj1.Gamma         = Gamma1
# obj1.Gamma.frozen  = True
obj1.scaling       = scaling1
obj1.scaling.frozen = False
print(obj1)

obj2 = AstroVoigt1D()
obj2.cent          = obj1.cent + diff_bw_c1_c2
obj2.cent.frozen   = False
obj2.b_eff         = b_eff2
# obj2.b_eff.frozen  = True
obj2.Gamma         = Gamma2
# obj2.Gamma.frozen  = True
obj2.scaling       = scaling2
obj2.scaling.frozen = False
print(obj2)

# add objects to model
model1 += obj1
model1 += obj2

# MODEL 1 END

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# MODEL 2 BEGIN

file2 = '/data/DR3_fits/HD170740/RED_564/HD170740_w564_n9_20160612_U.fits'
xmin2 = 5885.
xmax2 = 5898.

hdu = fits.open(file2)
spec_flux = hdu[0].data
crval1 = hdu[0].header["CRVAL1"]
cdelt1 = hdu[0].header["CDELT1"]
nwave = len(spec_flux)
wave = np.arange(0, nwave, 1)
spec_wave = (wave) * cdelt1 + crval1

# create data subset
min_idx = (np.abs(spec_wave - xmin2)).argmin()
max_idx = (np.abs(spec_wave - xmax2)).argmin()
x2 = spec_wave[min_idx:max_idx]
y2 = spec_flux[min_idx:max_idx]

# Auto-generate continuum guess parameters - no need to do it manually

y_spline2, y_points2= generate_continuum((x2, y2), delta_v=1000, n_piece=n_piece)

# create cont object and define parameters
cont2 = Cont1D()
# always at least 2 points / 1 piece
if n_points >= 1:
    cont2.y1            = y_points2[0]
    cont2.y1.frozen     = False
if n_points >= 2:
    cont2.y2            = y_points2[1]
    cont2.y2.frozen     = False
if n_points >= 3:
    cont2.y3            = y_points2[2]
    cont2.y3.frozen     = False
if n_points >= 4:
    cont2.y4            = y_points2[3]
    cont2.y4.frozen     = False
if n_points >= 5:
    cont2.y5            = y_points2[4]
    cont2.y5.frozen     = False
if n_points >= 6:
    cont2.y6            = y_points2[5]
    cont2.y6.frozen     = False
if n_points >= 7:
    cont2.y7            = y_points2[6]
    cont2.y7.frozen     = False
if n_points >= 8:
    cont2.y8            = y_points2[7]
    cont2.y8.frozen     = False
cont2.n_piece           = n_piece
print(cont2)

# add cont to model
model2 += cont2

# find peaks  
peak_cutoff = 0.3
prominence = (np.max(y2) - np.min(y2)) * peak_cutoff
peaks, _ = find_peaks(-y2, prominence=prominence)


# create voigt lines

# params:
b_eff3 = 6.0
Gamma3 = 6.0e9
scaling3 = 270

b_eff4 = 5.5
Gamma4 = 6.2e9
scaling4 = 270

# central wavelengths from NIST
c3 = 5889.95095
c4 = 5895.92424


diff_bw_c1_c3 = c3 - c1
diff_bw_c1_c4 = c4 - c1

diff_bw_c3_c4 = c4 - c3

obj3 = AstroVoigt1D()
obj3.cent          = x2[peaks[0]] # obj1.cent + diff_bw_c1_c3
obj3.cent.frozen = False
obj3.b_eff         = b_eff3
obj3.Gamma         = Gamma3
obj3.scaling       = scaling3
obj3.scaling.frozen = False
print(obj3)

obj4 = AstroVoigt1D()
obj4.cent          =  obj3.cent + diff_bw_c3_c4 # x2[peaks[1]]
# obj4.cent.frozen = False
obj4.b_eff         = b_eff4
obj4.Gamma         = Gamma4
obj4.scaling       = scaling4
obj4.scaling.frozen = False
print(obj4)

# add objects to model
model2 += obj3
model2 += obj4

# MODEL 2 END

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# setup to fit / plot

d1 = Data1D('Data 1', x1, y1)
d2 = Data1D('Data 2', x2, y2)

dall = DataSimulFit('combined', (d1, d2))
mall = SimulFitModel('combined', (model1, model2))

# ==========================================
# Initial guesses

    # Dataset 1
dplot1 = DataPlot()
dplot1.prepare(d1)
dplot1.plot()

mplot1 = ModelPlot()
mplot1.prepare(d1, model1)
dplot1.plot()
mplot1.overplot()
plt.show()

    # Dataset 2
dplot2 = DataPlot()
dplot2.prepare(d2)
dplot2.plot()

mplot2 = ModelPlot()
mplot2.prepare(d2, model2)
dplot2.plot()
mplot2.overplot()
plt.show()

# =========================================
# Fitting happens here - don't break please
stat = LeastSq()
opt = LevMar()
print(opt)

vfit = Fit(dall, mall, stat=stat, method=opt)
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
fplot1 = FitPlot()
mplot1.prepare(d1, model1)
fplot1.prepare(dplot1, mplot1)
fplot1.plot()

    # residual
title = 'Data 1'
plt.title(title)
plt.plot(x1, y1-model1(x1))
plt.show()

    # Dataset 2
fplot2 = FitPlot()
mplot2.prepare(d2, model2)
fplot2.prepare(dplot2, mplot2)
fplot2.plot()

    # residual
title = 'Data 2'
plt.title(title)
plt.plot(x2, y2-model2(x2))
plt.show()

    # both datasets - no residuals
splot = SplitPlot()
splot.addplot(fplot1)
splot.addplot(fplot2)

plt.tight_layout()
plt.show()
print()
print()


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Making database

star_info = []
list_of_models = [obj1, obj2, obj3, obj4]
line_params = []

for i in list_of_models:
    each_line = [i.cent.val, i.b_eff.val, i.Gamma.val, i.scaling.val]
    line_params.append(each_line)


star_info = [star_name, number_of_lines, line_params]

print(star_info)

# somehow add the star info to database here
