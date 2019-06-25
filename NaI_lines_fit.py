import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

from edibles.functions.continuum_guess import generate_continuum
from edibles.models import Cont1D, VoigtAbsorptionLine
from edibles.functions.load_fits_range import load_fits_range


from sherpa.data import Data1D, DataSimulFit
from sherpa.stats import LeastSq
from sherpa.optmethods import LevMar, MonCar, NelderMead
from sherpa.fit import Fit, SimulFitModel
from sherpa.plot import DataPlot, ModelPlot, FitPlot, SplitPlot


import time
now = time.time()


# ===========
# Main script
# ===========

star_name = 'HD170740'
number_of_lines = 4


file1 = '/data/DR3_fits/HD170740/BLUE_346/HD170740_w346_n6_20160612_B.fits'
xmin1 = 3300.
xmax1 = 3305.

file2 = '/data/DR3_fits/HD170740/RED_564/HD170740_w564_n9_20160612_U.fits'
xmin2 = 5885.
xmax2 = 5898.


# set SOME initial params
n_points = 4
n_piece = n_points - 1


peak_cutoff = 0.3


# params:
b_1 = 3.4
b_2 = 3.5
b_3 = 6.0
b_4 = 5.5


# central wavelengths from NIST
c1 = 3302.369
c2 = 3302.979
c3 = 5889.95095
c4 = 5895.92424
diff_bw_c1_c2 = c2 - c1
diff_bw_c1_c3 = c3 - c1
diff_bw_c1_c4 = c4 - c1
diff_bw_c3_c4 = c4 - c3



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# MODEL 1 BEGIN


data = load_fits_range(file1, xmin1, xmax1)
x1, y1 = data


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
# print(cont1)

# add cont to model
model1 = cont1

# find peaks
prominence = (np.max(y1) - np.min(y1)) * peak_cutoff
peaks, _ = find_peaks(-y1, prominence=prominence)


# create voigt lines

obj1 = VoigtAbsorptionLine()
obj1.lam_0          = x1[peaks[0]]
obj1.lam_0.frozen   = False
# obj1.b              = b_1
# obj1.b.frozen     = False
# obj1.d              = d1
# obj1.d.frozen     = False
# obj1.tau_0          = tau_01
obj1.tau_0.frozen = False
# print(obj1)

obj2 = VoigtAbsorptionLine()
obj2.lam_0          = obj1.lam_0 + diff_bw_c1_c2
obj2.lam_0.frozen   = False
# obj2.b              = obj1.b
# obj2.b.frozen       = False
# obj2.d              = d2
# obj2.d.frozen       = False
obj2.tau_0          = obj1.tau_0
obj2.tau_0.frozen   = False
# print(obj2)



# add objects to model
model1 *= obj1
model1 *= obj2

# MODEL 1 END

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# MODEL 2 BEGIN


data2 = load_fits_range(file2, xmin2, xmax2)

x2, y2 = data2


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
# print(cont2)

# add cont to model
model2 = cont2

# find peaks  
prominence = (np.max(y2) - np.min(y2)) * peak_cutoff
peaks, _ = find_peaks(-y2, prominence=prominence)


# create voigt lines

obj3 = VoigtAbsorptionLine()
obj3.lam_0          = x2[peaks[0]]
obj3.lam_0.frozen   = False
obj3.b              = b_3
# obj3.b.frozen       = True
# obj3.d              = d1
# obj3.d.frozen       = True
# obj3.tau_0          = tau_01
obj3.tau_0.frozen   = False
# print(obj3)

obj4 = VoigtAbsorptionLine()
obj4.lam_0          = obj3.lam_0 + diff_bw_c3_c4
obj4.lam_0.frozen   = False
# obj4.b              = obj3.b
# obj4.b.frozen       = True
# obj4.d              = d2
# obj4.d.frozen       = True
obj4.tau_0          = obj3.tau_0
obj4.tau_0.frozen   = False
# print(obj4)

# add objects to model
model2 *= obj3
model2 *= obj4

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
# dplot1.plot()

mplot1 = ModelPlot()
mplot1.prepare(d1, model1)
# dplot1.plot()
# mplot1.overplot()
# plt.show()

    # Dataset 2
dplot2 = DataPlot()
dplot2.prepare(d2)
# dplot2.plot()

mplot2 = ModelPlot()
mplot2.prepare(d2, model2)
# dplot2.plot()
# mplot2.overplot()
# plt.show()

# =========================================
# Fitting happens here - don't break please
stat = LeastSq()

# opt = LevMar()
opt = NelderMead()

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
# fplot1.plot()

#     # residual
# title = 'Data 1'
# plt.title(title)
# plt.plot(x1, y1-model1(x1))
# plt.show()

    # Dataset 2
fplot2 = FitPlot()
mplot2.prepare(d2, model2)
fplot2.prepare(dplot2, mplot2)
# fplot2.plot()

#     # residual
# title = 'Data 2'
# plt.title(title)
# plt.plot(x2, y2-model2(x2))
# plt.show()

    # both datasets - no residuals
splot = SplitPlot()
splot.addplot(fplot1)
splot.addplot(fplot2)

plt.tight_layout()
plt.show()




# # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # Making database

# star_info = []
# list_of_models = [obj1, obj2, obj3, obj4]
# line_params = []

# for i in list_of_models:
#     each_line = [i.lam_0.val, i.b.val, i.d.val, i.tau_0.val]
#     line_params.append(each_line)


# star_info = [star_name, number_of_lines, line_params]

print()
print()
duration = time.time() - now
print('Time taken: ' + str(duration))
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print('RESULTS FOR ' + star_name)
print('Line #    cent           b             d           tau_0')
print('1         {:.5f}     {:.5f}       {:.5f}     {:.5f}'.format(obj1.lam_0.val, obj1.b.val, obj1.d.val, obj1.tau_0.val))
print('2         {:.5f}     {:.5f}       {:.5f}     {:.5f}'.format(obj2.lam_0.val, obj2.b.val, obj2.d.val, obj2.tau_0.val))
print('3         {:.5f}     {:.5f}       {:.5f}     {:.5f}'.format(obj3.lam_0.val, obj3.b.val, obj3.d.val, obj3.tau_0.val))
print('4         {:.5f}     {:.5f}       {:.5f}     {:.5f}'.format(obj4.lam_0.val, obj4.b.val, obj4.d.val, obj4.tau_0.val))
print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
print()
