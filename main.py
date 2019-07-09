import numpy as np
from scipy.signal import find_peaks
import astropy.constants as cst

from edibles.fit import fit
from edibles.create_model import *
from edibles.functions.load_fits_range import load_fits_range
from edibles.functions.find_f_known import AtomicLines
from edibles.edibles_spectrum import edibles_spectrum
from edibles.edibles_settings import *


# file params

# star_name = 'HD170740'
# file = '/data/DR3_fits/HD170740/RED_860/HD170740_w860_n20_20140916_L.fits'
# xmin = 7662.
# xmax = 7670.

star_name = 'HD170740'
file = '/data/DR3_fits/HD170740/BLUE_346/HD170740_w346_n20_20140916_B.fits'
xmin = 3301.
xmax = 3305.


data = load_fits_range(file, xmin, xmax)
wave, flux = data


# sp = edibles_spectrum(datadir+"/HD170740/BLUE_346/HD170740_w346_n20_20140916_B.fits")
# print("Barycentric Velocity is", sp.v_bary)
# wave,flux = sp.GetSpectrum()
# plt.plot(wave, flux)
# axes = plt.gca()
# axes.set_xlim([7660,7705])
# axes.set_ylim([0,160])
# plt.vlines((7667.021,7701.093), 0, 160, linestyles='dashed', colors='r')
# plt.show()









# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Continuum
n_points = 5
cont = createCont(data, n_points)

# ==========
model = cont

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Lines
peak_cutoff = 0.3
prominence = (np.max(flux) - np.min(flux)) * peak_cutoff
peaks, _ = find_peaks(-flux, prominence=prominence)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# list_of_lines = []
# for i in range(len(peaks)):
#     name    = 'line' + str(i)
#     lam_0   = wave[peaks[i]]
#     b       = 2.6
#     d       = .005
#     tau_0   = 0.1
#     line = createLine(name, lam_0, b, d, tau_0)
#     # ===========
#     model *= line
#     list_of_lines.append(line)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# lab_list = [3302.369, 3302.978]


# wave = (wave - lab_list[0])/wave * cst.c.to('km/s').value

# data = wave, flux


# v_cloud_list = [1, 50]
# lam_0 = lab_list[0] / (1. - v_cloud_list[0]/cst.c.to('km/s').value)
# # print(v_cloud)
# print(lam_0)

# list_of_lines = []
# for i in range(len(peaks)):
#     name    = 'line' + str(i)
#     v_cloud = v_cloud_list[i]
#     lab_lam_0 = lab_list[0]
#     b       = 1
#     d       = .005
#     N       = 0.1
#     f_known = 6.82e-01
#     line = createKnownVelocityLine(name, v_cloud, b, d, N, f_known, lab_lam_0)
#     # ===========
#     model *= line
#     list_of_lines.append(line)
#     #  fit b params will not be accurate for telluric lines in KI region

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# ion0 = 'Na I'
# wave0 = wave[peaks[0]]
# ion1 = 'Na I'
# wave1 = wave[peaks[1]]

# AtomicLineList = AtomicLines()
# f0 = AtomicLineList.get_f_known(ion0, wave0)
# f1 = AtomicLineList.get_f_known(ion1, wave1)

# name    = ['line0', 'line1']
# lam_0   = [wave[peaks[0]], wave[peaks[1]]]
# b       = [2.0, 2.0]
# d       = [0.005, 0.005]
# N       = [0.14, 0.14]
# f_known = [f0, f1]

# cloud = createKnownCloud(name=name, num_lines=2, lam_0=lam_0, b=b, d=d, N=N, f_known=f_known)
# model *= cloud

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ion0 = 'Na I'
wave0 = wave[peaks[0]]
ion1 = 'Na I'
wave1 = wave[peaks[1]]

AtomicLineList = AtomicLines()
f0 = AtomicLineList.get_f_known(ion0, wave0)
f1 = AtomicLineList.get_f_known(ion1, wave1)

name    = ['line0', 'line1']
v_cloud = 20
b       = [2.0, 2.0]
d       = [0.005, 0.005]
N       = 0.14
f_known = [f0, f1]
lab_list = [3302.369, 3302.978]

cloud = createKnownVelocityCloud(name=name, num_lines=2, v_cloud=v_cloud, b=b, d=d, N=N, f_known=f_known, lab_lam_0=lab_list)
model *= cloud

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fit_model = fit(star_name, data, model)

for line in model:
    print(line)
