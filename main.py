import numpy as np
from scipy.signal import find_peaks

from edibles.fit import fit
from edibles.create_model import createLine, createKnownLine, createKnownCloud, createCont
from edibles.functions.load_fits_range import load_fits_range
from edibles.functions.find_f_known import find_F


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


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Continuum
n_points = 4
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

# list_of_lines = []
# for i in range(len(peaks)):
#     name    = 'line' + str(i)
#     lam_0   = wave[peaks[i]]
#     b       = 1
#     d       = .005
#     N       = 0.14
#     f_known = 6.82e-01
#     line = createKnownLine(name, lam_0, b, d, N, f_known)
#     # ===========
#     model *= line
#     list_of_lines.append(line)
#     #  fit b params will not be accurate for telluric lines in KI region

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ion0 = 'Na I'
wave0 = wave[peaks[0]]
ion1 = 'Na I'
wave1 = wave[peaks[1]]

f0 = find_F(ion0, wave0)
f1 = find_F(ion1, wave1)




name    = ['line0', 'line1']
lam_0   = [wave[peaks[0]], wave[peaks[1]]]
b       = [2.0, 2.0]
d       = [0.005, 0.005]
N       = [0.14, 0.14]
f_known = [f0, f1]

cloud = createKnownCloud(name=name, num_lines=2, lam_0=lam_0, b=b, d=d, N=N, f_known=f_known)
model *= cloud

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fit_model = fit(star_name, data, model)

for line in cloud:
    print(line)
