from edibles.fit.fit import fit
from edibles.fit.models.create_model import *
from edibles.functions.atomic_line_tool import AtomicLines
from edibles.functions.edibles_spectrum import EdiblesSpectrum
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
ion0 = 'Na I'
wave0 = 3302.3
ion1 = 'Na I'
wave1 = 3302.9

sp = EdiblesSpectrum(datadir+"/HD170740/BLUE_346/HD170740_w346_n20_20140916_B.fits")
data = sp.getSpectrum(xmin,xmax)
wave, flux = data

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Continuum
n_points = 3
cont = createCont(data, n_points)

# ==========
model = cont

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


AtomicLineList = AtomicLines()
f0 = AtomicLineList.get_f_known(ion0, wave0)
f1 = AtomicLineList.get_f_known(ion1, wave1)

lab_wave0 = AtomicLineList.getLabWavelength(ion0, wave0)
lab_wave1 = AtomicLineList.getLabWavelength(ion1, wave1)
lab_list = [lab_wave0, lab_wave1]

name    = ['line0', 'line1']
v_cloud = 18
b       = 2
d       = 0.005
N       = 21
f_known = [f0, f1]

cloud = createKnownVelocityCloud(name=name, num_lines=2, v_cloud=v_cloud, b=b, d=d, N=N, f_known=f_known, lab_lam_0=lab_list)

for line in cloud:
    model *= line

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fit_model = fit(star_name, data, model)

for line in model:
    print(line)
