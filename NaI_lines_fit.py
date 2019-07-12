import numpy as np
from scipy.signal import find_peaks

from edibles.fit import fit, multifit
from edibles.create_model import createKnownVelocityCloud, createCont
from edibles.functions.atomic_line_tool import AtomicLines
from edibles.edibles_spectrum import EdiblesSpectrum
from edibles.edibles_settings import datadir


star_name = 'HD170740'
file1 = '/HD170740/BLUE_346/HD170740_w346_n6_20160612_B.fits'
xmin1 = 3300.
xmax1 = 3305.
ion1_1 = 'Na I'
wave1_1 = 3302.3
ion1_2 = 'Na I'
wave1_2 = 3302.9
file2 = '/HD170740/RED_564/HD170740_w564_n9_20160612_U.fits'
xmin2 = 5885.
xmax2 = 5898.
ion2_1 = 'Na I'
wave2_1 = 5890.0
ion2_2 = 'Na I'
wave2_2 = 5896.0

sp1 = EdiblesSpectrum(datadir+file1)
data1 = sp1.getSpectrum(xmin1,xmax1)
wave1, flux1 = data1
sp2 = EdiblesSpectrum(datadir+file2)
data2 = sp2.getSpectrum(xmin2,xmax2)
wave2, flux2 = data2



# ======================================
# ======================================
# Cont parameters
n_points = 4

# ======================================
# ======================================
# Line parameters
v_cloud = -19
b       = 2
d       = 0.001
N       = 14



AtomicLineList = AtomicLines()
f1_1 = AtomicLineList.get_f_known(ion1_1, wave1_1)
f1_2 = AtomicLineList.get_f_known(ion1_2, wave1_2)
f2_1 = AtomicLineList.get_f_known(ion2_1, wave2_1)
f2_2 = AtomicLineList.get_f_known(ion2_2, wave2_2)
f_list = [f1_1, f1_2, f2_1, f2_2]

lab_wave1_1 = AtomicLineList.getLabWavelength(ion1_1, wave1_1)
lab_wave1_2 = AtomicLineList.getLabWavelength(ion1_2, wave1_2)
lab_wave2_1 = AtomicLineList.getLabWavelength(ion2_1, wave2_1)
lab_wave2_2 = AtomicLineList.getLabWavelength(ion2_2, wave2_2)
lab_list = [lab_wave1_1, lab_wave1_2, lab_wave2_1, lab_wave2_2]


names = ['NaI_3302.3', 'NaI_3302.9', 'NaI_5890', 'NaI_5896']


cloud = createKnownVelocityCloud(name=names, num_lines=4, v_cloud=v_cloud, b=b, 
                                    d=d, N=N, f_known=f_list, lab_lam_0=lab_list)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# MODEL 1 BEGIN

cont1 = createCont(data1, n_points)
# ==========
model1 = cont1

model1 *= cloud[0]
model1 *= cloud[1]

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# MODEL 2 BEGIN

cont2 = createCont(data2, n_points)
# ==========
model2 = cont2

model2 *= cloud[2]
model2 *= cloud[3]

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_list = [data1, data2]
model_list = [model1, model2]

fit_model_list = multifit(star_name, data_list, model_list)
print()
print()




# only printing the model parts we need
lines = []
for model in fit_model_list:
    for line in model:
        if line.name[0] == 'N':
            lines.append(line)


for line in lines:
    print(line.name, line.v_cloud.val, line.b.val, line.d.val, line.N.val)
