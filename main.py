import numpy as np
from scipy.signal import find_peaks


from edibles.fit import fit
from edibles.create_model import create_line, create_cont
from edibles.functions.load_fits_range import load_fits_range


# file params
star_name = 'HD170740'
file = '/data/DR3_fits/HD170740/RED_860/HD170740_w860_n20_20140915_L.fits'
xmin = 7661.
xmax = 7670.


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

=======
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Continuum
n_points = 4
cont = create_cont(data, n_points)

# ==========
model = cont

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Lines
peak_cutoff = 0.35
prominence = (np.max(flux) - np.min(flux)) * peak_cutoff
peaks, _ = find_peaks(-flux, prominence=prominence)


list_of_lines = []
for i in range(len(peaks)):
    name    = 'line' + str(i)
    lam_0   = wave[peaks[i]]
    b       = 2.6
    d       = .005
    tau_0   = 0.1
    line = create_line(name, lam_0, b, d, tau_0)

    # ===========
    model *= line
    list_of_lines.append(line)


model = fit(star_name, data, model)


for obj in list_of_lines:
    print(obj)
