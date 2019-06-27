import numpy as np
from scipy.signal import find_peaks


from edibles.fit import fit
from edibles.create_model import create_line, create_cont
from edibles.functions.load_fits_range import load_fits_range


# file params
star_name = 'HD170740'
file = '/data/DR3_fits/HD147084/BLUE_346/HD147084_w346_n10_20160419_B.fits'
xmin = 3300.
xmax = 3305.


data = load_fits_range(file, xmin, xmax)
wave, flux = data


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
