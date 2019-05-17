from edibles.new_fit.fit_voigt import fitter

# test 1    NaI 5890/5896
# =================================================

# file = '/data/DR3_fits/HD170740/RED_564/HD170740_w564_n9_20160612_U.fits'
# xmin = 5885.
# xmax = 5898.
# scaling = 100.


# test 2    NaI 3300
# =================================================

# file = '/data/DR3_fits/HD170740/BLUE_346/HD170740_w346_n6_20160612_B.fits'
# xmin = 3300.
# xmax = 3305.
# scaling = 50.


# test 3    KI 7665
# =================================================

file = '/data/DR3_fits/HD170740/RED_860/HD170740_w860_redl_20140915_O12.fits'
xmin = 7662
xmax = 7670
scaling = 1.


peak_cutoff = 0.5
n_points = 5
alpha = 0.05
gamma = 0.005




fitter(file, xmin=xmin, xmax=xmax, peak_cutoff=peak_cutoff, n_points=n_points, 
                alpha=alpha, gamma=gamma, scaling=scaling)

