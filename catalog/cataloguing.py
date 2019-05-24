from edibles.new_fit.fit_voigt import fitter








# load data from file

# # NaI 5890/5896
# file = '/data/DR3_fits/HD170740/RED_564/HD170740_w564_n9_20160612_U.fits'
# xmin = 5885.
# xmax = 5898.
# scaling = 250.

# NaI 3000
file = '/data/DR3_fits/HD170740/BLUE_346/HD170740_w346_n6_20160612_B.fits'
xmin = 3300.
xmax = 3305.
scaling = 30.

# # KI 7665
# file = '/data/DR3_fits/HD170740/RED_860/HD170740_w860_redl_20140915_O12.fits'
# xmin = 7662
# xmax = 7670
# scaling = 1.


peak_cutoff = 0.5
n_points = 5

b_eff=6.47
Gamma=6.064e9


fitter(file, xmin=xmin, xmax=xmax, peak_cutoff=peak_cutoff, n_points=n_points, 
                b_eff=b_eff, Gamma=Gamma, scaling=scaling)



def catalog_maker(star, num_peaks, data):

    catalog_line = []

    catalog_line.append(star)
    catalog_line.append(num_peaks)

    for i in range(len(num_peaks)):

        catalog_line.append(peak_data)