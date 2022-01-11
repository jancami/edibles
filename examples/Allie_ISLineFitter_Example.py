# This example scratch is designed for Allie as an example on
# how to use ISLineFitter and work with the Na 3302 doublet in
# some sight lines. It also contains information on how to store
# information into the SightLine class of IsLineFitter 2.0.

#####################################################################
# Step 0, the import statements
#####################################################################
# Standard packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Qt5Agg', force=True)

# EDIBLES related
from edibles.utils.edibles_oracle import EdiblesOracle  # file selection
from edibles.utils.edibles_spectrum import EdiblesSpectrum  # data I/O
from edibles.utils.ContinuumFitter import ContinuumFitter # data normalization
from edibles.utils.SightLineModel import SightLine  # model builder
from edibles.utils.ISLineFitter import ISLineFitter  # auto-fit for Na 3302 doublet


#####################################################################
# Step 1, load data
#####################################################################
# Select the right fits file for the star and wavelength region, read
# in the data, and display it on the screen.
# We will use HD23180 as example here. Feel free to explore other stars

# Step 1.1, load data around 3302
fits_log = EdiblesOracle()
file_all = fits_log.getFilteredObsList(object=["HD 23180"],
                                       OrdersOnly=True,
                                       Wave=3302)

# There are multiple files selected
# use the one for 2017-08-27 and order 11
for filename in file_all:
    if "20170827" in filename and "O11" in filename:
        break

# use EdiblesSpectrum to load fits file and use barycentric data
# we also want to focus around the Na 3302 doublet
sp_3302 = EdiblesSpectrum(filename)
wave_3302, flux_3302 = sp_3302.bary_wave, sp_3302.flux
idx = (wave_3302 > 3301) & (wave_3302 < 3304.5)
wave_3302, flux_3302 = wave_3302[idx], flux_3302[idx]


# Step 1.2, load data around Na D lines
# we can re-use the fits_log
file_all = fits_log.getFilteredObsList(object=["HD 23180"],
                                       OrdersOnly=True,
                                       Wave=5890)
# use the observation from the same day for order 4
for filename in file_all:
    if "20170827" in filename and "O4" in filename:
        break
sp_NaD = EdiblesSpectrum(filename)
wave_NaD, flux_NaD = sp_NaD.bary_wave, sp_NaD.flux

# The Na D lines are more separated, so we will break the data into smaller
# pieces to focus on each of the two lines
idx_D1 = (wave_NaD > 5888.5) & (wave_NaD < 5892)
wave_NaD1, flux_NaD1 = wave_NaD[idx_D1], flux_NaD[idx_D1]

idx_D2 = (wave_NaD > 5894.5) & (wave_NaD < 5898)
wave_NaD2, flux_NaD2 = wave_NaD[idx_D2], flux_NaD[idx_D2]



# Step 1.3, organize the data and take a look
# I do not want to repeat the same piece of code three time, so I put the three
# data segment into a dictionary and I can loop among the three pieces. To access
# the data：
# data_all[segment_name] = [segment_wave, segment_flux], or
# segment_wave = data_all[segment_name][0]
# segment_flux = data_all[segment_name][1]
data_all={"3302Doublet": [wave_3302, flux_3302],
          "Na_D1": [wave_NaD1, flux_NaD1],
          "Na_D2": [wave_NaD2, flux_NaD2]}

# For visualize, since I will reuse the code, it is easier to make it a function

def check_data(data):
    fig = plt.figure(figsize=[12, 6])
    for i, segment_name in enumerate(data.keys()):
        ax = fig.add_subplot(1, 3, i + 1)
        ax.plot(data[segment_name][0],
                data[segment_name][1],
                color="k", label=segment_name)
        ax.grid()
        ax.legend(loc="lower right")
        if i == 0:
            ax.set_ylabel("Raw Flux")
        if i == 1:
            ax.set_xlabel("Wavelength (AA)")
    plt.tight_layout()
    plt.show()

check_data(data_all)



#####################################################################
# Step 2, normalize the spectrum
#####################################################################
# We will use the ContinuumFitter that allows interactive plot
# Here is a how-to:
# 1. Left click on the start and end points of the continuum regions,
#    they are the parts of spectrum without absorption. The code is
#    smart and it is fine to click just "around" the point.
# 2. Right click to undo/remove your last selection.
# 3. Middle click (press the wheel/scroller) to finish selection.
# 4. If you cannot do middle-clicking, make sure the n_regions is the
#    number of regions in your mind. You have to select 2 * n_regions
#    points on the plot.
# 5. The code will ask if you are happy with the result, type Y to continue,
#    otherwise type N to redo.

# An auxiliary function to evaluate if the normalization is good
def make_test_plot(normalizer, continuum, anchor):
    fig = plt.figure(figsize=(10, 6.5))
    plt.gcf().subplots_adjust(hspace=0)
    spec = matplotlib.gridspec.GridSpec(ncols=1, nrows=2,
                             height_ratios=[4, 4])

    # Top panel for raw data and overall fitting
    ax0 = fig.add_subplot(spec[0])
    plt.gca().xaxis.set_visible(False)
    plt.step(normalizer.wave, normalizer.flux, color="0.5")
    plt.scatter(anchor.T[0], anchor.T[1], marker="x", s=80, color="r")
    plt.plot(normalizer.wave, continuum(normalizer.wave), color="orange")
    plt.ylabel("Raw Data")

    # Lower panel for normalized data and multi components
    ax1 = fig.add_subplot(spec[1])
    plt.step(normalizer.wave, normalizer.flux / continuum(normalizer.wave), color="0.5")
    plt.scatter(anchor.T[0], np.ones_like(anchor.T[1]), marker="x", s=80, color="r")
    plt.plot(normalizer.wave, np.ones_like(normalizer.wave), linestyle="--", color="orange")
    ax1.set_ylim([0.97, 1.03])
    plt.ylabel("Normalized Data")
    plt.xlabel('Wavelenght $\AA$')
    plt.show()

# I can repeat the code three times, but it is easier to do a loop among the
# three data segments. I stored the data in a dictionary.

for segment_name in data_all.keys():
    print("Now working on " + segment_name)
    normalizer = ContinuumFitter(wave=data_all[segment_name][0],
                                 flux=data_all[segment_name][1])
    while True:
        cont, anchor = normalizer.SplineManualRegion(n_anchors=6, n_regions=99)
        # cont is the continuum based on your selection.
        # try change n_anchors to see what happens

        make_test_plot(normalizer, cont, anchor)
        # Here we call the auxiliary function to see how good the fitting is.

        # This part of code ask your thoughts and permission to carry on
        while True:
            message = "Are you happy with the normalization and ready to proceed?[Y/N]"
            response = input(message)
            if response.upper() in ["Y", "N"]:
                break
                # ask new input if you did not type N, n, Y, or y

        # if you did not type Y or y, go back to line 155 to redo the fitting.
        if response.upper() == "Y":
            data_all[segment_name][1] = data_all[segment_name][1] / cont(data_all[segment_name][0])
            break

# Let's take a look at the normalized data
# we can use the function we defined earlier, rather than retype the code
check_data(data_all)



#####################################################################
# Step 3, auto-fit the Na-3302 doublet
#####################################################################
# we will use the data from the 3302 region a couple of times, so I decided
# to "take them out" from the dictionary first....
wave2fit, flux2fit = data_all["3302Doublet"][0], data_all["3302Doublet"][1]

# start ISLineFitter by feeding the data and tell it the data is normalized.
fitter = ISLineFitter(wave2fit, flux2fit, verbose=1, normalized=True)

# We can adjust the fitting with the following set-ups:
# What we are trying to fit? Na!
# What is the region we want to use? Everything within wave2fit! (i.e. min to max)
# What criteria do we use to tell if the fit is good? BIC!
best_result = fitter.fit(species="NaI",
                         WaveMin=np.min(wave2fit), WaveMax=np.max(wave2fit),
                         criteria="bic")

# Let's take a look at the result...
print(best_result.fit_report())
fitter.plotModel(which=-2)


print(a)
















# Step 2, build the model
# For a sight line, we can have one or multiple clouds (velocity components), and
# within each cloud, we can have one or multiple species. The SightLine class is
# designed to store the information. We start with a new and empty Sightline.
model_builder = SightLine()

# For your project, we will be focusing on Na. We will build a simple model with just
# one cloud, and compare it to the data you just see. So we add one cloud to the model.
# We will call it Cloud_0 with velocity offset of -5.0 km/s
model_builder.AddCloud(cloud_name="Cloud_0", v_cloud=13.15)
model_builder.AddSpecies(cloud="Cloud_0", species="NaI", known_N=3.945e13, known_b=0.37)

model_builder.AddCloud(cloud_name="Cloud_1", v_cloud=10.99)
model_builder.AddSpecies(cloud="Cloud_1", species="NaI", known_N=5.542e13, known_b=2.91)

model_builder.AddROI(wave, resolution=2.8)
flux_model = model_builder(wave, plot=False)

plt.plot(wave, flux, color="k")
plt.plot(wave, flux_model, color="r")
plt.show()


# Add NaD lines
file_all = fits_log.getFilteredObsList(object=["HD 23180"],
                                       OrdersOnly=True,
                                       Wave=5890)
for filename in file_all:
    if "20170827" in filename and "O4" in filename:
        break

sp2 = EdiblesSpectrum(filename)
wave2, flux2 = sp2.bary_wave, sp2.flux
idx2 = (wave2 > 5887.5) & (wave2 < 5898.5)
wave2, flux2 = wave2[idx2], flux2[idx2]

plt.plot(wave2, flux2)
plt.grid()
plt.show()


normalizer = ContinuumFitter(wave2, flux2)
while True:
    cont, anchor = normalizer.SplineManualRegion(n_anchors=6, n_regions=99)
    make_test_plot(normalizer, cont, anchor)
    while True:
        message = "Are you happy with the normalization and ready to proceed?[Y/N]"
        response = input(message)
        if response.upper() in ["Y", "N"]:
            break
    if response.upper() == "Y":
        flux2 = flux2 / cont(wave2)
        break


model_builder.AddROI(wave2, resolution=2.8)
flux_model2 = model_builder(wave2, plot=False)
plt.plot(wave2, flux2, color="k")
plt.plot(wave2, flux_model2, color="r")
plt.show()
