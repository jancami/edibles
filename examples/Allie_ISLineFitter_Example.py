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
fitter = ISLineFitter(wave2fit, flux2fit, verbose=0, normalized=True)

# We can adjust the fitting with the following set-ups:
# What we are trying to fit? Na!
# What is the region we want to use? Everything between min and max of wave2fit!
# What criteria do we use to tell if the fit is good? BIC!
auto_result = fitter.fit(species="NaI", WaveMin=np.min(wave2fit), WaveMax=np.max(wave2fit),
                         criteria="bic", n_components_min=1, n_components_max=1)
# Latest update on Jan-18: I update ISLineFitter so you can add min and max number of velocity components
# For this example I'm using 1 for both, that is, ISLineFitter will only try to use 1 velocity component.

# Let's take a look at the result, and load the result to the SightLine class, a data structure
print(auto_result.fit_report())
fitter.plotModel(which=-1)

# Use SightLine class to store the auto-fit parameters
# The result is at auto_result.params["par_name"].value
model_builder = SightLine()
# V_off is stored to a cloud
model_builder.AddCloud(cloud_name="Cloud_0", known_v_cloud=auto_result.params["V_off_Cloud0"].value)
# N and b are with the species
model_builder.AddSpecies(cloud="Cloud_0", species="NaI",
                         known_N=auto_result.params["N_Cloud0"].value,
                         known_b=auto_result.params["b_Cloud0"].value)
# Add resolution information
model_builder.AddROI(data_all["3302Doublet"][0], resolution=2.8)
model_builder.AddROI([data_all["Na_D1"][0], data_all["Na_D1"][0]], resolution=2.8)

# Take a look on what we have in model_builder
# we will define a new function to compare data and model for furture use
def CompareModelData(model_builder, data_all):
    figure = plt.figure(figsize=[8, 4], dpi=200)
    for i, key in enumerate(["3302Doublet", "Na_D1", "Na_D2"]):
        ax = figure.add_subplot(1, 3, i + 1)

        ax.plot(data_all[key][0], data_all[key][1], color="k", label="data")
        ax.plot(data_all[key][0], model_builder(data_all[key][0], plot=False), color="r", label="model")
        ax.set_xlabel(key)
        if i == 0:
            ax.set_ylabel("Normalized Flux")
        ax.grid()
        # ax.legend()
    plt.show()

print(model_builder)
CompareModelData(model_builder, data_all)

#####################################################################
# Step 4, Play with the parameters
#####################################################################
import copy
class TuneParameter():
    """
    This class is ued to collect what is from the fitted model (info stored in model_builder) and what
    do you think that should be changed. The code will combine the two pieces of data and create a new
    model based on your input.

    It is an infinite loop unless you tell it to stop. The loop is like this:
    1. The code ask you if you want to quit. If you type q or Q, it will quit and eventually end the
       entire script. Note your normalization data will be lost too! Other wise, it will proceed.
    2. It will ask you if you want to change one of the four free parameters: column density, velocity
       offset, Gaussian broadening parameter (b), and the resolving power around Na D lines. The code
       will report the current value, if want to change it, type in the new number and hit enter,
       otherwise, hit enter without typing anything. If you did not type a number, the code will ask
       you to do it again.
    3. Although you can type any value for each parameters, there could be some "safe boundaries" to
       think about:
       N: should be positive, and not "too different" from the current value (plus minus ~25% at most)
       b: my guess is 0.1 - about 5.0. Very small b will slow down the code
       v_off: I would say -100 to 100. v_off move around the lines, you would not want to chang it too much
       resolving power: 1 to 15. I doubt if it would be much smaller than 2.8 though.
    """

    def __init__(self, model_builder, data_all):
        self.old_model = model_builder
        self.data_all = data_all
        self.parameter = {}
        self.parameter["v_off"] = self.old_model.clouds["Cloud_0"].v_cloud.value
        self.parameter["b"] = self.old_model.clouds["Cloud_0"].species["NaI"].b.value
        N = self.old_model.clouds["Cloud_0"].species["NaI"].N.value
        N_mag = np.floor(np.log10(N))
        self.parameter["N_mag"] = int(N_mag)
        self.parameter["N"] = N / (10 ** N_mag)
        self.parameter["Resolution"] = self.old_model.ROI.all_regions[1]["Resolution"]

    def main(self):
        while True:
            # At the begining of the loop, check if you want to quit
            message = "Press any key to continue, or Q to quit.\n"
            response = input(message)
            if response.upper() == "Q":
                break


            # make a copy of the fitted model, nothing has changed by far
            new_model = copy.deepcopy(self.old_model)

            # check if you want to use different N, v, b, and resolution
            # the processes are similar though: create a message for you, collect your response,
            # and if you want to change anything, assign the new value to the parameter.

            # N
            message = "Change column density N (in the unit of 10E%i)? Currently %.2f" % \
                      (self.parameter["N_mag"], self.parameter["N"])
            N = self.__parseFloatInput(message)
            if N is not None:
                new_model.clouds["Cloud_0"].species["NaI"].N.value = N * 10**self.parameter["N_mag"]

            # v
            message = "Change velocity offset in km/s? Currently %.2f" % self.parameter["v_off"]
            v_off = self.__parseFloatInput(message)
            if v_off is not None:
                new_model.clouds["Cloud_0"].v_cloud.value = v_off

            # b
            message = "Change broadening parameter b in km/s? Currently %.2f" % self.parameter["b"]
            b = self.__parseFloatInput(message)
            if b is not None:
                new_model.clouds["Cloud_0"].species["NaI"].b.value = b

            # resolution?
            message = "Change resolving power around Na D lines in km/s? Currently %.2f" \
                      % self.parameter["Resolution"]
            resolution = self.__parseFloatInput(message)
            if resolution is not None:
                new_model.ROI.all_regions[1]["Resolution"] = resolution

            # Draw!
            print("Ok, based on your new parameter, the model spectrum looks like this...")
            CompareModelData(new_model, self.data_all)

    def __parseFloatInput(self, message):
        # this method is used to gather your input. The conditions are:
        # 1. you type nothing and it will return None, the code will not change anything
        # 2. you type a number and the code will use it to update the parameter value
        # 3. you type something that does not look like a number, the code will ask you to input again.

        message = message + "\n"
        while True:
            response = input(message)
            if response == "":
                response = None
                break
            else:
                try:
                    response = float(response)
                    break
                except:
                    print("Invalid input, please type a float number")

        return response


# Everything is coded inside the TuneParameter class, please see the comment within
turner = TuneParameter(model_builder, data_all)
turner.main()
# note when you type q to end the session, you finish the script and lose the normalized data
