# The code below measures the ISM velocity using Na doublet at 3302 in HD 49787
# This is a brief demonstration on some functions currently available in the software

import edibles.src.edibles_spectrum as eS
import astropy.constants as cst
import edibles.src.model as eM

# two spectra from HD49787 that will be used, as specified below
# the name of the spectra file can be found in the obs-log.
spec_name1 = "HD49787/BLUE_346/HD49787_w346_blue_20180212_O12.fits"
spec_name2 = "HD49787/BLUE_346/HD49787_w346_blue_20180215_O12.fits"
date = ["20180212","20180213"]

# Before start, please check the edibles_setting.py that records the directories
# Make sure you put the fits file under the datadir
# when a EdiblesSpectrum instance is initialized, it will load the edibles_linelist to a built-in dictionary
# for important species and their wavelengths. You can customize you own list of course.

# Initialize and visualize the data
sp = eS.EdiblesSpectrum([spec_name1,spec_name2], panel_name=date)
sodium_wavelength = sp.linelist["Na"][0:2]
sp.showSpectrum()
# You can also add more spectra by
# sp.loadSpectrum(filename,panel_name="123")

# Bary-correction and cut the spectra, to focus on the sodium doublets
sp.baryCorrection()
sp.cutSpectrum(xmin=3301.5, xmax=3304)
sp.showSpectrum()
# Since we loaded two spectra, i.e. 2 panels in edibles-spectrum, many command can be specified to one or several
# panels if you include e.g., panels = [0,1]. If nothing was given, all panels will be included.

# Continuum fitting but before that, we can mask the lines so the fitting is easier
sp.addMask()
sp.fitContinuum(mode="polynomial", n=3, apply_mask=True, min_sigma=0.5)
sp.resetMask()
sp.showSpectrum()
# if you already know there are two pieces of spectra to be masked, you can use sp.addMask(n=2)
# the continuum fitting is an iteration process, so don't worry if you miss a few point for the masking
# currently there are two continuum-fitting modes: polynomial and spline. You can use "p" and "s" for short
# n is the max-order for polynomial, or the number of points for spoine
# min_sigma is related to "what points are on the continuum", by default it is set to 0.3 but the spectra is a little
# noisy in this case so I set it to 0.5



# Model building

# All lines were assumed to have Voigt profile in optical depth
# There are two types of lines eventually, those with known wavelength whose position is defined by velocity-offset
# The ohter kind of line has "flexible" wavelengths
# The more basic model is a "Cloud" containing several lines with known wavelengths, they share the same ISM velocity,
# and their b parameters, representing the Gaussian width, can be "linked"
# A more complicated model is called "Sightline", that contains several clouds at different velocities, and a list of
# lines with flexible wavelengths
# We are using a Sightline model containing the sodium doublet

# but before create the model, we can estimate the position (v-offset) and strength (tau_0) of the lines
# let's start by convert the spectra into velocity frame
sp.converToVelocity(center=sodium_wavelength[0])
sp.showSpectrum()

# we can also get the numbers using the interacting figures
sp.getPosition_visual(n=4,wavelength=sodium_wavelength,tau_0=True,panels=0)
# n is the number of points to be taken
# wavelength if set, will be used to return relative velocity rather than wavelength
# tau_0, if set, will return the estimated tau_0 (using default b and d) for the model

model_working = eM.Sightline(name="HD49787")
model_working.addLine_toclouds(sodium_wavelength,velocity=28.8, tau_0=[0.005,0.002])
model_working.addLine_toclouds(sodium_wavelength,velocity=34.3, tau_0=[0.007,0.004])
# note since the two components are blending, I somewhat reduced the tau_0 from the earlier reported values
# you can also use clouds_idx, cloud_name to specific add lines to which cloud
# but use velocity (that is dfferent from existing ones in Sightline) can add a new cloud to sightline
# or you can also creat a Cloud instance and import that into sightline
# .addline_tolines will allow you to add flexible wavelength lines


# and we need to consider the instrumental profile
kernel = sp.getKernel(panels=0)
model_working.importInstrumental(kernel)
# The FFT in the convolution could shift the model around (trying to fix though), this effect depends only on the
# number of points of the input data and the kernel so can be monitored.

# note model_working is not the model, we need to compile it, and import it to edibles-spectrum
# Since we have multiple components in this model which makes it more complex, it is better to "fix" the offset from
# convolution. This offset can be calculated as follow, and included in the parameter when building the model. It will
# be corrected again when we read the result. A undesired detour but it works :-(
conv_correction = sp.conv_offset * sp.header[1]["CDELT1"] / sodium_wavelength[0] * cst.c.to('km/s').value
model2fit = model_working.compileModel("link_b", conv_correction=conv_correction)
sp.importModel(model2fit)

# Fit the data and get the result (only v_offset in our case)
sp.fitModel(apply_mask=True,panels=0)
# panels=0 means only use data in panel 0; the dashed orange lines in the residual plot are estimated noise
name,v = sp.raiseParameter(par="v_offset")
for item in v:
    v_correct = item + conv_correction
    print(v_correct)

# since the FFT/convolution introduced offset, the v cannot be directly used so a correction is given here
# and this offset depends on how many data points are there so you have to remake the model for each panel
model_working = eM.Sightline(name="HD49787")
model_working.addLine_toclouds(sodium_wavelength,velocity=28.8, tau_0=[0.005,0.002])
model_working.addLine_toclouds(sodium_wavelength,velocity=34.3, tau_0=[0.007,0.004])
kernel = sp.getKernel(panels=1)
model_working.importInstrumental(kernel)
conv_correction = sp.conv_offset * sp.header[1]["CDELT1"] / sodium_wavelength[0] * cst.c.to('km/s').value
sp.importModel(model_working.compileModel("link_b", conv_correction=conv_correction))

sp.fitModel(apply_mask=True,panels=1)

name,v = sp.raiseParameter(par="v_offset")
for item in v:
    v_correct = item + conv_correction
    print(v_correct)

# and you can see the log for all these process, and output everything:
sp.showLog()
sp.outputResult(filename="DoubleComponent")