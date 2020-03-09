# The aim of this script is to fit two of the six telluric lines between 7670 AA and 7685AA
# The general process includes: load the data, cut the spectrum, fit the continuum and
# fits the telluric lines with Voigt profiles

# import the settings
import edibles.src.edibles_spectrum as eS
import edibles.src.model as eM
import numpy as np
import matplotlib.pyplot as plt
from sherpa import models

## Step 1: load and cut the data
filename = '/HD169454/RED_860/HD169454_w860_redl_20160808_O12.fits'
sp = eS.EdiblesSpectrum(filename)
# and let's take a glance at the spectrum
sp.showSpectrum()

# the telluric lines lie in the very red (large wavelength) section so we can drop some data
sp.cutSpectrum(xmin=7667, xmax=7687)
sp.showSpectrum()
# There they are, but it's easier to fit a pair of them at one time, so cut them again
# Let's cut it visually to ~7672.5 to 7681.5
sp.cutSpectrum_visual()
sp.showSpectrum()

## Step 2: fit the continuum:
# before we fit the continuum it is better to mask the telluric region, which might affect the fitting
sp.addMask(n=1)
# fit the continuum with method "fitContinuum", set n=5 and everything else as default, but feel free to explore
sp.fitContinuum(n=5, apply_mask=True)
# you can also use "addMask" to remove some "bad" region from fitting, where n is number of regions
#sp.addMask(n=n)
#sp.fitContinuum(n=5, apply_mask=True)
sp.showSpectrum()

## Step 3: Define the model
# a Voigt profile is set up by four parameters, and it makes everything easier if you can guess the central wavelength
# you can read the plot, or call "searchPeaks", with the number of peaks.
peaks = sp.searchPeaks(n=2)

# the construction of the model begins with a constant continuum level (since it has been normalized)
CstCont = models.Const1D()
# we can initialize the "Sightline" class in edibles_model with this continuum
model_working = eM.Sightline('HD169454', CstCont)
# and add two Voigt lines, named telluric 1 and 2, to the model
model_working.addLine('Telluric1',peaks[0], tau_0 = 0.1)
model_working.addLine('Telluric2',peaks[1], tau_0 = 0.1)
# now this (fairly simple) model is ready
model2fit = model_working.model

## Step 5: model fitting
# import the model to edible_spectrum and do the fitting:
sp.importModel(model2fit)
sp.fitModel()
# while the default setting is recommended, you can still mask some "bad" region:
#sp.addMask(n=n)
#sp.fitModel(apply_mask=True)

#you can withdraw the fitted parameters with "raiseParameter" if the fitting is good
# this method will pop all parameters within all modules that includ the given keywords
telluric_wavelengths = sp.raiseParameter(par="lam")
# a = sp.raiseParameter(module="Telluric1") will give you everything for telluric line #1
# b = sp.raiseParameter(module="Telluric2", par="lam") will give you the lam_0 of telluric line #2

############################################
#additional tips:
# you can use "showSpectrum" and "statusReport" to see what's going on for the data
# with "resetSpectrum" you can start over without reloading the data
# "resetMask" will remove the mask you've set up


