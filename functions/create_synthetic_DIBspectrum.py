def generate_DIBspectrum(dib_catalog, start, end, step=0.01):

 # TBD: clean up import list.
 import numpy as np
 from StringIO import StringIO
 from matplotlib import *
 import matplotlib.pyplot as plt
 from numpy import loadtxt
 import matplotlib.cm as cm
 import pylab as plt
 from scipy import *
 import scipy.stats
 from scipy.interpolate import splrep,splev
 from scipy.optimize import curve_fit
 import scipy as sp
 import pyfits
 from scipy.ndimage import gaussian_filter
 from scipy.integrate import romb
 from scipy.integrate import simps
 from scipy import integrate
 from scipy import convolve

 l_dib, fwhm, ew = (np.loadtxt(dib_catalog, delimiter=",").T)
 sigma = fwhm / 2.35482

 wave_dib = np.arange(start,end,step) #np.linspace(start, end, step)
 dibspec  = np.ones_like(wave_dib)

 def dib_depth(wave, central, fwhm, ew):
   sigma = fwhm / 2.35482
   depth = (ew/1000.0) / ( sigma * np.sqrt(2.0 * np.pi) )   #ew = depth * sigma * np.sqrt(2.0 * np.pi)
   # CD = W(/1.571 * FWHM) ???
   #dib = 1.0 - ( depth * np.exp( (-1.0*((wave-central)**2.0))/(2.0*(sigma**2.0)) ) )
   dib_depth = -1.0 * ( ( depth * np.exp( (-1.0*((wave-central)**2.0))/(2.0*(sigma**2.0)) ) ) )
   return dib_depth

 for i in np.arange(len(l_dib)):
     dibspec = dibspec + dib_depth(wave_dib,central=l_dib[i],fwhm=fwhm[i],ew=ew[i])


 for (dib in np.range(len(l_dib))):
   min = (l_dib[dib] - (10.0*sigma[dib] )
   max = (l_dib[dib] + (10.0*sigma[dib] )
   idx = where( (wavelength_array >= min) && (wavelength_array <= max) )
   sub_wave_array = wavelength_array[idx]

   for (lambda_index in np.range(len(sub_wave_array))):
        depth = (ew[dib]/(2.5066*sigma[dib])) * exp( (-1.0*((wavelength - l_dib[)**2.0)) / (2.0 (sigma[**2.0)))
        flux_array[lambda_index] = flux_array[lambda_index] - depth

 return wavelength_array, flux_array
