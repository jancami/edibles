def generate_DIBspectrum(dib_catalog, start, end, step=0.01):

 import numpy as np
 from StringIO import StringIO
 from matplotlib import *
 from pylab import *
 from numpy import *
 import matplotlib.pyplot as plt
 from numpy import loadtxt
 import matplotlib.cm as cm
 import pylab as plt
 from scipy import *
 import scipy.stats
 from scipy.interpolate import splrep,splev
 import sys
 import os
 from scipy.optimize import curve_fit
 import scipy as sp
 from scipy.ndimage import gaussian_filter
 from scipy.integrate import romb
 from scipy.integrate import simps
 from scipy import integrate
 from scipy import convolve

 l_dib, fwhm_dib, ew_dib = (np.loadtxt(dib_catalog, skiprows=47, usecols=(0,1,2)).T)
 sigma_dib = fwhm_dib / 2.35482

 start = 3000.0
 end   = 11000.0
 #step  = 0.01

 wave_dib = np.arange(start,end,step) #np.linspace(start, end, step)
 dibspec  = np.ones_like(wave_dib)

for (dib,sigma,ew) in zip(l_dib,sigma_dib,ew_dib):
    print dib, sigma, ew
    xmin=dib-(10.0*sigma)
    xmax=dib+(10.0*sigma)
    lambda_index = where( (wave_dib >= xmin) & (wave_dib <= xmax) )
    for idx in lambda_index[0]:
        #print idx
        depth = (ew/1000.0) / ( sigma * np.sqrt(2.0 * np.pi) )
        dib_depth = -1.0 * ( ( depth * np.exp( (-1.0*((wave_dib[idx]-dib)**2.0))/(2.0*(sigma**2.0)) ) ) )
        dibspec[idx] = dibspec[idx] + dib_depth


"""

 def dib_depth(wave, central, fwhm, ew):
   sigma = fwhm / 2.35482
   depth = (ew/1000.0) / ( sigma * np.sqrt(2.0 * np.pi) )   #ew = depth * sigma * np.sqrt(2.0 * np.pi)
   # CD = W(/1.571 * FWHM) ???
   #dib = 1.0 - ( depth * np.exp( (-1.0*((wave-central)**2.0))/(2.0*(sigma**2.0)) ) )
   dib_depth = -1.0 * ( ( depth * np.exp( (-1.0*((wave-central)**2.0))/(2.0*(sigma**2.0)) ) ) )
   return dib_depth

 #for i in np.arange(len(l_dib)):
 #    dibspec = dibspec + dib_depth(wave_dib,central=l_dib[i],fwhm=fwhm[i],ew=ew[i])

 #dibspec = dibspec - len(c60_pos) +1
"""


d=(wave_dib,dibspec)
np.savetxt("DIBspectrum_Hobbs2009_HD183143.dat",np.array(d).T)
