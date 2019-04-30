from __future__ import print_function
import os, glob, sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import wofz
from scipy.stats import chi2
from scipy import interpolate, signal, stats
import warnings
import edibles.fit.initial_suggest as ini_guess
import edibles.fit.err_est as err
import edibles.fit.mpfit_3 as mpfit
import edibles.fit.cont_est as cont_est
import edibles.fit.line_properties as line_properties
from edibles.functions.astro_wrapper import voigt_astro
import edibles.fit.ref_index as ref_index
warnings.simplefilter('ignore')







def cont_func(p, fjac=None, x=None, y=None, err=None):



    # Parameter values are passed in "p"
    # If fjac==None then partial derivatives should not be
    # computed.  It will always be None if MPFIT is called with default
    # flag.
    model = F(x, p)
    # Non-negative status value means MPFIT should continue, negative means
    # stop the calculation.
    status = 0
    return [status, (y-model)/err]


   # EXAMPLE

# import mpfit
x = np.arange(100, dtype=np.float64)
p0 = [5.7, 2.2, 500., 1.5, 2000.]
y = ( p[0] + p[1]*[x] + p[2]*[x**2] + p[3]*np.sqrt(x) + p[4]*np.log(x))
fa = {'x':x, 'y':y, 'err':err}
m = mpfit(cont_func, p0, functkw=fa)
print('status = ', m.status)
if (m.status <= 0):
    print('error message = ', m.errmsg)
print('parameters = ', m.params)