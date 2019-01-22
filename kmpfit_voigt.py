#!/usr/bin/env python
#------------------------------------------------------------
# Script which demonstrates how to find the best-fit
# parameters of a Voigt line-shape model
# 
# Vog, 26 Mar 2012
#------------------------------------------------------------
import numpy
from matplotlib.pyplot import figure, show, rc
from scipy.special import wofz
from kapteyn import kmpfit
ln2 = numpy.log(2)

def voigt(x, y):
   # The Voigt function is also the real part of 
   # w(z) = exp(-z^2) erfc(iz), the complex probability function,
   # which is also known as the Faddeeva function. Scipy has 
   # implemented this function under the name wofz()
   z = x + 1j*y
   I = wofz(z).real
   return I

def Voigt(nu, alphaD, alphaL, nu_0, A, a=0, b=0):
   # The Voigt line shape in terms of its physical parameters
   f = numpy.sqrt(ln2)
   x = (nu-nu_0)/alphaD * f
   y = alphaL/alphaD * f
   backg = a + b*nu 
   V = A*f/(alphaD*numpy.sqrt(numpy.pi)) * voigt(x, y) + backg
   return V

def funcV(p, x):
    # Compose the Voigt line-shape
    alphaD, alphaL, nu_0, I, a, b = p
    return Voigt(x, alphaD, alphaL, nu_0, I, a, b)

def funcG(p, x):
   # Model function is a gaussian
   A, mu, sigma, zerolev = p
   return( A * numpy.exp(-(x-mu)*(x-mu)/(2*sigma*sigma)) + zerolev )

def residualsV(p, data):
   # Return weighted residuals of Voigt
   x, y, err = data
   return (y-funcV(p,x)) / err

def residualsG(p, data):
   # Return weighted residuals of Gauss
   x, y, err = data
   return (y-funcG(p,x)) / err


# Data from simulated MUSE cube
x = numpy.array([854.05,854.18,854.31,854.44,854.57,854.7,854.83,854.96,\
                 855.09,855.22,855.35,855.48,855.61,855.74,855.87,856.0,\
                 856.13,856.26,856.39,856.52,856.65,856.78,856.91])
y = numpy.array([6.31683382764,6.41273839772,6.43047296256,6.37437933311,\
                 6.34883451462,6.30711287633,6.24409954622,6.09241716936,\
                 5.75421549752,5.20381929725,4.18020502292,3.64663145132,\
                 4.25251198746,5.23945118487,5.76701752096,6.06587703526,\
                 6.15751018003,6.25985588506,6.35063433647,6.41795488447,\
                 6.42002335563,6.35883554071,6.36915982142])
N = len(y)
err = numpy.ones(N)
A = -2
alphaD = 0.5
alphaL = 0.5
a = 6
b = 0
nu_0 = 855
p0 = [alphaD, alphaL, nu_0, A, a, b]

# Do the fit
fitter = kmpfit.Fitter(residuals=residualsV, data=(x,y,err))
fitter.parinfo = [{}, {}, {}, {}, {}, {'fixed':True}]  # Take zero level fixed in fit
fitter.fit(params0=p0)

print ("\n========= Fit results Voigt profile ==========")
print ("Initial params:", fitter.params0)
print ("Params:        ", fitter.params)
print ("Iterations:    ", fitter.niter)
print ("Function ev:   ", fitter.nfev )
print ("Uncertainties: ", fitter.xerror)
print ("dof:           ", fitter.dof)
print ("chi^2, rchi2:  ", fitter.chi2_min, fitter.rchi2_min)
print ("stderr:        ", fitter.stderr)   
print ("Status:        ", fitter.status)

alphaD, alphaL, nu_0, I, a_back, b_back = fitter.params
c1 = 1.0692
c2 = 0.86639
hwhm = 0.5*(c1*alphaL+numpy.sqrt(c2*alphaL**2+4*alphaD**2))
print ("\nFWHM Voigt profile:     ", 2*hwhm)
f = numpy.sqrt(ln2)
Y = alphaL/alphaD * f
amp = I/alphaD*numpy.sqrt(ln2/numpy.pi)*voigt(0,Y)
print ("Amplitude Voigt profile:", amp)
print ("Area under profile:     ", I)

# Fit the Gaussian model
p0 = [-3, 855, 0.5, 6.3]
fitterG = kmpfit.Fitter(residuals=residualsG, data=(x,y,err))
#fitterG.parinfo = [{}, {}, {}, {}, {}]  # Take zero level fixed in fit
fitterG.fit(params0=p0)
print ("\n========= Fit results Gaussian profile ==========")
print ("Initial params:", fitterG.params0)
print ("Params:        ", fitterG.params)
print ("Iterations:    ", fitterG.niter)
print ("Function ev:   ", fitterG.nfev )
print ("Uncertainties: ", fitterG.xerror)
print ("dof:           ", fitterG.dof)
print ("chi^2, rchi2:  ", fitterG.chi2_min, fitterG.rchi2_min)
print ("stderr:        ", fitterG.stderr   )
print ("Status:        ", fitterG.status)

fwhmG = 2*numpy.sqrt(2*numpy.log(2))*fitterG.params[2]
print ("FWHM Gaussian: ", fwhmG)

# Plot the result
rc('legend', fontsize=6)
fig = figure()
frame1 = fig.add_subplot(1,1,1)
xd = numpy.linspace(x.min(), x.max(), 200)
frame1.plot(x, y, 'bo', label="data")
label = "Model with Voigt function"
frame1.plot(xd, funcV(fitter.params,xd), 'g', label=label)
label = "Model with Gaussian function"
frame1.plot(xd, funcG(fitterG.params,xd), 'm', ls='--', label=label)
offset = a_back+b_back*nu_0
frame1.plot((nu_0-hwhm,nu_0+hwhm), (offset+amp/2,offset+amp/2), 'r', label='fwhm')
frame1.plot(xd, a_back+b_back*xd, "y", label='Background')
frame1.set_xlabel("$\\nu$")
frame1.set_ylabel("$\\phi(\\nu)$")
vals = (fitter.chi2_min, fitter.rchi2_min, fitter.dof)
title = "Profile data with Voigt- vs. Gaussian model"
frame1.set_title(title, y=1.05)
frame1.grid(True)
leg = frame1.legend(loc=3)
show()