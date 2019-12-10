import numpy as np
from astro.utilities import percentile,indexnear
from astro.spec import plotlines
import inspect
import minuit
import matplotlib.pyplot as pl
import lm

Ckms = 299792.458                 # speed of light km/s, exact
wlya = 1215.6701

def polyfitr(x, y, order=2, clip=6, xlim=None, ylim=None,
             mask=None, debug=False):
    """ Fit a polynomial to data, rejecting outliers.

    Fits a polynomial f(x) to data, x,y.  Finds standard deviation of
    y - f(x) and removes points that differ from f(x) by more than
    clip*stddev, then refits.  This repeats until no points are
    removed.

    Inputs
    ------
    x,y:
        Data points to be fitted.  They must have the same length.
    order: int (2)
        Order of polynomial to be fitted.
    clip: float (6)
        After each iteration data further than this many standard
        deviations away from the fit will be discarded.
    xlim: tuple of maximum and minimum x values, optional
        Data outside these x limits will not be used in the fit.
    ylim: tuple of maximum and minimum y values, optional
        As for xlim, but for y data.
    mask: sequence of pairs, optional
        A list of minimum and maximum x values (e.g. [(3, 4), (8, 9)])
        giving regions to be excluded from the fit.
    debug: boolean, default False
        If True, plots the fit at each iteration in matplotlib.

    Returns
    -------
    coeff, x, y:
        x, y are the data points contributing to the final fit. coeff
        gives the coefficients of the final polynomial fit (use
        np.polyval(coeff,x)).

    Examples
    --------
    >>> x = np.linspace(0,4)
    >>> np.random.seed(13)
    >>> y = x**2 + np.random.randn(50)
    >>> coeff, x1, y1 = polyfitr(x, y)
    >>> np.allclose(coeff, [1.05228393, -0.31855442, 0.4957111])
    True
    >>> coeff, x1, y1 = polyfitr(x, y, order=1, xlim=(0.5,3.5), ylim=(1,10))
    >>> np.allclose(coeff, [3.23959627, -1.81635911])
    True
    >>> coeff, x1, y1 = polyfitr(x, y, mask=[(1, 2), (3, 3.5)])
    >>> np.allclose(coeff, [1.08044631, -0.37032771, 0.42847982])
    True
    """

    x = np.asanyarray(x)
    y = np.asanyarray(y)
    isort = x.argsort()
    x, y = x[isort], y[isort]

    keep = np.ones(len(x), bool)
    if xlim is not None:
        keep &= (xlim[0] < x) & (x < xlim[1])
    if ylim is not None:
        keep &= (ylim[0] < y) & (y < ylim[1])
    if mask is not None:
        badpts = np.zeros(len(x), bool)
        for x0,x1 in mask:
            badpts |=  (x0 < x) & (x < x1)
        keep &= ~badpts

    x,y = x[keep], y[keep]
    if debug:
        fig = pl.figure()
        ax = fig.add_subplot(111)
        ax.plot(x,y,'.')
        ax.set_autoscale_on(0)
        pl.show()

    coeff = np.polyfit(x, y, order)
    if debug:
        pts, = ax.plot(x, y, '.')
        poly, = ax.plot(x, np.polyval(coeff, x), lw=2)
        pl.show()
        raw_input('Enter to continue')
    norm = np.abs(y - np.polyval(coeff, x))
    stdev = np.std(norm)
    condition =  norm < clip * stdev
    y = y[condition]
    x = x[condition]
    while norm.max() > clip * stdev:
        if len(y) < order + 1:
            raise Exception('Too few points left to fit!')
        coeff = np.polyfit(x, y, order)
        if debug:
            pts.set_data(x, y)
            poly.set_data(x, np.polyval(coeff, x))
            pl.show()
            raw_input('Enter to continue')
        norm = np.abs(y - np.polyval(coeff, x))
        stdev = norm.std()
        condition =  norm < clip * stdev
        y = y[condition]
        x = x[condition]

    return coeff,x,y

def wleastsq(x, y, ysig=None):
    """ Calculate the line of best fit with weights.

    Input
    -----
      x : sequence of floats
          input x data
      y : sequence of floats, len(x)
          input y data, f(x).
      ysig : sequence of floats, len(x), optional
         If the y none sigma errors are given, points are weighted
         by their inverse variance.

    Returns
    -------
      (b, a),(b_sigma, a_sigma)
        The fitted parameters and their one sigma errors.  The fitted
        line has equation y = a + b*x.
    """
    if ysig == None:
        ysig = np.ones(len(x))
    yvar = ysig * ysig   # variance

    s = np.sum(1. / yvar)
    sx = np.sum(x / yvar)
    sxx = np.sum(x*x / yvar)
    sy = np.sum(y / yvar)
    sxy = np.sum(x*y / yvar)

    # See NR section 15.2 for a derivation of the below solutions for
    # the best fit values of a and b.
    # 
    # y = a + b*x 

    temp = s*sxx - sx*sx
    a = (sxx*sy - sx*sxy) / temp
    b = (s*sxy - sx*sy) / temp
    sig_a = np.sqrt(sxx / temp)
    sig_b = np.sqrt(s / temp)

    return (b,a),(sig_b,sig_a)


def fitgauss(x, y, sig, fwhm, x0, height, fix=None):
    """ Fits a gaussian to data with errors using chi**2 minimisation.

    Inputs
    ------
    x, y:
        data values.
    sig:
        One sigma errors in y data values.
    fwhm, x0, height:
        Initial guess values for the gaussian parameters.
    fix: optional
        Name of parameters to fix (e.g. 'height x0').

    Returns
    -------
    results:
        Best fitting height, fwhm and x position, and the Minuit class
        structure which contains all the info about the fit
        (covariance matrix etc).
    """

    fwhm_on_sigma = 2*np.sqrt(2*np.log(2))
    x,y,sig = map(np.asarray, (x,y,sig))
    def g(x, fwhm, x0, height):
        """ Gaussian."""
        sigma = fwhm / fwhm_on_sigma
        return height * np.exp(-0.5 * ((x-x0) / sigma) ** 2)

    guesses = fwhm, x0, height
    res = minchi2(x, y, sig, g, guesses, fix=fix)
    print 'Reduced chi2 = %s' % (res.findchisq(res.par) / (len(x) - 3))
    return res.par[0], res.par[1], res.par[2], res

def fitgaussold(x, y, sig, fwhm, x0, height):
    """ Fits a gaussian to data with errors using chi**2 minimisation.

    Input
    -----
    x:
        x values.
    y:
        y = f(x) data values.
    sig:
        One sigma errors in data values.
    fwhm, x0, height:
        Initial guess values of the gaussian parameters.

    Returns
    -------
    a,fwhm,c,m:
        Best fitting height, fwhm and x position, and the Minuit class
        structure which contains all the info about the fit
        (covariance matrix etc).
    """

    fwhm_on_sigma = 2*np.sqrt(2*np.log(2))
    x,y,sig = map(np.asarray,(x,y,sig))
    def g(x,a,b,c):
        """ Gaussian."""
        return a*np.exp(-(x-c)**2/(2*b*b))

    def chi2(a,b,c):
        ymodel = g(x,a,b,c)
        chi2_i = (ymodel - y) / sig
        return (chi2_i*chi2_i).sum()

    guesses = dict(a=height, b=fwhm/fwhm_on_sigma, c=x0)
    m = minuit.Minuit(chi2)
    m.values.update(guesses)
    m.migrad()         # find parameter values that minimise chi2
    m.hesse()
    print 'Reduced chi2 = %s' % (m.fval / (len(x) - 3))
    a,b,c = (m.values[k] for k in 'abc')
    return a, b*fwhm_on_sigma, c, m

def fitfunc(x, y, sig, func, guess, fixed=None, debug=False):
    """ Fits a model to data with errors using chi**2 minimisation
    and minuit's migrad minimisation routine.

    Input
    -----
    x:
        independent data values
    y:
        y = f(x) dependent data values. Same length as x.
    sig:
        One sigma errors in data values. Same length as x and y.
    func:
        The model to fit to the data. The first function argument must
        be the sequence of x values, and then one argument for each
        parameter in the model.
    guess:
        Dictionary giving the starting guesses for the model
        parameters. Keys are the parameter arguments from func.

    Returns
    -------
    names, values, m, nchi2:  
        The best fitting parameter values and their names (same order
        as they appear as arguments in the input function), the Minuit
        class object containing all the info about the fit and the
        chi^2 per degree of freedom (degrees of freedom = number of
        data points - number of parameters fitted in the model) . The
        minuit object can be used to find the covariance matrix and
        errors, adjust starting guesses, step-sizes, stopping criteria
        and re-run the fit.

    Examples
    --------
    Fit a gaussian to some simulated data. First generate the data.

    >>> x = np.linspace(-5,5,100)
    >>> def g(x,a,x0,sigma):
    ...    return a*np.exp(-0.5*((x-x0)/sigma)**2)
    ...
    >>> np.random.seed(101)
    >>> y = g(x,5,0,1) + np.random.randn(100)*0.5
    >>> sig = np.ones(100)*0.5

    Plot data and the model used to generate data.

    >>> fig = pl.figure()
    >>> p = pl.plot(x, g(x,5,0,1), label='original model')
    >>> p = pl.plot(x, y, label='fake data')

    Now fit with fitfunc and plot the fit. Give starting guesses close
    to the expected values.

    >>> guess = dict(a=8, x0=0.5, sigma=2)
    >>> names,vals,m,nchi2 = fitfunc(x, y, sig, g, guess, debug=True)  # doctest: +ELLIPSIS
    Parameters: a, x0, sigma
    Initial guesses: {...}
    Reduced chi2 = 1.1026136...
    >>> p = pl.plot(x, g(x, *vals), label='fitted model')
    >>> p = pl.legend()

    m.errors gives the 1 sigma errors calculated from the covariance
    matrix. m.covariance gives all the covariance terms.

    Print the best fitting values and errors.

    >>> print zip(names,vals,[m.errors[k] for k in names])  # doctest: +SKIP

    Plot the normalised covariance matrix (AKA correlation matrix).

    >>> p = pl.figure()
    >>> p = pl.pcolor(np.array(m.matrix(correlation=True)),vmin=-1,vmax=1)
    >>> p = pl.xlabel('a, x0, sigma')
    >>> p = pl.ylabel('a, x0, sigma')
    >>> p = pl.colorbar()

    You can see the parameter values are not completely
    independent. The height and width are anti-correlated, which means
    an increase in height tends to be accompanied by a decrease in
    width, and vice-versa.
    """
    x,y,sig = map(np.asarray,(x,y,sig))

    def chi2(*args):
        ymodel = func(x, *args)
        chi2_i = (ymodel - y) / sig
        return (chi2_i*chi2_i).sum()

    parnames = inspect.getargs(func.func_code)[0][1:]
    temp = ', '.join(parnames)
    if debug: print 'Parameters:', temp
    m = minuit.Minuit(eval('lambda %s: chi2(%s)' % (temp,temp), locals()))
    m.values.update(guess)
    if debug: print 'Initial guesses: %s' % m.values
    if fixed is not None:
        if debug: print 'Fixed values: %s' % m.fixed
        m.fixed.update(fixed)
    m.migrad()         # find parameter values that minimise chi2
    m.hesse()          # calculate covariance matrix

    # this corrects for case where x is the output of np.meshgrid,
    # i.e. when fitting a 2-d function.
    #print x.size, x.shape
    if len(x.shape) == 3:
        npts = x.size / 2.0
    else:
        npts = len(x)
    nfittedpar = len(parnames)
    if fixed is not None:
        nfittedpar -= len(x for x in fixed.values() if x)
    nchi2 = m.fval / (npts - nfittedpar)
    if debug:  print 'Reduced chi2 = %s' % nchi2
    bestfitpar = tuple(m.values[k] for k in parnames)
    return parnames, bestfitpar, m, nchi2

def fitfunc_r(x, y, func, guess, fixed=None, clip=6, debug=False):
    """ Fits a function to data without errors, removing outliers.

    func  : Function to be fitted. Must have signature
            func(x, parameter1, parameter2, ...)
            (first parameter needn't be called x, but it must be the
            ordinate values)
    guess : Dictionary mapping parameter names to guesses
    fixed : Dictionary mapping parameter names to True or False

    returns:
      #parameter names in the same order as function signature

      best fitting parameter values (in the same order as they appear
      in func signature)

      #minuit fitting object
      #remaining x values fitted
      #remaining y values fitted

    """
    x,y = map(np.asarray,(x,y))

    if debug:
        a = pl.axes()
        pl.plot(x,y,'.')
        a.set_autoscale_on(0)
        pl.draw()

    def squares(*args):
        ymodel = func(x,*args)
        diff = ymodel - y
        return (diff*diff).sum()

    parnames = inspect.getargs(func.func_code)[0][1:]
    temp = ', '.join(parnames)
    print 'Parameters:', temp
    m = minuit.Minuit(eval('lambda %s: squares(%s)' % (temp,temp), locals()))
    m.values.update(guess)
    print 'Initial guesses: %s' % m.values
    if debug:
        pts, = pl.plot(x, y, '.')
        fit, = pl.plot(x, func(x,*[m.values[k] for k in parnames]), lw=2)
        pl.draw()
        raw_input('Enter to continue')

    if fixed is not None:
        m.fixed.update(fixed)
    m.migrad()         # find parameter values that minimise chi2

    bestfitpar = tuple(m.values[k] for k in parnames)

    norm = np.abs(y - func(x, *bestfitpar))
    stdev = np.std(norm)
    condition =  norm < clip * stdev
    y = y[condition]
    x = x[condition]
    while norm.max() > clip * stdev:
        if len(y) < len(parnames) + 1:
            raise Exception('Too few points left to fit!')
        if debug:
            pts.set_data(x, y)
            fit.set_data(x, func(x,*bestfitpar))
            pl.draw()
            raw_input('Enter to continue')
        m.migrad()
        bestfitpar = tuple(m.values[k] for k in parnames)
        norm = np.abs(y - func(x, *bestfitpar))
        stdev = norm.std()
        condition =  norm < clip * stdev
        y = y[condition]
        x = x[condition]

    return bestfitpar #,m,x,y


class InterpCubicSpline:
    """Interpolate a cubic spline through a set of points.

    Instantiate the class with two arrays of points: x and
    y = f(x).

    Inputs:

    x                 : array of x values
    y                 : array of y values (y = f(x))
    firstderiv = None : derivative of f(x) at x[0]
    lastderiv  = None : derivative of f(x) at x[-1]

    After initialisation, the instance can be called with an array of
    values xp, and will return the cubic-spline interpolated values
    yp = f(xp).

    The spline can be reset to use a new first and last derivative
    while still using the same initial points by calling the set_d2()
    method.

    If you want to calculate a new spline using a different set of x
    and y values, you'll have to instantiate a new class.

    The spline generation is loosely based on the Numerical recipes
    routines.

    Examples
    --------

    """
    def __init__(self,x,y,firstderiv=None,lastderiv=None):
        x = np.asarray(x).astype('float')
        lx = list(x)
        # check all x values are unique
        if any((lx.count(val) > 1) for val in set(lx)):
            raise Exception('non-unique x values were found!')
        y = np.asarray(y).astype('float')
        cond = x.argsort()                    # sort arrays
        self.x = x[cond]
        self.y = y[cond]
        self.npts = len(x)
        self.set_d2(firstderiv,lastderiv)

    def __call__(self,xp):
        """ Given an array of x values, returns cubic-spline
        interpolated values yp = f(xp) using the derivatives
        calculated in set_d2().
        """
        x = self.x;  y = self.y;  npts = self.npts;  d2 = self.d2

        # make xp into an array
        if not hasattr(xp,'__len__'):  xp = (xp,)
        xp = np.asarray(xp)

        # for each xp value, find the closest x value above and below
        i2 = np.searchsorted(x,xp)

        # account for xp values outside x range
        i2 = np.where(i2 == npts, npts-1, i2)
        i2 = np.where(i2 == 0, 1, i2)
        i1 = i2 - 1

        h = x[i2] - x[i1]
        a = (x[i2] - xp) / h
        b = (xp - x[i1]) / h
        temp = (a**3 - a)*d2[i1] +  (b**3 - b)*d2[i2]
        yp = a * y[i1] + b * y[i2] + temp * h*h/6.

        return yp

    def _tridiag(self,temp,d2):
        x, y, npts = self.x, self.y, self.npts
        for i in range(1,npts-1):
            ratio = (x[i]-x[i-1]) / (x[i+1]-x[i-1])
            denom = ratio * d2[i-1] + 2.       # 2 if x vals equally spaced
            d2[i] = (ratio - 1.) / denom       # -0.5 if x vals equally spaced
            temp[i] = (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1])
            temp[i] = (6.*temp[i]/(x[i+1]-x[i-1]) - ratio * temp[i-1]) / denom
        return temp

    def set_d2(self, firstderiv=None, lastderiv=None, verbose=False):
        """ Calculates the second derivative of a cubic spline
        function y = f(x) for each value in array x. This is called by
        __init__() when a new class instance is created.

        optional inputs:

        firstderiv = None : 1st derivative of f(x) at x[0].  If None,
                             then 2nd derivative is set to 0 ('natural').
        lastderiv  = None : 1st derivative of f(x) at x[-1].  If None,
                             then 2nd derivative is set to 0 ('natural').
        """
        if verbose:  print 'first deriv,last deriv',firstderiv,lastderiv
        x, y, npts = self.x, self.y, self.npts
        d2 = np.empty(npts)
        temp = np.empty(npts-1)

        if firstderiv is None:
            if verbose:  print "Lower boundary condition set to 'natural'"
            d2[0] = 0.
            temp[0] = 0.
        else:
            d2[0] = -0.5
            temp[0] = 3./(x[1]-x[0]) * ((y[1]-y[0])/(x[1]-x[0]) - firstderiv)

        temp = self._tridiag(temp,d2)

        if lastderiv is None:
            if verbose:  print "Upper boundary condition set to 'natural'"
            qn = 0.
            un = 0.
        else:
            qn = 0.5
            un = 3./(x[-1]-x[-2]) * (lastderiv - (y[-1]-y[-2])/(x[-1]-x[-2]))

        d2[-1] = (un - qn*temp[-1]) / (qn*d2[-2] + 1.)
        for i in reversed(range(npts-1)):
            d2[i] = d2[i] * d2[i+1] + temp[i]

        self.d2 = d2

def spline_continuum(wa, fl, er, edges, minfrac=0.01, nsig=3.0,
                      resid_std=1.3, debug=False):
    """ Given a section of spectrum, fit a continuum to it very
    loosely based on the method in Aguirre et al. 2002.

    Parameters
    ----------
    wa               : Wavelengths.
    fl               : Fluxes.
    er               : One sigma errors.
    edges            : Wavelengths giving the chunk edges.
    minfrac = 0.01   : At least this fraction of pixels in a single chunk
                       contributes to the fit.
    nsig = 3.0       : No. of sigma for rejection for clipping.
    resid_std = 1.3  : Maximum residual st. dev. in a given chunk.
    debug = False    : If True, make helpful plots.

    Returns
    -------
    Continuum array, spline points, first derivative at first and last
    spline points

    Examples
    --------
    """

    # Overview:

    # (1) Calculate the median flux value for each wavelength chunk.

    # (2) fit a 1st order spline (i.e. series of straight line
    # segments) through the set of points given by the central
    # wavelength for each chunk and the median flux value for each
    # chunk.

    # (3) Remove any flux values that fall more than nsig*er below
    # the spline.

    # Repeat 1-3 until the continuum converges on a solution (if it
    # doesn't throw hands up in despair! Essential to choose a
    # suitable first guess with small enough chunks).

    if len(edges) < 2:
        raise ValueError('must be at least two bin edges!')

    wa,fl,er = (np.asarray(a) for a in (wa,fl,er))

    if debug:
        ax = pl.gca()
        ax.cla()
        ax.plot(wa,fl)
        ax.plot(wa,er)
        ax.axhline(0, color='0.7')
        good = ~np.isnan(fl) & ~np.isnan(er)
        ymax = 2*sorted(fl[good])[int(len(fl[good])*0.95)]
        ax.set_ylim(-0.1*ymax, ymax)
        ax.set_xlim(min(edges), max(edges))
        ax.set_autoscale_on(0)
        pl.draw()

    npts = len(wa)
    mask = np.ones(npts, bool)
    oldco = np.zeros(npts, float)
    co = np.zeros(npts, float)

    # find indices of chunk edges and central wavelengths of chunks
    indices = wa.searchsorted(edges)
    indices = [(i0,i1) for i0,i1 in zip(indices[:-1],indices[1:])]
    if debug:  print ' indices',indices
    wavc = [0.5*(w1 + w2) for w1,w2 in zip(edges[:-1],edges[1:])]

    # information per chunks
    npts = len(indices)
    mfl = np.zeros(npts, float)     # median fluxes at chunk centres
    goodfit = np.zeros(npts, bool)  # is fit acceptable?
    res_std = np.zeros(npts, float) # residuals standard dev
    res_med = np.zeros(npts, float) # residuals median
    if debug:
        print 'chunk centres',wavc
        cont, = ax.plot(wa,co,'k')
        midpoints, = ax.plot(wavc,mfl,'rx',mew=1.5,ms=8)

    # loop that iterative fits continuum
    while True:
        for i,(j1,j2) in enumerate(indices):
            if goodfit[i]:  continue
            # calculate median flux
            #print i,j1,j2
            w,f,e,m = (item[j1:j2] for item in (wa,fl,er,mask))
            ercond = e > 0
            cond = m & ercond
            chfl = f[cond]
            chflgood = f[ercond]
            if len(chflgood) == 0: continue
            #print len(chfl), len(chflgood)
            if float(len(chfl)) / len(chflgood) < minfrac:
                f_cutoff = percentile(chflgood[ercond], minfrac)
                cond = ercond & (f >= f_cutoff)
            if len(f[cond]) == 0:  continue
            mfl[i] = np.median(f[cond])
        # calculate the spline. add extra points on either end to give
        # a nice slope at the end points.
        extwavc = ([wavc[0]-(wavc[1]-wavc[0])] + wavc +
                   [wavc[-1]+(wavc[-1]-wavc[-2])])
        extmfl = ([mfl[0]-(mfl[1]-mfl[0])] + list(mfl) +
                  [mfl[-1]+(mfl[-1]-mfl[-2])])
        co = np.interp(wa,extwavc,extmfl)
        if debug:
            cont.set_ydata(co)
            midpoints.set_xdata(wavc)
            midpoints.set_ydata(mfl)
            pl.draw()

        # calculate residuals for each chunk
        for i,(j1,j2) in enumerate(indices):
            if goodfit[i]:  continue
            ercond = er[j1:j2] > 0
            cond = ercond & mask[j1:j2]
            chfl = fl[j1:j2][cond]
            chflgood = fl[j1:j2][ercond]
            if len(chflgood) == 0:  continue
            if float(len(chfl)) / len(chflgood) < minfrac:
                f_cutoff = percentile(chflgood[ercond], minfrac)
                cond = ercond & (fl[j1:j2] > f_cutoff)
            #print len(co), len(fl), i1, j1, j2
            residuals = (fl[j1:j2][cond] - co[j1:j2][cond]
                         ) / er[j1:j2][cond]
            res_std[i] = residuals.std()
            if len(residuals) == 0:
                continue
            res_med[i] = np.median(residuals)
            # If residuals have std < 1.0 and mean ~1.0, we might have
            # a reasonable fit.
            if res_std[i] <= resid_std:
                goodfit[i] = True

        if debug:
            print 'median and st. dev. of residuals by region - aiming for 0,1'
            for i,(f0,f1) in  enumerate(zip(res_med, res_std)):
                print '%s %.2f %.2f' % (i,f0,f1)
            raw_input('Enter...')

        # (3) Remove flux values that fall more than N*sigma below the
        # spline fit.
        cond = (co - fl) > nsig * er
        if debug:
            print np.nanmax(np.abs(co - oldco)/co)
        # Finish when the biggest change between the new and old
        # medians is smaller than the number below.
        if np.nanmax(np.abs(co - oldco)/co) < 4e-3:
            break
        oldco = co.copy()
        mask[cond] = False

    # finally fit a cubic spline through the median values to
    # get a smooth continuum.
    d1 = (mfl[1] - mfl[0]) / (wavc[1]-wavc[0])
    d2 = (mfl[-1] - mfl[-2]) / (wavc[-1]-wavc[-2])
    final = InterpCubicSpline(wavc, mfl, firstderiv=d1, lastderiv=d2)

    return final(wa), zip(wavc,mfl), (d1,d2)

def median_continuum(flux, error, numsig=1.5, plot=False):
    """ Estimates the continuum using a median and sigma clipping.

    Given the fluxes and one sigma errors for the section, calculates
    the flux median. Then rejects all flux values less than numsig*sig
    lower than the median, where sig is the median one sigma error
    value of the flux values. Repeat this process, only retaining the
    not-rejected flux values each time, until no flux values are
    rejected. Once this condition is reached, take the current median
    value as the continuum.

    Returns the continuum value.
    """

    if plot:
        pl.plot(flux)
    while True:
        medfl = np.median(flux)
        meder = np.median(error)
        if plot:
            l = pl.axhline(medfl)
        cond = (flux > (medfl - meder*numsig))
        badflux = flux[~cond]
        if len(badflux) == 0:
            return medfl
        flux = flux[cond]
        error = error[cond]

def splice(co0, co1, i, j, forced=None):
    """ Join two overlapping curves smoothly using a cubic spline.

    co0, co1: arrays of shape (N,)
      The two curves to be joined. They must have the same length and
      overlap completely.

    i, j: int
      Roughly speaking, co0 values will be retained below i, and co1 values
      will be retained above j.

    forced: int, optional
      The number of pixels and continuum values between i and j that
      continuum will be forced to pass through.

    Returns:  the new continuum array with shape (N,).
    """
    # go npix to either side of the joining point, and measure slopes
    newco = np.empty_like(co0)

    # derivatives
    d1 = co0[i] - co0[i-1]
    d2 = co1[j] - co1[j-1]

    if forced is not None:
        indices = [i] + list(zip(*forced)[0]) + [j]
        covals = [co0[i]] + list(zip(*forced)[1])+ [co1[j]]
    else:
        indices = [i,j]
        covals = [co0[i],co1[j]]

    spline = InterpCubicSpline(indices, covals,firstderiv=d1, lastderiv=d2)
    newco[:i] = co0[:i].copy()
    newco[i:j] = spline(range(i,j))
    newco[j:] = co1[j:].copy()

    return newco

class InteractiveCoFit(object):
    help_message = """
'a'        : add a new continuum point
'd'        : delete the nearest point
'b'        : add a break in the continuum
'r'        : remove a break in the continuum
'q'        : quit
"""
    def __init__(self, wa, fl, er, contpoints, co=None,
                 nbin=8, redshift=None, atmos=None, fig=None):
        """ Initialise figure, plots and variables.

        Inputs
        ------
        wa:   Wavelengths
        fl:   Fluxes
        er:   One sigma errors
        nbin: int (8)
            Number of pixels to bin arrays in wavelength. Default 8.
        contpoints: list of x,y tuple pairs (None)
            The points through which a cubic spline is passed,
            defining the continuum.
        redshift: float (None)
            Redshift used to plot reference emission lines.
        atmos: list of wavelength pairs (None)
            Regions of atmospheric absorption to plot.

        Updates
        -------
        self.spec:  Dictionary of wa, fl, er.
        self.contpoints:  Points used to define the continuum.
        self.nbin:  The input nbin value.
        self.markers:  Dictionary of matplotlib plotting artists.
        self.connections:  Callback connections.
        self.fig:  The plotting figure instance.
        """
        #setup
        #print co
        self.spec = dict(wa=wa, fl=fl, er=er, co=co)
        self.nbin = nbin
        self.breaks = [wa[0], wa[-1]] # wavelengths of breaks in the continuum
        self.contpoints = list(contpoints)
        self.markers = dict()
        self.fig = (pl.figure() if fig is None else fig)
        # disable any existing key press callbacks
        cids = list(fig.canvas.callbacks.callbacks['key_press_event'])
        for cid in cids:
            fig.canvas.callbacks.disconnect(cid)
        self.connections = []
        self.continuum = None
        self.finished = False
        self.redshift = redshift
        self.atmos = atmos
        self.makefig()
        self.updatefig()
        self.modifypoints()
        pl.show()

    def makefig(self):
        """ Set up the figure and do initial plots.

        Updates
        -------
        self.markers
        """
        wa,fl,er = [self.spec[k][0:-1:self.nbin] for k in 'wa fl er'.split()]
        if self.spec['co'] is not None:
            co = self.spec['co'][0:-1:self.nbin]
        # axis for spectrum & continuum
        a0 = self.fig.add_axes((0.05,0.1,0.9,0.6))
        a0.set_autoscale_on(0)
        # axis for residuals
        a1 = self.fig.add_axes((0.05,0.75,0.9,0.2),sharex=a0)
        a1.set_autoscale_on(0)
        a1.axhline(0,color='g')
        a1.axhline(1,color='g',alpha=0.5)
        a1.axhline(-1,color='g',alpha=0.5)
        a1.axhline(2,color='g',linestyle='dotted',alpha=0.5)
        a1.axhline(-2,color='g',linestyle='dotted',alpha=0.5)
        m0, = a1.plot([0],[0],'ok',ms=3,alpha=0.2)
        a1.set_ylim(-4,4)
        a0.axhline(0, color='0.7')
        if self.spec['co'] is not None:
            a0.plot(wa,co, color='0.7', lw=1, ls='dashed')
        a0.plot(wa,er, color='orange', alpha=0.8)
        a0.plot(wa,fl, 'b', lw=1, linestyle='steps-mid')
        m1, = a0.plot([0],[0],'r',alpha=0.7)
        m2, = a0.plot([0],[0],'o',mfc='None', mew=1, ms=8, mec='r', picker=5,
                      alpha=0.7)
        a0.set_xlim(min(wa), max(wa))
        good = ~np.isnan(er) & ~np.isnan(fl)
        ymin = -np.median(er[good])*5
        ymax = sorted(fl[good])[int(len(fl[good])*0.95)]*2.0
        a0.set_ylim(ymin, ymax)
        if self.redshift is not None:
            junk = plotlines(self.redshift, a0, atmos=self.atmos)
        self.fig.canvas.draw()
        self.markers.update(contpoints=m2, cont=m1, resid=m0)

    def updatefig(self):
        """ Calculates the new continuum, residuals and updates the plots.

        Updates
        -------
        self.markers
        self.continuum
        """
        wa,fl,er = (self.spec[key] for key in 'wa fl er'.split())
        co = np.empty(len(wa))
        co.fill(np.nan)
        for b0,b1 in zip(self.breaks[:-1], self.breaks[1:]):
            cpts = [(x,y) for x,y in self.contpoints if b0 <= x <= b1]
            if len(cpts) == 0:
                continue 
            spline = InterpCubicSpline(*zip(*cpts))
            i,j = wa.searchsorted([b0,b1])
            co[i:j] = spline(wa[i:j])
                
        resid = (fl - co) / er
        self.markers['contpoints'].set_data(zip(*self.contpoints))
        nbin = self.nbin
        self.markers['cont'].set_data(wa[::nbin], co[::nbin])
        self.markers['resid'].set_data(wa[::nbin], resid[::nbin])
        self.continuum = co
        self.fig.canvas.draw()

    def on_keypress(self, event):
        """ Add or remove a continuum point.

        Updates
        -------
        self.contpoints
        """
        if event.key == 'q':
            print """
Disconnected. Continuum is in .continuum, continuum points are in .contpoints.
To resume adjusting points use .modifypoints().
"""
            for item in self.connections:
                self.fig.canvas.mpl_disconnect(item)
            self.finished = True
            return
        if event.inaxes != self.fig.axes[0]:  return
        if event.key == 'a':
            # add a point to contpoints
            x,y = event.xdata,event.ydata
            if x not in zip(*self.contpoints)[0]:
                self.contpoints.append((x,y))
                self.updatefig()
        elif event.key == 'd':
            # remove a point from contpoints
            contx,conty = zip(*self.contpoints)
            sep = np.hypot(event.xdata - contx, event.ydata - conty)
            self.contpoints.remove(self.contpoints[sep.argmin()])
            self.updatefig()
        elif event.key == 'b':
            # Add a break to the continuum.
            self.breaks.append(event.xdata)
            self.breaks.sort()
            self.updatefig()
        elif event.key == 'r':
            # remove a break
            i = indexnear(self.breaks, event.xdata)
            if i not in (0, len(self.breaks)-1):
                self.breaks.remove(self.breaks[i])
            self.updatefig()
        elif event.key == '?':
            print self.help_message

    def modifypoints(self):
        """ Add/remove continuum points."""
        print self.help_message
        id1 = self.fig.canvas.mpl_connect('key_press_event',self.on_keypress)
        self.connections.extend([id1])

def fitqsocont(wa, fl, er, redshift, oldco=None, nbin=1, divmult_forest=None,
               divmult=1, atmos=True, debug=False):

    """ divmult=3 works well for R~40000, S/N~10, z=3 QSO spectrum.

    nbin bins the data for plotting and continuum fitting (obsolete)
    """

    if divmult_forest is None:
        divmult_forest = divmult
    
    # choose initial reference continuum points.  Increase divmult for
    # fewer initial continuum points (generally needed for poorer S/N
    # spectra).

    # atmospheric absorption bands to avoid, from Britt

    zp1 = 1 + redshift
    #reflines = np.array([1025.72, 1215.6701, 1240.14, 1398.0,
    #                     1549.06, 1908,      2800            ])

    # generate the edges of wavelength chunks to send to fitting routine

    # these edges and divisions are generated by trial and error

    # for S/N = 15ish and resolution = 2000ish
    ndiv = np.array([20, 4, 8, 5, 5, 2, 5, 5, 10, 8, 5, 5, 15, 25, 20])
    ndiv[:2] = np.ceil(ndiv[:2] * divmult_forest)
    ndiv[2:] = np.ceil(ndiv[2:] * divmult)    
    edges0 = np.array([800,  1190, 1215, 1263, 1290, 1340, 1370, 1410, 1515,
                       1600, 1800, 1900, 1940, 2240, 3000, 4000])
    edges1 = zp1 * edges0
    temp = [np.linspace(w0,w1,n+1)[:-1] for w0,w1,n in
            zip(edges1[:-1], edges1[1:], ndiv)]
    edges2 = np.concatenate(temp)
    i0,i1,i2 = edges2.searchsorted([wa.min(), 1210*zp1, wa.max()])
    contpoints = []
    #pl.figure(1)
    #pl.clf()
    if i1 - i0 > 3:
        # forest
        co,cp,deriv = spline_continuum(wa, fl, er, edges2[i0:i1], debug=debug,
                                        nsig=2.0, resid_std=1.)
        contpoints.extend(cp)
        # outside forest
        co,cp,deriv = spline_continuum(wa, fl, er, edges2[i1:i2], debug=debug)
    else:
        co,cp,deriv = spline_continuum(wa, fl, er, edges2[i0:i2], debug=debug)
    contpoints.extend(cp)
    fig = pl.figure(2)
    fig.clf()
    wrapper = InteractiveCoFit(wa, fl, er, contpoints, co=oldco, nbin=nbin,
                               redshift=redshift, fig=fig, atmos=atmos)
    while True:
        if wrapper.finished: break
        pl.waitforbuttonpress()

    return wrapper.continuum, wrapper.contpoints

c0 = 0.0224
c1 = -1.56
c2 = 0.0033
c3 = 1073
c4 = 9
c5 = 0.0023
c6 = 1123
c7 = 9
c8 = 0.021
c9 = 1216
c10 = 29
# only use closest region for now, to minimise flux calibration
# problems
powlaw_regions = (1330, 1360), #(1440, 1490), (1650, 1680)

def lyafcont(wa, fl):
    """ Using the model of Bernardi et al, guess the continuum for the
    lya forest. Only works for flux-calibrated spectra. Continuum is
    only valid over the lya forest region!

    Inputs
    ------
    wa, fl: arrays, shape (N,)

    Returns
    -------
    continuum: array, shape (N,)
    """
    # work out c0 from (hopefully) clean regions outside forest
    guesses = []
    for w1,w2 in powlaw_regions:
        cond = (w1 <= wa) & (wa <= w2)
        if not np.any(cond):
            print 'no good regions found!'
            continue
        print 'mean wa', np.mean(wa[cond])
        print 'median fl', np.median(fl[cond])
        c0guess = np.median(fl[cond]) / (c0 * (np.mean(wa[cond]) / wlya) ** c1)
        guesses.append(c0guess)

    mult =  np.median(guesses) if len(guesses) > 0 else 1.

    print 'multiplier=',mult
    return mult * ( c0 * (wa / wlya) ** c1 +
                    c2 * np.exp(-0.5 * ((wa - c3) / c4)**2) +
                    c5 * np.exp(-0.5 * ((wa - c6) / c7)**2) +
                    c8 * np.exp(-0.5 * ((wa - c9) / c10)**2) )
    
def poly_scale_spec(s, sref, mask=None, order=9, clip=None, debug=True):
    """ Remove the large-scale variations between two spectra by
    dividing them by one another and then fitting a low-order
    polynomial.

    s           : Spectrum object to be scaled to match reference
    sref        : Reference spectrum object
    order = 9   : Order of polynomial to fit
    mask = None : Input to polyfitr, masks wavelength regions.
    clip = None : Number of sigma to clip after each iteration.

    Returns
    -------
    out: array
        sr * out matches s
    """

    # Assume reference spectrum covers a larger wavelength range.
    s1 = sref.rebin(wstart=s.wa[0], dw=s.dw, dv=s.dv, npts=len(s.wa))

    c0 = (s1.er > 0) & (s.er > 0)
    c1 = ~np.isnan(s1.fl) & ~np.isnan(s.fl) & (s1.fl!=0) & (s.fl!=0)
    good = c0 & c1
    scale = s1.fl[good] / s.fl[good]
    wa = s.wa[good]
    pl.plot(wa, scale, '.')
    if clip is not None:
        coeff,x,y = polyfitr(wa, scale, order=order, mask=mask,
                             debug=debug, clip=clip)
    else:
        coeff,x,y = polyfitr(wa, scale, order=order, mask=mask, debug=debug)

    pl.plot(s.wa, np.polyval(coeff,s.wa))

    return np.polyval(coeff,s.wa)
