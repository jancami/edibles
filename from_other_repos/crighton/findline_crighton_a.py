""" Module that identifies absorption features in a spectrum.

Generally you just want to use findline.findlines().  Other routines
included are all necessary for findlines().

plotfeatures takes the output of findlines() and plots the detected
features.
"""
# Given a spectrum with an error array, fitted continuum, and
# resolution, search for absorption lines.

from astro.spec import find_wa_edges
from astro.plot import axvlines
from astro.utilities import find_edges_true_regions
from math import sqrt
import numpy as np
import pdb

def find_ew_per_resel(wa, fl, er, co, resolution, weighted_aperture=True):
    """ Finds the ew over a resolution element for each pixel.

    Input
    -----
    wa,fl,er,co : arrays of floats
        The wavelength, flux, one sigma error and continuum arrays for
        the spectrum. All must be the same length.

    resolution : float or array of floats
        Resolution of the spectrum = wav / dwav = c / dvel.  It can be
        an array with the same length as wa.

    weighted_aperture : bool (True)
        If True, weight the equivalenet width per resolution element
        by the instrumental spread function, assumed to be a Gaussian
        with FWHM determined using the resolution.
    
    Returns
    -------
    sp : Record array
        Same length as the input arrays. Fields of `sp` give the
        following information:

        wa      :  wavelength
        fl      :  flux
        er      :  one sigma error
        co      :  continuum
        dw      :  wavelength width per pixel
        ew      :  equivalent width per pixel (positive is absorption)
        ewer    :  error in the equivalent width per pixel
        ewres   :  equivalent over a resolument element per pixel
        ewreser :  error in the equivalent width over one resolution element
        good    :  Is the pixel good? (error > 0, flux, continuum defined)
        fwhmpix :  The fwhm of the instrumental spread function per pixel in
                   units of pixels
    """

    good = (er > 0) & ~np.isnan(fl) & ~np.isnan(co)
    if not np.any(good):
        raise ValueError('No good pixels!')
    edges = find_wa_edges(wa)            # pixel edges
    dw = edges[1:] - edges[:-1]          # width per pixel
    nfl = fl / co                        # normalised flux
    ner = er / co                        # normalised 1 sigma error
    dec = 1 - nfl                        # flux decrement
    decer = ner                          # error in decrement
    # equivalent width per pixel
    ew = dw * dec
    ewer = dw * decer                   # equivalent width error
    fwhmwa = wa / resolution            # resolution fwhm in wavelength units
    fwhmpix = fwhmwa / dw               # resolution fwhm in pixels
    # number of pixels per FWHM res. element (round up)
    npixels = np.ceil(fwhmpix)

    sp = np.rec.fromarrays([wa, dw, ew, ewer, good],
                           names='wa,dw,ew,ewer,good')
    npts = len(wa)
    ewres = []
    ewres_err = []
    for i,(npix,fwhm) in enumerate(zip(npixels, fwhmwa)):
        i0 = max(0, i - npix)
        i1 = min(npts, i + npix + 1)
        s = sp[ int(i0):int(i1) ]
        if weighted_aperture:
            # find the instrumental spread function
            sigma = fwhm / 2.35
            temp = (s.wa - wa[i]) / sigma
            prof0 = np.exp(-temp*temp)
            prof1 = prof0 / prof0.sum()
            #raw_input(str(prof1))
        if not np.any(s.good):
            ewres.append(np.nan)
            ewres_err.append(np.nan)
            continue
        if weighted_aperture:
            weighted_ew = s.ew * prof1
            weighted_ew_err = s.ewer * prof1
            # Normalisation suggested by Schneider, but
            # I can't see why.  Read Horne 86?
            #norm = (prof1 ** 2)[s.good].sum()
            #
            # This normalisation allows direct comparison to the
            # non-weighted version
            norm = 1. / len(s[s.good])
            ewres.append( weighted_ew[s.good].sum() / norm )
            ewres_err.append( sqrt((weighted_ew_err[s.good]**2).sum()) / norm )
        else:
            s = s[s.good]
            ewres.append(s.ew.sum())
            ewres_err.append(sqrt((s.ewer**2).sum()))

    # spectrum
    sp = np.rec.fromarrays(
        [wa, fl, er, co, dw, ew, ewer, ewres, ewres_err, good, fwhmpix],
        names='wa,fl,er,co,dw,ew,ewer,ewres,ewreser,good,fwhmpix')

    return sp

def _remove_close(indices, signif, minsep, nmin=1, joinfactor=1.2):
    """ Remove a line from indices if it is within fwhmpix pixels of
    another, more significant line, or it is not differentiated
    sufficiently from a nearby more significant component.

    indices must be sorted.

    joinfactor (>1) determines how easily a feature is broken up into
    components. Increasing it means fewer components in blends,
    decreasing means more components.
    """
    newindices = [indices[0]]
    for i1,i2 in zip(indices[1:-1], indices[2:]):
        i0 = newindices[-1]
        sig1 = signif[i1]
        # if the component we're looking at, i1, is less than 1
        # resolution elemnt away from a more significant component,
        # then drop i1.
        if i1 - i0 < minsep[i1]:
            if sig1 < signif[i0]:
                continue
        elif i2 - i1 < minsep[i1]:
            if sig1 < signif[i2]:
                continue
        minsig0 = signif[i0:i1 + 1].min()
        # if (1) i1 is one of multiple components inside a region
        # above significance of nmin, (2) a more significant component
        # than i1 exists inside the same region and (3) i1 is not a
        # well-defined component, then drop it.
        if minsig0 > nmin and sig1 < signif[i0] and \
               sig1 / minsig0 < joinfactor:
            continue
        minsig1 = signif[i1:i2+1].min()
        if minsig1 > nmin and sig1 < signif[i2] and \
               sig1 / minsig1 < joinfactor:
            continue
        newindices.append(i1)
    if len(indices) > 2:
        newindices.append(indices[-1])
    elif len(indices) == 2:
        # special case.
        i0, i1 = indices
        sig0, sig1 = signif[i0], signif[i1]
        # if the component we're looking at, i1, is less than 1
        # resolution elemnt away from a more significant component,
        # then drop i1.
        if i1 - i0 < minsep[i1]:
            if sig1 < sig0:
                newindices = [i0]
            else:
                newindices = [i1]

        minsig = signif[i0:i1 + 1].min()
        # if (1) i1 is one of multiple components inside a region
        # above significance of nmin, (2) a more significant component
        # than i1 exists inside the same region and (3) i1 is not a
        # well-defined component, then drop it.
        if minsig > nmin:
            if sig1 < sig0 and sig1 / minsig < joinfactor:
                newindices = [i0]
            elif sig0 / minsig < joinfactor:
                newindices = [i1]
        else:
            newindices.append(i1)

    return newindices

def find_feature_edges(centres, above_nmin, signif):
    """ Given the line centres and significance array, add the line
    edges, and whether it is part of a blend.

    Returns indices of left, centre, right edges and whether is it a
    blend (1) or not(0).

    TODO: change so that this creates one list of feature edges
    (ignoring blends), and another list of component positions.
    
    """
    # now find feature edges
    centres = np.asarray(centres)
    ind0,ind1 = find_edges_true_regions(above_nmin)
    cond =  ind1 - ind0 >= 2
    ind0 = ind0[cond]
    ind1 = ind1[cond]
    feature = []
    for i0,i1 in zip(ind0,ind1):
        cond = (i0 < centres) & (centres < i1)
        if not np.any(cond):
            continue
        cen = centres[cond]
        if len(cen) == 1:
            # not blended
            feature.append([i0-1, cen[0], i1+1])
            continue
        #pl.vlines(cen, 0, 10, colors='y')
        right = cen[0] + signif[cen[0]:cen[1]+1].argmin()
        feature.append([i0-1, cen[0], right])
        #pl.vlines([i0-1,right[-1]], 0, 10)
        #pl.vlines([cen[0]], 0, 10, colors='r')
        #raw_input()
        for j in range(1, len(cen)-1):
            left = right
            c0,c1 = cen[j:j+2]
            right = c0 + signif[c0:c1+1].argmin()
            feature.append([left, c0, right])
            #pl.vlines([left[-1],right[-1]], 0, 10)
            #pl.vlines([c0], 0, 10, colors='r')
            #raw_input()
        left = right
        feature.append([left, cen[-1], i1+1])
        #pl.vlines([left[-1],right[-1]], 0, 10)
        #pl.vlines([cen[-1]], 0, 10, colors='r')
        #raw_input()

    # make sure the start and end indices are not negative or longer
    # than the array
    feature[0][0] = max(0, feature[0][0])
    feature[-1][2] = min(len(signif)-1, feature[-1][2])

    return zip(*feature)

def find_feature_indices(ew, ew_err, fwhmpix, ndetect=5, nmin=1,
                         joinfactor=1.2):
    """ Find the start, peak significance and end pixel indices of
    absorption and emission features.

    Note ew must be positive for absorption."""

    good = ~np.isnan(ew) & (ew_err > 0)

    # note should really be using the interpolated error here, without
    # it features will be more significant than they should be.
    signif = ew / ew_err

    # gradually decrease the detection level from the significance of
    # the maximum feature to the minimum level, adding features as the
    # rise above the detection level.
    detectlevels = np.sort(signif[good])

    lines = set()
    for level in reversed(detectlevels):
        if level < ndetect: break
        condition = (signif >= level) & good
        # find the edges of contiguous regions
        ind0,ind1 = find_edges_true_regions(condition)
        # add to lines where we have the tip of a feature
        lines.update(ind0[ind0-ind1 == 0])

    if not lines:
        return None
    lines = sorted(lines)
    # remove lines that are within one resolution element of another
    # more significant line
    lines = _remove_close(lines, signif, fwhmpix, joinfactor=joinfactor)
    left,indices,right = find_feature_edges(lines, (signif >= nmin) & good, signif)
    return zip(left,lines,right)

def find_feature_fwhm(wa, fl, co):
    """ Rough estimate of full width at half maximum for a feature."""
    nfl = fl / co
    # find half the maximum flux decrement
    halfmin = 0.5 * max(1-nfl)

    # move in from each edge until we reach half the maximum
    imax = -1
    while 1 - nfl[imax] < halfmin:
        imax -= 1

    imin = 0
    while 1 - nfl[imin] < halfmin:
        imin += 1

    if imax == imin: imax += 1

    return wa[imax] - wa[imin]

def find_feature_centre(wa, ew):
    """ Estimates the wavelength centre of a feature and its
    error. From Churchill QSO textbook."""
    ewtot = sum(ew)
    wa_ewtot = sum(wa * ew)
    try:
        centre = wa_ewtot / ewtot
    except:
        print wa, ew
        raise
    return centre

def findlines(wa, fl, er, co, resolution, ndetect=5,
              cull_atmos=False, joinfactor=1.2):
    """ Finds features in a spectrum, inspired by the Schneider et
    al. method (1993, ApJS, 87:45-62), although it doesn't use an
    interpolated error.

    resolution : Resolution of the spectrum = wav / dwav = c / dvel
                 Can be an array with the same length as wa.

    The instrumental spread function is assumed to be a gaussian
    function with fwhm, dwav, given by the resolution according to the
    above relation.

    ndetect is the significance level required before something is
    considered a feature.

    if cull_atmos is True, lines that are inside regions of known
    atmospheric absorption are removed.

    Returns the detected features with estimates for their wavelength,
    width, equivalent width
    """
    assert resolution > 100, 'resolution seems too low, %s' % resolution
    # equivalent width for one resolution element per pixel
    #sp1 = find_ew_per_resel(wa, dw, ew0, ew0er, resolution)
    sp = find_ew_per_resel(wa, fl,er,co, resolution, weighted_aperture=1)

    indices = find_feature_indices(sp.ewres, sp.ewreser, sp.fwhmpix,
                                   ndetect=ndetect, joinfactor=joinfactor)
    if indices is None:
        return None, sp

    # estimate line centres, widths, equivalent widths, and
    # significance.
    ewidth = []
    ewidth_sig = []
    centres = []
    fwhm = []
   
    for i0,i,i1 in indices:
        s = sp[i0:i1+1]
        cond = s.good
        if not np.any(cond):
            ewidth.append(np.nan)
            ewidth_sig.append(np.nan)
            centres.append(np.nan)
            fwhm.append(np.nan)
            continue
        s = s[cond]
        ewidth.append( sum(s.ew) )
        ewidth_sig.append( sqrt(sum(s.ewer**2)) )
        centres.append(find_feature_centre(s.wa, s.ew))
        fwhm.append( find_feature_fwhm(s.wa, s.fl, s.co) )

    idens = ['?          '] * len(indices)
    redshifts = np.zeros(len(indices))
    num = range(1, len(indices)+1)
    ind = zip(*indices)
    ind0 = np.array(ind[0])
    ind2 = np.array(ind[2])
    features = np.rec.fromarrays(
        [num, ewidth, ewidth_sig, fwhm, redshifts, idens, sp.wa[ind0],
         centres, sp.wa[ind2], ind[0], ind[1], ind[2]],
        names='num,ew,ewer,fwhm,z,iden,wa0,wac,wa1,i0,ic,i1')

    return features, sp

def plotfeatures(lines, sp, redshift=None, fig=None,show=True):
    """ plots the spectrum and features found using findlines.  Takes
    output of findlines as input.  Requires matplotlib."""
    try:
        import matplotlib.pyplot as pl
    except ImportError:
        print('**** Pylab not found; plotfeatures() unavailable. ****')
        raise
    else:
        from astro.utilities import indexnear
        from astro.spec import plotlines

    wa,fl,er,co = sp.wa, sp.fl, sp.er, sp.co
    if fig is None:
        fig = pl.gcf()
    fig.clf()
    a = fig.add_subplot(211)
    a.plot(wa,fl,'b', alpha=0.7, ls='steps-mid')
    a.plot(wa,er,'orange', alpha=0.7)
    a.plot(wa,co,'r', alpha=0.7)
    a.axhline(0, color='gray')
    sortfl = np.sort(fl[sp.good])
    if lines is not None:
        ref = sortfl[int(len(sortfl)*0.95)] -  sortfl[int(len(sortfl)*0.05)]
        i = np.concatenate( [lines.i0, lines.i1] )
        #ref = np.median(fl[sp.good])
        axvlines(wa[i], ax=a, colors='k', alpha=0.5)
        temp = co[np.array([indexnear(wa,wav) for wav in lines.wac])]
        ymin = temp + ref*0.2
        ymax = ymin + ref*0.4
        a.vlines(lines.wac, ymin, ymax, colors='r')
    a.set_ylim(0, 1.5)#sortfl[int(0.95*len(sortfl))]*1.2)
    if redshift is not None:
        plotlines(redshift, a, atmos=True)
    a1 = fig.add_subplot(212, sharex=a)
    a1.axhline(0, color='gray')
    a1.plot(wa, sp.ewres, 'k', lw=2, alpha=0.7)
    a1.plot(wa, sp.ewreser, 'g', lw=2, alpha=0.7)
    cond = sp.good & (sp.ewres > 0)
    sortew = np.sort(sp.ew[cond])
    ymax = 1.0 #np.sort(sortew[int(0.99*len(sp[cond]))]*5)
    a1.set_ylim(-0.2, ymax)
    a1.set_xlabel('Wavelength')
    a1.set_ylabel('EW per pixel over res. FWHM')
    a1.set_xlim(wa.min(), wa.max())
    if show:
        pl.show()
    return fig
