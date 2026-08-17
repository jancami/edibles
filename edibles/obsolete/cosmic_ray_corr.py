import numpy as np
import matplotlib.pylab as plt
import scipy


# %%%%%%%%%% Gaussian function %%%%%%%%%%%%%%%%%
def gauss(xw, H, A, x0, sigma):
    return H + A * np.exp(-(xw - x0) ** 2 / (2 * sigma ** 2))


# %%%%%%%%%% fitting function %%%%%%%%%%%%%%%%%
def gaussian_fit(x, y, mean=None, sigma=None, bounds=False,
                       linefwhm=None, miny=None, maxy=None):
    if mean == None:
        mean = sum(x * y) / sum(y)
    else: mean = mean
    if sigma == None:
        sigma = np.sqrt(sum(y * (x - mean) ** 2) / sum(y))
    else: sigma = sigma
    if miny is None: miny = min(y)
    if maxy is None: maxy = max(y)
    try:
        if bounds == False:
            popt, pcov = scipy.optimize.curve_fit(gauss, x, y, p0=[miny, maxy, mean, sigma],
            bounds=((-np.inf, -np.inf, mean-0.1, -np.inf), (np.inf, np.inf, mean+0.1, np.inf)))
        else:
            popt, pcov = scipy.optimize.curve_fit(gauss, x, y, p0=[miny, maxy, mean, linefwhm],
            bounds=((-np.inf, -np.inf, mean-2*sigma, 0), (np.inf, np.inf, mean+2*sigma, 1.5*linefwhm)) )
    except:
        popt = [1, 1, 1, 1]
    return popt


# %%%%%%%%%% cosmic ray correction %%%%%%%%%%%%%%%%%
def remove_spike(x, y, method='sigmaclipping', savename=None highclip=3, res=None):

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #   spike removing using median
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if method == 'gauss':

        # make a copy from original data
        xclean = x.copy()
        yclean = y.copy()

        # measure resolving power
        if res is None:
            line_width = np.diff(x)[0]
        else:
            dlambda = np.mean(x)/res
            line_width = dlambda / (2 * np.sqrt(2. * np.log(2.)))


        # apply median filtering for computing the first sigma
        ym = scipy.ndimage.median_filter(y, size=5)
        sigma = np.std(y-ym)
        sig0 = sigma



        # points more than 1-sigma away
        xspikes = []
        yspikes = []
        yyt = []
        idx = [xx for xx, e in enumerate(y-ym) if (e > sigma) or (e < -sigma)]
        for looptt in range(len(idx)):
            xspikes.append(x[idx[looptt]])
            yspikes.append(y[idx[looptt]])
            yyt.append(y[idx[looptt]]-ym[idx[looptt]])



        # %%%%%%%%%%%%%% Iteration %%%%%%%%%%%%%%%
        idx_line = []
        idx_spike_temp = []
        for iteration in range(10):

            # fit Gaussian to spikes
            for loop in range(len(xspikes)):
                try:
                    l0 = xspikes[loop]
                    idxt = [xx for xx, ee in enumerate(x) if (ee >= l0-1) & (ee <= l0+1)]
                    xx = []
                    yy = []
                    for looplis in range(len(idxt)):
                        xx.append(x[idxt[looplis]])
                        yy.append(y[idxt[looplis]])

                    H1, A1, x01, sigma1 = gaussian_fit(xx, yy, mean=l0, sigma=line_width)
                    yg = gauss(xx, *gaussian_fit(xx, yy, mean=x01, sigma=sigma1))


                    # --------------
                    # select spikes
                    # --------------
                    if (np.abs(sigma1) <= line_width) & (A1 >= 0):

                        idxt = [xx for xx, ee in enumerate(x) if (ee >= l0-3*np.abs(sigma1)) &
                            (ee <= l0+3*np.abs(sigma1))]
                        idx_spike_temp.extend(idxt)

                    # -------------
                    # select lines
                    # -------------
                    else:
                        # to avoid selecting the bad fitting with wrong large widths as actual lines
                        if (np.abs(sigma1) < 10*line_width) & (np.abs(sigma1) > line_width):
                            idxt = [xx for xx, ee in enumerate(x) if (ee >= l0-5*np.abs(sigma1)) &
                                (ee <= l0+5*np.abs(sigma1))]
                            if len(idxt) > 0:
                                idx_line.extend(idxt)
                except:
                    pass

            # remove lines and spikes from spectrum
            idx_line = [*set(idx_line)]
            idx_nonline = list(set(np.arange(len(x))).difference(set(idx_line)))
            idx_spike_temp = [*set(idx_spike_temp)]
            idx_noline_nospike = list(set(idx_nonline).difference(set(idx_spike_temp)))


            xx = []; yy = []
            for looplis in range(len(idx_noline_nospike)):
                xx.append(x[idx_noline_nospike[looplis]])
                yy.append(y[idx_noline_nospike[looplis]])

            # measure UPDATED sigma from CLEANED spectrum
            ym = scipy.ndimage.median_filter(yy, size=5)
            sigma = np.std(yy-ym)


            if np.abs(sig0 - sigma) > 1.e-4:
                # apply sigma-clipping based on new sigma on actual spectrum.
                ym = scipy.ndimage.median_filter(y, size=5)
                idx = [xx for xx, ee in enumerate(y-ym) if ee >= highclip*sigma]
                # Remove lines from selected spikes
                idx_spike = list(set(idx).difference(set(idx_line)))
                idx_spike.sort()

                xspikes = []
                yspikes = []
                for looptt in range(len(idx_spike)):
                    xspikes.append(x[idx_spike[looptt]])
                    yspikes.append(y[idx_spike[looptt]])
                sig0 = sigma

            else:
                ym = scipy.ndimage.median_filter(y, size=15)

                # replace spikes
                for loopid in range(len(idx_spike_temp)):
                    idxt = idx_spike_temp[loopid]
                    # meausre sigma from 2 pixels on both sides
                    sigma = np.abs(np.std(ym[idxt-3:idxt+2]))
                    noise = np.random.normal(0, sigma, 1)
                    yclean[idxt] = ym[idxt] + noise[0]

                break


        plt.plot(x, y, 'gray')
        plt.plot(xclean, yclean, 'k')
        plt.xlabel('Wavelength ($\AA$)')
        plt.ylabel('Flux')
        plt.title('final result')
        plt.show()
    return xclean, yclean



# example
x, y = np.loadtxt('HD63804_3302.ascii', unpack=True)
xclean, yclean = remove_spike(x, y, highclip=3, res=71500, method='gauss')
