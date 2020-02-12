import numpy as np
from astropy.io import fits
import astropy.constants as cst
import matplotlib.pyplot as plt
import edibles_Working.edibles_function as eF
datadir = '/media/EDIBLES/EDIBLES_DATADIR'
#filename = '/HD169454/RED_860/HD169454_w860_redl_20160808_O12.fits'

class EdiblesSpectrum:
    # This object will contain a spectrum from EDIBLES,
    # and a set of methods to operate on the data.

    def loadSpectrum(self):
        # Assume the file is a DR3 product here.
        hdu = fits.open(self.filename)
        self.header = hdu[0].header
        self.target = self.header["OBJECT"]
        self.date = self.header["DATE-OBS"]
        self.flux = hdu[0].data
        self.flux_units = "arbitrary"
        crval1 = self.header["CRVAL1"]
        cdelt1 = self.header["CDELT1"]
        nwave = len(self.flux)
        grid = np.arange(0, nwave, 1)
        self.wave = (grid) * cdelt1 + crval1
        self.wave_units = "AA"
        self.reference_frame = "geocentric"
        self.v_bary = self.header["HIERARCH ESO QC VRAD BARYCOR"]
        #self.bary_wave = self.wave + (self.v_bary/cst.c.to('km/s').value)*self.wave

        spec_name = self.filename.split('/')[-1]
        self.spec_name = spec_name.split('.')[0]
        self.__resetLabels__()

        return None

    def __init__(self, filename):
        """
        Filename is relative to the DR3 directory
        """
        self.filename = datadir + filename
        self.loadSpectrum()

    def statusReport(self):
        print("Spectrum name: {spname}".format(spname=self.spec_name))

        print("Current reference frame: {ref_frame}".format(ref_frame=self.reference_frame))

        if self.wave_units=="AA":
            print("Wavelength is in AA, between {wav_min:.2f} and {wav_max:.2f}"
                  .format(wav_min=np.min(self.wave_working),wav_max=np.max(self.wave_working)))
        else:
            print("Wavelength in km/s")
            print("centers at {wav_center:.2f}, and spans over between "
                  "{wav_min:.2f}km/s and {wav_max:.2f}km/s"
                  .format(wav_center=self.velocity_center,
                          wav_min=np.min(self.wave_working), wav_max=np.max(self.wave_working)))

        if self.normalized:
            print("The spectrum is normalized")
        else:
            print("The spectrum is NOT normalized")

        if self.bary_corr == True:
            print("Shifted for {bary_shift:.2f}km/s for barycentric correction "
                  "plus {other_shift:.2f}km/s by user"
                  .format(bary_shift=self.v_bary, other_shift=self.total_offset))
        else:
            if self.total_offset != 0:
                print("Not been barycentric corrected, "
                      "but {other_shift:.2f}km/s by user"
                      .format(other_shift=self.total_offset))
            else:
                print("The spectrum is in geocentric frame and never shifted")

        return None

    def showSpectrum(self):
        fig, ax = plt.subplots(1,1)
        ax.plot(self.wave_working,self.flux_working, color='k')
        self.__XYLabel__(ax)
        self.__auxiliaryPlot__(ax)
        plt.show()
        return None

    def __XYLabel__(self, ax=None):
        # This is all about axis and labels
        if ax is None: ax = plt.gca()
        if self.velocity_center is None:
            ax.set_xlabel('AA')
        else:
            ax.set_xlabel('km/s')

        if self.optical_depth:
            ax.set_ylabel('Optical Depth')
            ax.invert_yaxis()
        else:
            if self.normalized:
                ax.set_ylabel('Normalized Flux')

            else:
                ax.set_ylabel('Flux')

        return None

    def __auxiliaryPlot__(self, ax=None):
        # trivial matters, continuum at unity and model for now
        if ax is None: ax = plt.gca()

        if self.normalized:
            ax.plot(self.wave_working, np.ones_like(self.wave_working), linestyle='--', color='orange')

        if self.model is not None:
            ax.plot(self.wave_working, self.model(self.wave_working), color='red')

        return None

    def resetSpectrum(self):
        message = 'All handling of spectrum will be abort, continue?[Y/N]'
        go_flag = eF.go_message(message)

        if go_flag:
            self.__resetLabels__
        else:
            pass

        return None

    def __resetLabels__(self):
        self.wave_working = np.copy(self.wave)
        self.flux_working = np.copy(self.flux)
        self.bary_corr = False
        self.normalized = False
        self.velocity_center = None
        self.total_offset = 0
        self.optical_depth = False
        self.reference_frame = "geocentric"
        self.wave_units = "AA"
        self.SNR = None
        self.continuum = None
        self.model = None
        self.model_fit = False
        self.masked = np.zeros_like(self.wave_working)
        return None

    def baryCorrection(self):
        assert self.bary_corr is False, "Already in barycentric reference frame"
        if self.velocity_center is None:
            self.wave_working = [self.wave_working[i] * (1 + self.v_bary / cst.c.to("km/s").value)
                                 for i in range(len(self.wave_working))]
        else:
            self.wave_working = [self.wave_working[i] + self.v_bary
                                 for i in range(len(self.wave_working))]
        self.bary_corr = True
        self.reference_frame = "barycentric"

        return (self.wave_working, self.flux_working)

    def cutSpectrum(self, xmin=None, xmax=None, center=None, span=None):
        xmin_lst = [np.min(self.wave_working)]
        xmax_lst = [np.max(self.wave_working)]
        if xmin is not None: xmin_lst.append(xmin)
        if xmax is not None: xmax_lst.append(xmax)
        if center is not None and span is not None:
            xmin_lst.append(center - 0.5 * span)
            xmax_lst.append(center + 0.5 * span)

        xmin = np.max(xmin_lst)
        xmax = np.min(xmax_lst)
        if not (np.min(self.wave_working) <= xmin <= np.max(self.wave_working)):
            print("Warning, x_min outside working region, rest to current minimal")
            xmin = np.min(self.wave_working)
        if not (np.min(self.wave_working) <= xmax <= np.max(self.wave_working)):
            print("Warning, x_max outside working region, rest to current maximal")
            xmax = np.max(self.wave_working)
        assert xmin < xmax, "xmin must be less than xmax"

        idx = (self.wave_working > xmin) * (self.wave_working < xmax)
        self.wave_working = self.wave_working[idx]
        self.flux_working = self.flux_working[idx]
        self.masked = self.masked[idx]
        return (self.wave_working, self.flux_working)

    def cutSpectrum_visual(self):
        # an interactive backends will be needed
        import matplotlib
        matplotlib.use('Qt5Agg', warn=False, force=True)
        import matplotlib.pyplot as tmp_plt

        fig1, ax = tmp_plt.subplots(1, 1)
        ax.plot(self.wave_working, self.flux_working)
        ax.set_xlim(1.1 * np.min(self.wave_working) - 0.1 * np.max(self.wave_working),
                    1.1 * np.max(self.wave_working) - 0.1 * np.min(self.wave_working))
        self.__XYLabel__(ax)

        print("Select the boundaries of working region")
        points = tmp_plt.ginput(2, mouse_add=1, mouse_pop=3, mouse_stop=2)

        # now let's switch back to the normal backend
        tmp_plt.close(fig1)
        matplotlib.use("module://backend_interagg")
        import matplotlib.pyplot as plt

        # cut the data via method "cutSpectrum"
        (a,b) =self.cutSpectrum(xmin=points[0][0], xmax=points[1][0])
        return a,b

    def addMask(self, n=None):
        self.resetMask()

        if n is None:
            n = input('How many pieces of spectrum should be masked?')
            n = np.int(n)

        # an interactive backends will be needed
        import matplotlib
        matplotlib.use('Qt5Agg', warn=False, force=True)
        import matplotlib.pyplot as tmp_plt

        fig1, ax = tmp_plt.subplots(1, 1)
        ax.plot(self.wave_working, self.flux_working)
        ax.set_xlim(1.1*np.min(self.wave_working)-0.1*np.max(self.wave_working),
                    1.1*np.max(self.wave_working)-0.1*np.min(self.wave_working))
        self.__XYLabel__(ax)

        print("Select the boundaries of the regions to mask")
        print("left to add; right to pop; middle to stop")
        points = tmp_plt.ginput(2*n, mouse_add=1, mouse_pop=3, mouse_stop=2)

        # now let's switch back to the normal backend
        tmp_plt.close(fig1)
        matplotlib.use("module://backend_interagg")
        import matplotlib.pyplot as plt

        n_pairs = int(np.floor(len(points)/2))
        boundary_left = np.array([])
        boundary_right = np.array([])
        for i in range(n_pairs):
            boundary_left = np.append(boundary_left, points[i*2][0])
            boundary_right = np.append(boundary_right, points[i*2+1][0])


        idx_mask = [i for i in range(len(self.wave_working))
                    if not eF.not_within_boundaries(self.wave_working[i], boundary_left, boundary_right)]

        fig2, ax = plt.subplots(1,1)
        ax.plot(self.wave_working, self.flux_working)
        idx_mask2 = eF.continous_idx(idx_mask)
        for i in range(len(idx_mask2)):
            ax.plot(self.wave_working[idx_mask2[i]], self.flux_working[idx_mask2[i]]
                    ,marker='x', markersize=3,color='red')
        self.__XYLabel__(ax)
        plt.show()

        message = 'Remove these points from continuum fitting?[Y/N]'

        go_flag = eF.go_message(message)
        if go_flag:
            self.masked = np.zeros_like(self.wave_working)
            self.masked[idx_mask] = 1
            return idx_mask
        else:
            return None

    def resetMask(self):
        self.masked = np.zeros_like(self.wave_working)

    def fitContinuum(self, mode='spline', n=3,
                     lsigma=1, usigma=2, iterates=30, min_sigma = 0.1,
                     silence=False, apply_mask=True):
        data_tuple = (self.wave_working, self.flux_working)
        if apply_mask: mask = self.masked
        else: mask = np.zeros_like(self.wave_working)
        (cont, idx) = eF.iterate_continuum(data_tuple,mode=mode, n=n,
                                           lsigma=lsigma, usigma=usigma,
                                           iterates=iterates, min_sigma=min_sigma,
                                           mask = mask)

        # making plot
        go_flag = True
        if not silence:
            fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, sharex=True)
            # origional
            ax1.plot(self.wave_working, self.flux_working, color='k')
            ax1.plot(self.wave_working, cont(self.wave_working), color='blue')
            ax1.set_ylabel("Original")

            # normalized
            ax2.plot(self.wave_working, self.flux_working / cont(self.wave_working), color='k')
            ax2.plot(self.wave_working, np.ones_like(self.wave_working), linestyle='--', color='blue')
            ax2.set_ylabel("Normalized")

            # highlights
            idx2 = eF.continous_idx(idx)
            for i in range(len(idx2)):
                idx_plot = idx2[i]
                ax1.plot(self.wave_working[idx_plot], self.flux_working[idx_plot],
                         linestyle='', marker='o', markersize=3, color='red')
                ax2.plot(self.wave_working[idx_plot], self.flux_working[idx_plot] / cont(self.wave_working[idx_plot]),
                         color='red')

            # residual
            std_res = self.flux_working - cont(self.wave_working)
            ax3.plot(self.wave_working, std_res, color='k')
            std_res = np.std(std_res[idx])
            ax3.plot(self.wave_working, np.ones_like(self.wave_working) * std_res, linestyle='--', color='blue')
            ax3.plot(self.wave_working, np.ones_like(self.wave_working) * (-std_res), linestyle='--', color='blue')
            ax3.set_ylim([-std_res * 2, std_res * 2])
            ax3.set_ylabel("Residual")

            ax1.grid()
            ax2.grid()
            ax3.grid()
            plt.show()

            message = 'Keep this continuum?[Y/N]'
            go_flag = eF.go_message(message)

        if go_flag:
            self.flux_working = self.flux_working / cont(self.wave_working)
            self.normalized = True
            #self.SNR = 1 / np.std(self.flux_working[idx])
            self.continuum = cont
            return (self.wave_working, self.flux_working)
        else:
            return None

    def estimateSNR(self):
        assert self.normalized, "Please normalize the data first!"

        std_res = np.std(np.abs(np.ones_like(self.flux_working) - self.flux_working))
        mask = np.zeros_like(self.flux_working)
        mask_idx = [i for i in range(len(self.flux_working)) if np.abs(self.flux_working[i] - 1) >= 0.5 * std_res]
        mask[mask_idx] = 1
        data_tuple = (self.wave_working, self.flux_working)
        (cont, idx) = eF.iterate_continuum(data_tuple, mode='polynomial', n=0, min_sigma = 0.1, mask=mask)
        #plt.plot(self.wave_working, self.flux_working)
        #plt.plot(self.wave_working[idx], self.flux_working[idx], marker='o', markersize=2, color='blue')
        #plt.plot(self.wave_working[mask_idx], self.flux_working[mask_idx], marker='x', markersize=2, color='red')
        #plt.show()
        SNR = (1/np.std(self.flux_working[idx]))
        self.SNR = SNR
        return SNR

    def converToVelocity(self,center=None):
        assert self.velocity_center is None, 'Already in velocity frame'
        if center is None:
            center = input('Type center wavelength in AA:')
            center = float(center)

        go_flag = True
        if not np.min(self.wave_working) <= center <= np.max(self.wave_working):
            message = 'Warning! Center is outside working region, continue?[Y/N]'
            go_flag = eF.go_message(message)

        if go_flag:
            self.wave_working = [(self.wave_working[i] - center) / center * cst.c.to('km/s').value
                                 for i in range(len(self.wave_working))]
            self.velocity_center = center
            self.wave_units = 'km/s'
            return (self.wave_working, self.flux_working)
        else:
            return (None, None)

    def converToWavelength(self):
        assert self.velocity_center is not None, 'Already in wavelength frame'
        self.wave_working = [(1 + self.wave_working[i] / cst.c.to('km/s').value) * self.velocity_center
                             for i in range(len(self.wave_working))]
        self.velocity_center = None
        self.wave_units = 'AA'

        return (self.wave_working, self.flux_working)

    def shiftSpectrum(self,v_offset = None):
        if v_offset is None:
            v_offset = input('Type velocity offset in km/s:')
            v_offset = float(v_offset)

        if self.velocity_center is None:
            self.wave_working = self.wave_working * (1 + v_offset / cst.c.to('km/s').value)
        else:
            self.wave_working = self.wave_working + v_offset

        self.total_offset = self.total_offset + v_offset
        return (self.wave_working, self.flux_working)

    def convertToOpticalDepth(self):
        assert self.normalized == True, 'Normalize the spectrum first!'

        if np.min(self.flux_working) <= 0:
            print('Negative points detected and removed!')
            idx = [i for i in range(len(self.flux_working)) if self.flux_working[i] > 0]
            self.wave_working = self.wave_working[idx]
            self.flux_working = self.flux_working[idx]

        self.flux_working = -1 * np.log(self.flux_working)
        self.optical_depth = True
        return (self.wave_working, self.flux_working)

    def convertToFlux(self):
        assert self.optical_depth == True, 'Flux not in optical depth'
        self.flux_working = np.exp(-1 * self.flux_working)
        self.optical_depth = False
        return (self.wave_working, self.flux_working)

    def searchPeaks(self, n=None, prominence=None):
        from scipy.signal import find_peaks, peak_prominences

        x = - self.flux_working

        if prominence is None:
            if self.SNR is None:
                SNR = self.estimateSNR()
            prominence = 3 / self.SNR

        peaks, _ = find_peaks(x, prominence=prominence)
        prominences = peak_prominences(x, peaks)[0]

        if n is not None:
            threshold = np.sort(prominences)[-n]
            idx = [i for i in range(len(peaks)) if prominences[i] >= threshold]
            peaks = peaks[idx]
            prominences = prominences[idx]

        # visualization
        fig, ax = plt.subplots(1,1)
        ax.plot(self.wave_working,self.flux_working, color="blue")
        self.__XYLabel__(ax)
        self.__auxiliaryPlot__(ax)
        ax.plot(self.wave_working[peaks], self.flux_working[peaks], marker="x", markersize=10, color="red", linestyle="")
        #xspan = ax.get_xlim()[1] - ax.get_xlim()[0]
        #yspan = ax.get_ylim()[1] - ax.get_ylim()[0]
        #for i in range(len(peaks)):
        #    idx = peaks[i]
        #   ax.text(self.wave_working[idx] - xspan/20, self.flux_working[idx] - yspan/50
        #          , "{wav:.1f}".format(wav = self.wave_working[idx]))
        plt.show()

        return self.wave_working[peaks]

    def importModel(self, model):
        self.model = model
        print("New model imported!")
        self.showSpectrum()
        return None

    def fitModel(self, stat="LeastSq", opt="NelderMead", apply_mask=False):
        stat_lib = ['LeastSq']
        opt_lib = ['LevMar', 'NelderMead']
        assert stat in stat_lib, "The stat you choose is not available"
        assert opt in opt_lib, "The opt you choose is not available"

        from sherpa.stats import LeastSq
        from sherpa.optmethods import LevMar, NelderMead
        from sherpa.data import Data1D
        from sherpa.fit import Fit

        stat = eval(stat + "()")
        opt = eval(opt + "()")

        x, y = self.wave_working, self.flux_working
        if apply_mask:
            idx = [i for i in range(len(self.masked)) if self.masked[i]==0]
            x=x[idx]
            y=y[idx]

        data2fit = Data1D('data2fit', x, y)
        model2fit = self.model
        fit = Fit(data2fit, model2fit, stat=stat, method=opt)
        result = fit.fit()

        fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True)
        ax1.plot(self.wave_working, self.flux_working, color='k')
        ax1.plot(self.wave_working, model2fit(self.wave_working), color='red')
        self.__XYLabel__(ax1)

        ax2.plot(self.wave_working, model2fit(self.wave_working) - self.flux_working, color='k')
        ax2.plot(self.wave_working, np.zeros_like(self.wave_working), linestyle='--', color='orange')
        if self.SNR is None:
            SNR = self.estimateSNR()
        SNR = self.estimateSNR()

        ax2.plot(self.wave_working, np.ones_like(self.wave_working)/SNR, linestyle='--', color='blue')
        ax2.plot(self.wave_working, np.ones_like(self.wave_working) /(-SNR), linestyle='--', color='blue')
        self.__XYLabel__(ax2)
        ax2.set_ylabel('Residual')
        plt.show()

        message = "Keep this fit?[Y/N]"
        go_flag = eF.go_message(message)

        if go_flag:
            self.model = model2fit
            self.model_fit = True
        else:
            return None

    def raiseParameter(self, module="any", par='any'):
        if not self.model_fit:
            print("="*16)
            print("Warning! These values are not final!")
            print("=" * 16)

        module = module.lower()
        par = par.lower()

        names = []
        values = []

        for i in range(len(self.model.pars)):
            if not self.model.pars[i].hidden:
                tmp_fullname = self.model.pars[i].fullname.lower() + 'any'
                if module in tmp_fullname and par in tmp_fullname:
                    names.append(self.model.pars[i].fullname)
                    values.append(self.model.pars[i].val)
                    print((names[-1], values[-1]))

        return names, values

if __name__ == '__main__':
    filename = '/HD170740/RED_860/HD170740_w860_redu_20140916_O6.fits'
    filename = '/HD169454/RED_860/HD169454_w860_redl_20160808_O12.fits'
    sp = EdiblesSpectrum(filename)
    print("Barycentric Velocity is", sp.v_bary)
    print(sp.target)
    plt.plot(sp.wave, sp.flux, label='Geocentric')

    bary_data = sp.getSpectrum(xmin=7660, xmax=7705, bary=True)

    plt.plot(bary_data[0], bary_data[1], label='Barycentric')
    axes = plt.gca()
    axes.set_xlim([7660, 7705])
    axes.set_ylim([0, 160])
    plt.vlines((7667.021, 7701.093), 0, 160, linestyles='dashed', colors='r')
    plt.legend()
    plt.show()
