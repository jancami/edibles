import numpy as np
from astropy.io import fits
import astropy.constants as cst
import matplotlib.pyplot as plt
import os

import edibles.edibles_settings as Setting
import edibles.src.datahandling as DataHandling
datadir = '/media/EDIBLES/EDIBLES_DATADIR'
#filename = '/HD169454/RED_860/HD169454_w860_redl_20160808_O12.fits'

class EdiblesSpectrum:
    # This object will contain a spectrum from EDIBLES,
    # and a set of methods to operate on the data.

    def __init__(self, filename, panel_name="", isresult=False):
        """
        Filename is relative to the DR4 directory
        """
        self.__clearAll__()
        self.color_map = ["m", "g", "c", "y", "tan", "teal"]

        self.__loadLineList__()

        if not isresult:
            self.loadSpectrum(filename, panel_name=panel_name)

        if isresult:
            resultdir = Setting.resultdir
            filename = resultdir + filename
            self.loadResult(filename)

    def loadSpectrum(self, filename, clear_all=False, panel_name=""):
        if clear_all: self.__clearLabels__()

        filenames = DataHandling.parseInput(1, filename, checklen=False)
        panel_names = DataHandling.parseInput(len(filenames), panel_name)

        # basic info
        for i, filename in enumerate(filenames):
            hdu = fits.open(Setting.datadir + filename)
            header = hdu[0].header
            self.header.append(header)
            self.target.append(header["OBJECT"])
            self.date.append(header["DATE-OBS"])
            self.v_bary.append(header["HIERARCH ESO QC VRAD BARYCOR"])
            self.filename.append(filename)
            spec_name = filename.split('/')[-1]
            self.spec_name.append(spec_name.split('.')[0])
            self.star_name.append(spec_name.split("_")[0])
            self.panel_names.append(panel_names[i])
            # self.bary_wave = self.wave + (self.v_bary/cst.c.to('km/s').value)*self.wave

            # wave and flux
            flux = hdu[0].data
            crval1 = header["CRVAL1"]
            cdelt1 = header["CDELT1"]
            nwave = len(flux)
            grid = np.arange(0, nwave, 1)
            wave = (grid) * cdelt1 + crval1

            self.flux.append(flux)
            self.flux_working.append(flux)
            self.flux_units.append("Arbitrary")
            self.wave.append(wave)
            self.wave_working.append(wave)
            self.wave_working_velocity.append(np.zeros_like(wave))
            self.wave_units.append("AA")
            self.reference_frame.append("geocentric")
            self.masked.append(np.zeros_like(wave))

            # processing marker
            self.bary_corr.append(False)
            self.normalized.append(False)
            self.velocity_center.append(None)
            self.optical_depth.append(False)
            self.SNR.append(None)

            self.spec_count = self.spec_count + 1

        # finishing up
        log_str = ">>loadSpectrum: filename={filename}, clear_all={clear_all}, panel_name={panel_name}"\
                .format(filename=filenames, clear_all=clear_all, panel_name=panel_names)
        self.__addLog__(log_str)

        return None

    def __clearAll__(self):
        # basic info
        self.filename = []
        self.spec_name = []
        self.star_name = []
        self.header = []
        self.target = []
        self.date = []
        self.v_bary = []
        self.panel_names = []

        # key data
        self.flux = []
        self.flux_working = []
        self.flux_units = []
        self.wave = []
        self.wave_working = []
        self.wave_working_velocity=[]
        self.wave_units = []
        self.reference_frame = []
        self.masked = []

        # process marker
        self.bary_corr = []
        self.normalized = []
        self.velocity_center = []
        self.optical_depth = []
        self.SNR = []

        # other attributes that are not list (only one for the whole class)
        self.model = None
        self.model_fit = False
        self.line_marker = np.array([])
        self.log = ""
        self.spec_count = 0
        self.linelist = {}
        self.linelistspecies=[]
        return None

    def popPanels(self,panels=None):
        panels = self.__parsePanelsInput__(panels)
        panels[::-1].sort()
        go_flag = True
        if len(panels) == self.spec_count:
            go_flag = DataHandling.go_message("Pop ALL panels?[Y/N]")

        if go_flag:
            for panel in panels:
                self.target.pop(panel)
                self.date.pop(panel)
                self.v_bary.pop(panel)
                self.filename.pop(panel)
                self.spec_name.pop(panel)
                self.star_name.pop(panel)
                self.panel_names.pop(panel)

                # wave and flux
                self.flux.pop(panel)
                self.flux_working.pop(panel)
                self.flux_units.pop(panel)
                self.wave.pop(panel)
                self.wave_working.pop(panel)
                self.wave_working_velocity.pop(panel)
                self.wave_units.pop(panel)
                self.reference_frame.pop(panel)
                self.masked.pop(panel)

                # processing marker
                self.bary_corr.pop(panel)
                self.normalized.pop(panel)
                self.velocity_center.pop(panel)
                self.optical_depth.pop(panel)
                self.SNR.pop(panel)
                self.spec_count = self.spec_count - 1
        panels_str = np.array2string(panels, separator=", ")
        log_str = ">>popPanels: panels="
        self.__addLog__(log_str + panels_str)
        return None

    def duplicatePanels(self, ntimes=1, panels=None):
        panels = self.__parsePanelsInput__(panels)
        panels[::-1].sort()
        ntimes = DataHandling.parseInput(len(panels),ntimes)

        for i,panel in enumerate(panels):
            counter = 0
            while counter < ntimes[i]:
                self.target.insert(panel, self.target[panel])
                self.date.insert(panel, self.date[panel])
                self.v_bary.insert(panel, self.v_bary[panel])
                self.filename.insert(panel, self.filename[panel])
                self.spec_name.insert(panel, self.spec_name[panel])
                self.star_name.insert(panel, self.star_name[panel])
                self.panel_names.insert(panel, self.panel_names[panel])

                # wave and flux
                self.flux.insert(panel, self.flux[panel])
                self.flux_working.insert(panel, self.flux_working[panel])
                self.flux_units.insert(panel, self.flux_units[panel])
                self.wave.insert(panel, self.wave[panel])
                self.wave_working.insert(panel, self.wave_working[panel])
                self.wave_working_velocity.insert(panel, self.wave_working_velocity[panel])
                self.wave_units.insert(panel, self.wave_units[panel])
                self.reference_frame.insert(panel, self.reference_frame[panel])
                self.masked.insert(panel, self.masked[panel])

                # processing marker
                self.bary_corr.insert(panel, self.bary_corr[panel])
                self.normalized.insert(panel, self.normalized[panel])
                self.velocity_center.insert(panel, self.velocity_center[panel])
                self.optical_depth.insert(panel, self.optical_depth[panel])
                self.SNR.insert(panel, self.SNR[panel])
                self.spec_count = self.spec_count + 1
                counter+= 1

        panels_str = np.array2string(panels, separator=", ")
        log_str = ">>duplicatePanels: ntimes={ntimes}, panels=".format(ntimes=ntimes)
        self.__addLog__(log_str + panels_str)
        return None

    def combineSpectrum(self,new_spectrum):
        self.target+= new_spectrum.target
        self.date+= new_spectrum.date
        self.v_bary+= new_spectrum.v_bary
        self.filename+= new_spectrum.filename
        self.spec_name+= new_spectrum.spec_name
        self.star_name+= new_spectrum.star_name
        self.panel_names+= new_spectrum.panel_names

        # wave and flux
        self.flux+= new_spectrum.flux
        self.flux_working+= new_spectrum.flux_working
        self.flux_units+= new_spectrum.flux_units
        self.wave+= new_spectrum.wave
        self.wave_working+= new_spectrum.wave_working
        self.wave_working_velocity+= new_spectrum.wave_working_velocity
        self.wave_units+= new_spectrum.wave_units
        self.reference_frame+= new_spectrum.reference_frame
        self.masked+= new_spectrum.masked

        # processing marker
        self.bary_corr+= new_spectrum.bary_corr
        self.normalized+= new_spectrum.normalized
        self.velocity_center+= new_spectrum.velocity_cente
        self.optical_depth+= new_spectrum.optical_depth
        self.SNR+= new_spectrum.SNR

        # non-list attributes
        self.spec_count+= new_spectrum.spec_count

        self.__addLog__("*"*5 + "Start of Spectra Combining" + "*"*5)
        self.log+= new_spectrum.log
        self.__addLog__("*" * 5 + "End of Spectra Combining" + "*" * 5)

        if self.model is None and new_spectrum.model is not None:
            self.model = new_spectrum.model
        self.model_fit = False
        for item in new_spectrum.linelist.keys():
            if item not in self.linelistspecies:
                self.linelist[item] = new_spectrum.linelist[item]
                self.linelistspecies.append(item)
        return None

    def renamePanels(self, names="", panels=None):
        panels = self.__parsePanelsInput__(panels)
        names = DataHandling.parseInput(len(panels), names)
        for i, panel in enumerate(panels):
            self.panel_names[panel] = names[i]
        panels_str = np.array2string(panels, separator=", ")
        log_str = ">>renamePanels: panels="
        self.__addLog__(log_str + panels_str)
        return None

    def __loadLineList__(self, list=None, reset=False):
        if reset:
            self.linelist = {}

        if list is None:
            list = Setting.edibles_linelist

        with open(list) as f:

            while True:
                line = f.readline().replace("\n", "")
                if not line: break
                (specie, wavelength) = line.split(": ")
                self.linelist[specie] = np.fromstring(wavelength, sep=",")

        for item in self.linelist.keys():
            self.linelistspecies.append(item)

        self.__addLog__(">>addLineList(list='{list}', reset=reset".format(list=list, reset=reset))
        return None

    def addLinelist(self, species=None, wavelengths=None):
        if type(species) is str and wavelengths is not None:
            wavelengths = DataHandling.parseInput(1, wavelengths, checklen=False)
            self.linelist[species] = wavelengths
            self.linelistspecies.append(species)
            wavelengths_str = np.array2string(np.array(wavelengths), separator=",")
        self.__addLog__(">>addLinelist: species='{species}', wavelengths=".format(species=species) + wavelengths_str)

    def addComment(self):
        comment=input("please type your comment here")
        self.__addLog__(comment)

    def __addLog__(self, log_str):
        self.log+= log_str + "\n"

    def __parsePanelsInput__(self, panels):
        if panels is None:
            panels = np.arange(len(self.wave_working))
        else:
            panels = DataHandling.parseInput(1, panels, checklen=False)
            #if panels.__class__ is int:
                #panels = np.array([panels])
            #else:
                #panels = np.array(panels)
        panels = np.array(panels).astype("int32")
        panels = np.unique(panels)
        return panels

    def __parseX__(self, panel=0):
        if self.velocity_center[panel] is None:
            return self.wave_working[panel]
        else:
            return self.wave_working_velocity[panel]

    def addMask(self, n=None, panels=None):
        panels = self.__parsePanelsInput__(panels)
        n = DataHandling.parseInput(len(panels), n)

        for i in range(len(panels)):
            panel = int(panels[i])
            x = self.__parseX__(panel=panel)
            y = np.copy(self.flux_working[panel])
            masked = np.zeros_like(x)

            #######
            ##self.showSpectrum(panels = panel)
            #######

            while True:
                n_now = n[i]
                if n_now is None:
                    n_now = input('How many regions should be masked for Panel {panel}?'.format(panel=panel))
                    n_now = np.int(n_now)

                if n_now == 0:
                    break

                else:
                    # an interactive backends will be needed
                    import matplotlib
                    matplotlib.use('Qt5Agg', warn=False, force=True)
                    import matplotlib.pyplot as tmp_plt

                    fig1, ax = tmp_plt.subplots(1, 1)
                    ax.plot(x, y)
                    ax.set_xlim(1.1 * np.min(x) - 0.1 * np.max(x), 1.1 * np.max(x) - 0.1 * np.min(x))
                    self.__XYLabel__(ax=ax, panel=panel)

                    print("Select the boundaries of regions to mask\nLeft = add; Right = pop; Middle = stop")
                    points = tmp_plt.ginput(2*n_now, mouse_add=1, mouse_pop=3, mouse_stop=2)

                    # find the points on spectrum that are most close to the mouse-clicks
                    points_x, points_y = [],[]
                    for j in range(len(points)):
                        points_x.append(points[j][0])
                        points_y.append(points[j][1])

                    points_idx = DataHandling.nearest_point((points_x, points_y), (x, y))
                    points = x[points_idx]

                    # identify the region to be masked
                    n_pairs = int(np.floor(len(points)/2))
                    boundary_left = np.array([])
                    boundary_right = np.array([])
                    for i in range(n_pairs):
                        boundary_left = np.append(boundary_left, points[i*2])
                        boundary_right = np.append(boundary_right, points[i*2+1])

                    (masked, idx_masked) = DataHandling.within_boundaries(x, boundary_left, boundary_right)

                    # switch back to the normal backend and user interaction
                    tmp_plt.close(fig1)
                    matplotlib.use("module://backend_interagg")
                    import matplotlib.pyplot as plt

                    fig2, ax = plt.subplots(1,1)
                    ax.plot(x, y)
                    idx_masked_split = DataHandling.continous_idx(idx_masked)
                    for i in range(len(idx_masked_split)):
                        ax.plot(x[idx_masked_split[i]], y[idx_masked_split[i]],marker='x', markersize=3,color='red')
                        self.__XYLabel__(ax, panel=panel)
                    plt.show()
                    message = 'Remove these points from fitting?[Y/N]'
                    if DataHandling.go_message(message):
                        self.__addMaskValue__(boundary_left, boundary_right, panels=panel)
                        break
                    else:
                        n_now = None
                # if n_now > 0
            # while True
        # for panel in panels
        return None

    def __addMaskValue__(self,boundary_left, boundary_right, panels=0):
        self.resetMask(panels=panels)

        x = self.__parseX__(panel=panels)
        boundary_left, boundary_right = np.array(boundary_left), np.array(boundary_right)
        (masked, idx_masked) = DataHandling.within_boundaries(x, boundary_left, boundary_right)
        self.masked[panels] = masked

        boundary_left_str = np.array2string(boundary_left,precision=1, separator=", ")
        boundary_right_str = np.array2string(boundary_right,precision=2, separator=", ")
        log_str = ">>addMask:" \
                  " boundary_left=" + boundary_left_str+ \
                  ", boundary_right=" + boundary_right_str+ \
                  ", panels={panel}".format(panel=panels)
        self.__addLog__(log_str)
        return None

    def baryCorrection(self, panels=None):
        panels = self.__parsePanelsInput__(panels)
        for panel in panels:
            if self.bary_corr[panel] is True:
                print("Panel {panel} already in barycentric reference frame. Skipped.".format(panel=panel))
                panels = panels[(panels != panel)]

        for panel in panels:
            self.wave_working[panel] = DataHandling.velocityShift(self.wave_working[panel], self.v_bary[panel])
            self.wave_working_velocity[panel] = self.wave_working_velocity[panel] + self.v_bary[panel]
            self.bary_corr[panel] = True
            self.reference_frame[panel] = "Barycentric"

        panels_str = np.array2string(panels, separator=", ")
        self.__addLog__(">>baryCorrection: panels=" + panels_str)
        return None

    def converToVelocity(self, center=None, panels=None):
        panels = self.__parsePanelsInput__(panels)
        center = DataHandling.parseInput(len(panels), center)

        for panel in panels:
            if self.velocity_center[panel] is not None:
                print("Panel {panel} already in velocity frame. Skipped.".format(panel=panel))
                panels = panels[(panels != panel)]

        for i, panel in enumerate(panels):
            center_now = center[i]
            if center_now is None:
                center_now = float(input('Type center wavelength in AA:'))

            if not np.min(self.wave_working[panel]) <= center_now <= np.max(self.wave_working[panel]):
                message = "Warning! Center for Panel {panel} is outside working region, continue?[Y/N]".format(panel=panel)
                go_flag = DataHandling.go_message(message)
            else:
                go_flag = True

            if go_flag:
                self.wave_working_velocity[panel] = (self.wave_working[panel] - center_now) / center_now * cst.c.to('km/s').value
                self.velocity_center[panel] = center_now
                self.wave_units[panel] = 'km/s'

        panels_str = np.array2string(panels, separator=", ")
        log_str = ">>converToVelocity: " \
                  "center={center}, " \
                  "panels=".format(center=center)
        self.__addLog__(log_str+panels_str)
        return panels

    def converToWavelength(self, panels=None):
        panels = self.__parsePanelsInput__(panels)
        for panel in panels:
            if self.velocity_center[panel] is None:
                print("Panel {panel} already in wavelength frame. Skipped.".format(panel=panel))
                panels = panels[(panels != panel)]

        for panel in panels:
            self.wave_working[panel] = self.velocity_center[panel] * \
                                       (1 + self.wave_working_velocity[panel] / cst.c.to('km/s').value)
            self.velocity_center[panel] = None
            self.wave_units[panel] = 'AA'

        panels_str = np.array2string(panels, separator=", ")
        self.__addLog__(">>converToWavelength: panels="+panels_str)
        return None

    def convertToOpticalDepth(self, panels=None):
        panels = self.__parsePanelsInput__(panels)
        for panel in panels:
            if self.normalized[panel] == False:
                print("Panel {panel} is not normalized. Skipped.".format(panel=panel))
                panels = panels[panels != panel]
            if self.optical_depth[panel] == True:
                print("Panel {panel} already in optical depth. Skipped.".format(panel=panel))
                panels = panels[panels != panel]

        for panel in panels:
            if np.min(self.flux_working[panel]) <= 0:
                print('Negative points detected and removed in Panel {panel}!'.format(panel=panel))
                idx = [i for i in range(len(self.flux_working[panel])) if self.flux_working[panel][i] > 0]
                self.wave_working[panel] = self.wave_working[panel][idx]
                self.wave_working_velocity[panel] = self.wave_working_velocity[panel][idx]
                self.flux_working[panel] = self.flux_working[panel][idx]
                self.masked[panel] = self.maksed[panel][idx]

            self.flux_working[panel] = -1 * np.log(self.flux_working[panel])
            self.optical_depth[panel] = True

        panels_str = np.array2string(panels, separator=", ")
        self.__addLog__(">>convertToOpticalDepth: panels="+panels_str)

        return panels

    def convertToFlux(self, panels=None):
        panels = self.__parsePanelsInput__(panels)
        for panel in panels:
            if self.optical_depth[panel] == False:
                print("Panel {panel} is not in optical depth. Skipped.".format(panel=panel))
                panels = panels[panels != panel]

        for panel in panels:
            self.flux_working[panel] = np.exp(-1 * self.flux_working[panel])
            self.optical_depth[panel] = False

        panels_str = np.array2string(panels, separator=", ")
        self.__addLog__(">>convertToFlux: panels="+panels_str)
        return None

    def cutSpectrum(self, xmin=None, xmax=None, center=None, span=None, panels=None):
        panels = self.__parsePanelsInput__(panels)
        (xmin, xmax, center, span) = DataHandling.parseInput(len(panels), xmin, xmax, center, span)

        for i in range(len(panels)):
            panel = int(panels[i])
            xmin_now, xmax_now, center_now, span_now = xmin[i], xmax[i], center[i], span[i]
            x = self.__parseX__(panel=panel)

            if center_now is not None and span_now is not None:
                xmin_now = center_now - 0.5 * span_now
                xmax_now = center_now + 0.5 * span_now

            (xmin_now, xmax_now) = self.__getMinMax__(xmin_now, xmax_now, panel=panel)

            idx = (x > xmin_now) * (x < xmax_now)
            self.wave_working[panel] = self.wave_working[panel][idx]
            self.wave_working_velocity[panel] = self.wave_working_velocity[panel][idx]
            self.flux_working[panel] = self.flux_working[panel][idx]
            self.masked[panel] = self.masked[panel][idx]

            log_str = ">>cutSpectrum: " \
                      "xmin={xmin:.2f}, xmax={xmax:.2f}, panels={panel}".format(xmin=xmin_now, xmax=xmax_now,
                                                                                panel=panel)
            self.__addLog__(log_str)
        return None

    def cutSpectrum_visual(self, panels=None):
        panels = self.__parsePanelsInput__(panels)

        # an interactive backends will be needed
        import matplotlib
        matplotlib.use('Qt5Agg', warn=False, force=True)
        import matplotlib.pyplot as tmp_plt

        for panel in panels:
            x = self.__parseX__(panel=panel)
            fig1, ax = tmp_plt.subplots(1, 1)
            ax.plot(x, self.flux_working[panel])
            ax.set_xlim(1.1 * np.min(x) - 0.1 * np.max(x), 1.1 * np.max(x) - 0.1 * np.min(x))
            self.__XYLabel__(ax, panel=panel)

            print("Select the boundaries of working region")
            points = tmp_plt.ginput(2, mouse_add=1, mouse_pop=3, mouse_stop=2)

            tmp_plt.close(fig1)
            points_tuple = ([points[0][0], points[1][0]], [points[0][1], points[1][1]])
            cut_idx = DataHandling.nearest_point(points_tuple, (x, self.flux_working[panel]))
            (xmin, xmax) = (x[cut_idx[0]], x[cut_idx[1]])
            self.cutSpectrum(xmin=xmin, xmax=xmax, panels=panel)

        # now let's switch back to the normal backend
        matplotlib.use("module://backend_interagg")
        import matplotlib.pyplot as plt
        return None

    def cutSpectrumTo(self, cut=None, to=0, mode="wavelength"):
        cuts = self.__parsePanelsInput__(cut)
        tos, modes = DataHandling.parseInput(len(cuts), to, mode)

        for i, to in enumerate(tos):
            assert modes[i].lower() in ["wavelength", "velocity"], "mode should be 'wavelength' or 'velocity'"
            assert self.wave_working_velocity[to][0] < self.wave_working_velocity[to][-1], \
                "Panel {panel} does not have velocity grid".format(panel=to)

        for i, cut in enumerate(cuts):
            if modes[i] == "wavelength":
                xmin, xmax = np.min(self.wave_working[to[i]]), np.max(self.wave_working[to[i]])
                convert_trigger = False
                if self.velocity_center[cut] is not None:
                    v_center = np.copy(self.velocity_center[cut])
                    self.converToWavelength(panels=cut)
                    convert_trigger = True
                self.cutSpectrum(xmin=xmin, xmax=xmax, panels=cut)
                if convert_trigger:
                    self.converToVelocity(center = v_center)
            else:
                xmin, xmax = np.min(self.wave_working_velocity[to[i]]), np.max(self.wave_working_velocity[to[i]])
                convert_trigger = False
                if self.velocity_center[cut] is None:
                    self.converToVelocity(panels=cut)
                    convert_trigger = True
                self.cutSpectrum(xmin=xmin, xmax=xmax, panel=cut)
                if convert_trigger:
                    self.converToWavelength(panels=cut)
        return None

    def __getMinMax__(self, xmin, xmax, panel=0):
        assert xmin < xmax, "xmin must be smaller than xmax!"

        x = self.__parseX__(panel=panel)
        if not DataHandling.within_boundaries(xmin, np.min(x), np.max(x))[0]:
            xmin = np.min(x)
            print("Warning, x_min outside working region, reset to current min")

        if not DataHandling.within_boundaries(xmax, np.min(x), np.max(x))[0]:
            xmax = np.max(x)
            print("Warning, x_max outside working region, reset to current max")

        return (xmin, xmax)

    def estimateSNR(self, panels=None):
        panels = self.__parsePanelsInput__(panels)
        SNR_all = np.array([])
        for panel in panels:
            if self.normalized[panel] is False:
                print("Panel {panel} not normalized. Skipped.".format(panel=panel))
                panels = panels[panels != panel]

        for panel in panels:
            x = self.__parseX__(panel=panel)
            (SNR, idx_SNR) = DataHandling.estimateSNR(x, self.flux_working[panel], silenced=True)
            self.SNR[panel] = SNR
            SNR_all = np.append(SNR_all, SNR)
            self.__addLog__("SNR for panel {panel}: {SNR:.2f}".format(panel=panel, SNR=SNR))

        panels_str = np.array2string(panels, separator=", ")
        self.__addLog__(">>estimateSNR: panels="+panels_str)

        return SNR_all

    def fitContinuum(self, mode='spline', n=3,
                     lsigma=1, usigma=2, iterates=30, min_sigma = 0.1,
                     silence=False, apply_mask=False, panels=None):

        panels = self.__parsePanelsInput__(panels)
        (mode, n, lsigma, usigma, iterates, min_sigma, silence, apply_mask) = \
            DataHandling.parseInput(len(panels), mode, n, lsigma, usigma, iterates, min_sigma, silence, apply_mask)

        for i in range(len(panels)):
            panel = panels[i]
            x = self.__parseX__(panel=panel)
            data_tuple = (x, self.flux_working[panel])
            if apply_mask[i]: mask = np.copy(self.masked[panel])
            else: mask = np.zeros_like(x)
            (cont, idx_cont) = DataHandling.iterate_continuum(data_tuple,mode=mode[i], n=n[i],
                                           lsigma=lsigma[i], usigma=usigma[i],
                                           iterates=iterates[i], min_sigma=min_sigma[i],
                                           mask = mask)

            # making plot
            go_flag = True
            if not silence[i]:
                #fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, sharex=True)
                fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True)
                # origional
                ax1.plot(x, self.flux_working[panel], color='k')
                ax1.plot(x, cont(x), color='blue')
                ax1.set_ylabel("Original")

                # normalized
                ax2.plot(x, self.flux_working[panel] / cont(x), color='k')
                ax2.plot(x, np.ones_like(x), linestyle='--', color='blue')
                ax2.set_ylabel("Normalized")
                self.__XYLabel__(ax=ax2, x_label=False, panel=panel)

                # highlights
                idx_cont_split = DataHandling.continous_idx(idx_cont)
                for j in range(len(idx_cont_split)):
                    idx_plot = idx_cont_split[j]
                    ax1.plot(x[idx_plot], self.flux_working[panel][idx_plot],
                             linestyle='', marker='o', markersize=2, color='red')
                    ax2.plot(x[idx_plot], self.flux_working[panel][idx_plot] / cont(x[idx_plot]), color='red')

                # residual
                #std_res = self.flux_working - cont(self.wave_working)
                #ax3.plot(self.wave_working, std_res, color='k')
                #std_res = np.std(std_res[idx_cont])
                #ax3.plot(self.wave_working, np.ones_like(self.wave_working) * std_res, linestyle='--', color='blue')
                #ax3.plot(self.wave_working, np.ones_like(self.wave_working) * (-std_res), linestyle='--', color='blue')
                #ax3.set_ylim([-std_res * 2, std_res * 2])
                #ax3.set_ylabel("Residual")

                ax1.grid()
                ax2.grid()
                #ax3.grid()
                plt.show()

                go_flag = DataHandling.go_message('Keep this continuum?[Y/N]')

            if go_flag:
                self.flux_working[panel] = self.flux_working[panel] / cont(x)
                self.normalized[panel] = True

                log_str = ">>fitContinuum: " \
                          "mode='{mode}', n={n}, lsigma={lsigma}, usigma={usigma}, " \
                          "iterates={iterates}, min_sigma={min_sigma}, " \
                          "silence={silence}, apply_mask={apply_mask}, panels={panel}" \
                          .format(mode=mode[i], n=n[i], lsigma=lsigma[i], usigma=usigma[i],
                                  iterates=iterates[i], min_sigma=min_sigma[i],
                                  silence = silence[i], apply_mask=apply_mask[i], panel=panel)
                self.__addLog__(log_str)

                # add SNR only if fit with polynomial
                if mode[i].lower() == "polynomial" or mode[i].lower() == "p":
                    res_norm = np.ones_like(self.flux_working[panel]) - self.flux_working[panel]
                    std_res_norm = np.std(res_norm[idx_cont])
                    self.SNR[panel] = 1. / std_res_norm

        return None

    def fitModel(self, stat="LeastSq", opt="NelderMead", apply_mask=False, panels=None):
        stat_lib = ['LeastSq']
        opt_lib = ['LevMar', 'NelderMead']
        assert stat in stat_lib, "The stat you choose is not available"
        assert opt in opt_lib, "The opt you choose is not available"

        from sherpa.stats import LeastSq
        from sherpa.optmethods import LevMar, NelderMead
        from sherpa.data import Data1D
        from sherpa.fit import Fit

        stat_apply = eval(stat + "()")
        opt_apply = eval(opt + "()")

        panels = self.__parsePanelsInput__(panels)
        x,y,masked = np.array([]),np.array([]),np.array([])
        for panel in panels:
            x = np.append(x, self.wave_working[panel])
            y = np.append(y, self.flux_working[panel])
            masked = np.append(masked, self.masked[panel])

        x_uniq = np.unique(x)
        if len(x_uniq) < len(x):
            idx_uniq = np.array([]).astype("int64")
            for i in range(len(x_uniq)):
                idx_uniq = np.append(idx_uniq, (x == x_uniq[i]).nonzero()[0][0])
            x, y, masked = x[idx_uniq], y[idx_uniq], masked[idx_uniq]

        if apply_mask:
            idx = (masked == 0).nonzero()
            x=x[idx]
            y=y[idx]

        data2fit = Data1D('data2fit', x, y)
        model2fit = self.model
        fit = Fit(data2fit, model2fit, stat=stat_apply, method=opt_apply)
        result = fit.fit()

        # make result plots
        n_rows = len(panels)
        v_centers = []
        for panel in panels:
            v_centers.append(self.velocity_center[panel])

        if None in v_centers:
            fig, axs = plt.subplots(nrows=n_rows, ncols=2)
        else:
            fig, axs = plt.subplots(nrows=n_rows, ncols=2, sharex="row")

        for i in range(n_rows):
            panel = panels[i]
            x_plot = self.__parseX__(panel=panel)
            # result
            axs = fig.axes[i*2]
            axs.plot(x_plot, self.flux_working[panel], color="k")
            axs.plot(x_plot, np.ones_like(x_plot), linestyle='--', color="orange")
            axs.plot(x_plot, model2fit(self.wave_working[panel]), color="red")
            axs.grid()
            self.__XYLabel__(ax=axs, x_label=False, panel=panel)
            if apply_mask and 1 in self.masked[panel]:
                axs.scatter(x_plot[self.masked[panel] == 1], self.flux_working[panel][(self.masked[panel] == 1)],
                            marker="x", color="red")
            if i == n_rows - 1:
                axs.set_xlabel("Fitted Model")

            # residual
            axs = fig.axes[i*2 + 1]
            axs.plot(x_plot, model2fit(self.wave_working[panel]) - self.flux_working[panel], color='k')
            axs.plot(x_plot, np.zeros_like(x_plot), color="orange")
            if self.SNR[panel] is None:
                SNR = self.estimateSNR(panels = panel)
            axs.plot(x_plot, np.ones_like(x_plot)/self.SNR[panel], linestyle='--', color="orange")
            axs.plot(x_plot, -1 * np.ones_like(x_plot) / self.SNR[panel], linestyle='--', color="orange")
            axs.grid()
            if apply_mask and 1 in self.masked[panel]:
                axs.scatter(x_plot[self.masked[panel] == 1], model2fit( self.wave_working[panel][self.masked[panel] == 1] ),
                            marker="x", color="red")
            if i == n_rows - 1:
                axs.set_xlabel("Residuals")
        plt.show()

        if DataHandling.go_message("Keep this fit?[Y/N]"):
            self.model = model2fit
            self.model_fit = True

            panels_str = np.array2string(panels, separator=", ")
            log_str = ">>fitModel: " \
                      "stat='{stat}', opt='{opt}', apply_mask={apply_mask}, panels=" \
                      .format(stat=stat, opt=opt, apply_mask=apply_mask)
            self.__addLog__(log_str + panels_str)
            self.__addLog__(str(model2fit))
        else:
            return result

    def printHeader(self,panels=None):
        panels = sp.__parsePanelsInput__(panels)
        for i, panel in enumerate(panels):
            print("="*10 + " Header {all}/{now} ".format(all=len(panels), now=i) + "="*10)
            print(self.header[panel])

    def showLog(self):
        print(self.log)

    def showSpectrum(self, x_label=True, y_label=True, continuum=True,
                     model=True, masked=True,
                     highlight="All", panels=None,
                     save=False, savepath=None, filename=None):
        panels = self.__parsePanelsInput__(panels)
        self.__makeBasicPlots__(x_label=x_label, y_label=y_label,
                                continuum=continuum, model=model,
                                masked=masked, highlight=highlight,
                                panels=panels)
        if save:
            if savepath is None:savepath = Setting.plotdir
            if not os.path.exists(savepath):os.makedirs(savepath)
            if filename is None:filename = input("Please type file name to be used:")
            plt.savefig(savepath + filename+".png")
        plt.show()
        return None

    def __makeBasicPlots__(self,x_label=True, y_label=True, continuum=True, model=True, masked=True, highlight="All", panels=None):
        panels = self.__parsePanelsInput__(panels)
        n_plots = len(panels)

        x_label_backup = x_label
        v_centers = []
        for panel in panels:
            v_centers.append(self.velocity_center[panel])
        if None in v_centers:
            fig, axs = plt.subplots(nrows=n_plots)
        else:
            fig, axs = plt.subplots(nrows=n_plots, sharex=True)
            x_label = False
        if v_centers == [None]*n_plots: x_label=False

        x_labels, y_labels, continua, models, maskeds, highlights, x_label_backup = DataHandling.parseInput(n_plots, x_label, y_label,
                                                                                            continuum, model, masked,
                                                                                            highlight, x_label_backup)

        for i in range(n_plots):
            panel = panels[i]
            x = self.__parseX__(panel=panel)
            fig = plt.gcf()
            axs = fig.axes[i]
            axs.plot(x, self.flux_working[panel], color="k")
            axs.grid()
            self.__XYLabel__(ax=axs, x_label=x_labels[i], y_label=y_labels[i],panel=panel)
            self.__auxiliaryPlot__(ax=axs, continuum=continua[i], model=models[i], masked=maskeds[i],highlight=highlights[i], panel=panel)
            if i == n_plots-1 and (None not in v_centers or v_centers == [None]*n_plots):
                self.__XYLabel__(ax=axs, x_label=x_label_backup[i], y_label=y_labels[i],panel=panel)

        return None

    def __XYLabel__(self, ax=None, x_label=True, y_label=True, panel=0):
        # This is all about axis and labels
        if ax is None: ax = plt.gca()

        if x_label:
            if type(x_label) in [type("abc"), type(np.array(["abc"])[0])]:
                ax.set_xlabel(x_label)
            else:
                if self.velocity_center[panel] is None:ax.set_xlabel('AA')
                else:ax.set_xlabel('km/s')

        if y_label:
            if type(y_label) in [type("abc"), type(np.array(["abc"])[0])]:
                ax.set_ylabel(y_label)
                if self.optical_depth[panel]: ax.invert_yaxis()
            else:
                y_label_str = "Panel {panel}\n".format(panel=panel)
                if self.panel_names[panel] != "":
                    y_label_str = y_label_str + self.panel_names[panel] + "\n"

                if self.optical_depth[panel]:
                    y_label_str = y_label_str + "Opt. Dep."
                    ax.invert_yaxis()
                else:
                    if self.normalized[panel]:
                        y_label_str = y_label_str + "Norm. Flux"
                    else:
                        y_label_str = y_label_str + "Flux"
                ax.set_ylabel(y_label_str)
        return None

    def __auxiliaryPlot__(self, ax=None, continuum=True, model=True, masked=True, highlight="All", panel=0):
        # trivial matters, continuum at unity, model, masked; marked lines to be added
        if ax is None: ax = plt.gca()
        x = self.__parseX__(panel=panel)

        # continuum
        if self.normalized[panel] and continuum:
            if self.optical_depth[panel]:
                ax.plot(x, np.zeros_like(x), linestyle='--', color='orange')
            else:
                ax.plot(x, np.ones_like(x), linestyle='--', color='orange')
        # model
        if self.model is not None and model:
            ax.plot(x, self.model(self.wave_working[panel]), color='blue')
        # masked points
        if np.sum(self.masked[panel]) > 0 and masked:
            idx_masked = (self.masked[panel] == 1).nonzero()[0]
            ax.scatter(x[idx_masked], self.flux_working[panel][idx_masked], marker="x", color="red")

        # marked lines
        if highlight:
            if highlight.lower() == "all": highlight = self.linelistspecies
            xmin, xmax = np.min(self.wave_working[panel]), np.max(self.wave_working[panel])
            ymin, ymax = np.min(self.flux_working[panel]), np.max(self.flux_working[panel])
            for i, species in enumerate(highlight):
                if species not in self.linelistspecies:
                    print("Line list for " + species +
                          " not found. Available species include: {list}".format(list = self.linelistspecies))
                else:
                    wavelengths = self.linelist[species]
                    within, idx = DataHandling.within_boundaries(wavelengths, xmin, xmax)
                    if np.sum(within) > 0:
                        wavelengths = wavelengths[idx]
                        if self.velocity_center[panel] is not None:
                            wavelengths = (wavelengths - self.velocity_center[panel]) / self.velocity_center[panel]\
                                           * cst.c.to('km/s').value
                        for wavelength in wavelengths:
                            ax.plot([wavelength, wavelength], [ymax, ymin], color=self.color_map[i])
                        ax.text(np.mean(wavelengths),ymin,species, color=self.color_map[i])
        return None

    def searchPeaks(self, n=None, prominence=3, panels=None):
        panels = self.__parsePanelsInput__(panels)
        SNR = []
        for panel in panels:
            assert self.normalized[panel] is True,"Panel {panel} has not been normalized!".format(panel=panel)
            if self.SNR[panel] is None:
                self.estimateSNR(panels=panel)
            SNR.append(self.SNR[panel])
        SNR = np.max(SNR)

        x,y = np.array([]),np.array([])
        for panel in panels:
            x = np.append(x, self.wave_working[panel])
            y = np.append(y, self.flux_working[panel])

        x_uniq = np.unique(x)
        if len(x_uniq) < len(x):
            idx_uniq = np.array([]).astype("int64")
            for i in range(len(x_uniq)):
                idx_uniq = np.append(idx_uniq, (x == x_uniq[i]).nonzero()[0][0])
            x, y = x[idx_uniq], y[idx_uniq]

        peak_idx = DataHandling.searchPeak(y, n_peaks=n, normalized=True, SNR=SNR, prominence=prominence)
        peak_wavelengths = x[peak_idx]
        peak_fluxs = y[peak_idx]

        n_rows = len(panels)
        v_centers = []
        for panel in panels:
            v_centers.append(self.velocity_center[panel])

        if None in v_centers:
            fig, axs = plt.subplots(nrows=n_rows)
        else:
            fig, axs = plt.subplots(nrows=n_rows, sharex=True)

        for i, panel in enumerate(panels):
            x_plot = self.__parseX__(panel=panel)
            # result
            axs = fig.axes[i]
            axs.plot(x_plot, self.flux_working[panel], color="k")
            self.__XYLabel__(ax=axs, panel=panel)
            self.__auxiliaryPlot__(ax=axs, masked=False, model=False, highlight=False, panel=panel)

            xmin, xmax = np.min(self.wave_working[panel]), np.max(self.wave_working[panel])
            within, idx = DataHandling.within_boundaries(peak_wavelengths, xmin, xmax)
            if np.sum(within) > 0:
                peak_wavelengths_plot = peak_wavelengths[idx]
                peak_fluxs_plot = peak_fluxs[idx]
                if self.velocity_center[panel] is not None:
                    peak_wavelengths_plot = (peak_wavelengths_plot - self.velocity_center[panel]) / self.velocity_center[panel] \
                                  * cst.c.to('km/s').value
                axs.plot(peak_wavelengths_plot, peak_fluxs_plot, marker="x", markersize=10, color="red", linestyle="")

        return peak_wavelengths

    def importModel(self, model):
        self.model = model
        self.__addLog__(">>importModel: ")
        self.__addLog__(str(model))
        print("New model imported!")
        print(model)
        self.showSpectrum(masked=False)
        return None

    def resetMask(self, panels = None):
        panels = self.__parsePanelsInput__(panels)
        for panel in panels:
            self.masked[panel] = np.zeros_like(self.wave_working[panel])

        panels_str = np.array2string(panels, separator=", ")
        self.__addLog__(">>resetMask: panels=" + panels_str)

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

    def resetSpectrum(self, panels=None):
        panels = self.__parsePanelsInput__(panels)
        for i in range(len(panels)):
            panel = panels[i]
            message = "All handling of Panel {panel} will be abort, continue?[Y/N]".format(panel=panel)
            if DataHandling.go_message(message):
                self.wave_working[panel] = np.copy(self.wave[panel])
                self.wave_units[panel] = "AA"
                self.wave_working_velocity[panel] = np.zeros_like(self.wave_working[panel])
                self.flux_working[panel] = np.copy(self.flux[panel])
                self.flux_units[panel] = "Arbitrary"
                self.masked[panel] = np.zeros_like(self.wave_working[panel])

                self.reference_frame[panel] = "geocentric"
                self.bary_corr[panel] = False
                self.normalized[panel] = False
                self.velocity_center[panel] = None
                self.optical_depth[panel] = False
                self.SNR[panel] = None
                self.__addLog__(">>resetSpectrum: panels={panel}".format(panel=panel))
        return None

    def shiftSpectrum(self,v_offset = None, panels=None):
        panels = self.__parsePanelsInput__(panels)
        v_offset = DataHandling.parseInput(len(panels), v_offset)

        for i in range(len(panels)):
            panel = panels[i]
            v_offset_now = v_offset[i]
            if v_offset_now is None:
                v_offset_now = float(input('Type velocity offset for Panel {panel} in km/s:'.format(panel=panel)))

            self.wave_working_velocity[panel] = self.wave_working_velocity[panel] + v_offset_now
            self.wave_working[panel] = self.wave_working[panel] * (1 + v_offset_now / cst.c.to('km/s').value)

            self.reference_frame[panel] = "customized"

            log_str = ">>shiftSpectrum: v_offset={v_offset}, panels={panel}".format(v_offset=v_offset_now, panel=panel)
            self.__addLog__(log_str)
        return None

    def statusReport(self, panels=None):
        panels = self.__parsePanelsInput__(panels)
        for panel in panels:
            # header
            print("="*15 + "  Panels {panel}  ".format(panel=panel) + "="*15)

            # filename and panel name
            if self.panel_names[panel] != "":
                print("Panel Name: " + self.panel_names[panel])
            print(self.filename[panel])

            # X
            x_str = "X: "
            x = self.__parseX__(panel=panel)
            xmin, xmax = np.min(x), np.max(x)
            if self.velocity_center[panel] is None:
                x_str+= "Wavelength, between {xmin:.2f} and {xmax:.2f} AA".format(xmin=xmin, xmax=xmax)
            else:
                x_str+= "Velocity Grid, between {xmin:.2f} and {xmax:.2f} km/s, ".format(xmin=xmin, xmax=xmax)
                x_str+= "center at {center} AA".format(center = self.velocity_center[panel])
            print(x_str)

            #Y
            y_str = "Y: "
            if not self.normalized[panel]:
                y_str+= "Flux"
            else:
                if self.optical_depth:
                    y_str+= "Optical Depth"
                else:
                    y_str+= "Normalized Flux"
            print(y_str)

            # frame
            print("In "+ self.reference_frame[panel] + " frame")
            print("")

    def outputResult(self, path=None, filename=None, panels=None):
        panels = self.__parsePanelsInput__(panels)
        if path is None:
            path = Setting.resultdir

        if not os.path.exists(path):
            os.makedirs(path)

        if filename is None:
            filename = self.star_name[0] + "_Results.txt"

        f = open(path + filename, "w")
        f.write(self.log)
        f.write("=" * 20 + "\n")
        f.close()
        self.outputData(path=path, filename=filename, panels=panels)

        return None

    def outputData(self, path=None, filename=None, panels=None):
        panels = self.__parsePanelsInput__(panels)
        if path is None:
            path = Setting.resultdir

        if not os.path.exists(path):
            os.makedirs(path)

        if filename is None:
            filename = self.star_name[0] + "_Data.txt"

        f = open(path + filename, "a")

        for panel in panels:
            f.write("=" * 10 + "  Panel {panel}".format(panel=panel) + "=" * 10 + "\n")

            x = self.__parseX__(panel)
            n_column = 2
            data_array = np.c_[x, self.flux_working[panel]]
            if self.velocity_center[panel] is None:
                header = "Wavelength, Flux"
            else:
                header = "Velocity, Flux"
                f.write("=== X in velocity scale, centered at {center:.2f} AA ===\n".format(
                    center=self.velocity_center[panel]))

            if self.model is not None:
                n_column += 1
                header += ", Model"
                data_array = np.c_[data_array, self.model(self.wave_working[panel])]

            f.write(header + "\n")
            for i in range(data_array.shape[0]):
                data_str = np.array2string(data_array[i, :], separator=", ")[1:-1] + "\n"
                f.write(data_str)

            f.write("\n")
        f.close()
        return None

    def loadResult(self, filename):
        with open(filename) as f:
            leap_all = [0]
            while True:
                line=f.readline().replace("\n","")
                if line == "="*20: break
                if line == "*" * 5 + "Start of Spectra Combining" + "*" * 5:
                    leap_all.append(self.spec_count)
                if line == "*" * 5 + "End of Spectra Combining" + "*" * 5:
                    leap_all.pop(-1)

                if line[0:2] == ">>":
                    line = line.replace(">>", "")
                    if ": " in line:
                        command, parameter = line.split(": ")
                    else:
                        command, parameter = line, ""

                    command = self.__parseCommand__(command, parameter, leap=leap_all[-1])
                    if command is not None: eval("self." + command)

        return None

    def __parseCommand__(self,command, parameter, leap=0):
        sp_methods = ["fitModel", "importModel", "fitContinuum", "addMask"]
        if command in sp_methods:
            if command in ["fitModel", "importModel"]: command = None
            if command == "fitContinuum" and "silence=False" in parameter:
                parameter.replace("silence=False", "silence=True")
            if command == "addMask": command = "__addMaskValue__"

        if command is not None:
            if parameter == "":
                command+= "()"
            else:
                if leap and "panels=" in parameter:
                    parameter_1st, parameter_2nd = parameter.split("panels=")

                    if "[" in parameter_2nd:
                        parameter_2nd = parameter_2nd.replace("[", "").replace("]", "")
                        parameter_2nd = np.fromstring(parameter_2nd, dtype=int, sep=",") + leap
                        panels_str = np.array2string(parameter_2nd, separator=", ")
                    else:
                        parameter_2nd = int(parameter_2nd) + leap
                        panels_str = str(parameter_2nd)

                    parameter = parameter_1st + "panels=" + panels_str

                command = command + "(" + parameter + ")"
        return command

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
