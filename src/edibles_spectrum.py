# +
# NAME:
#     edibles_spectrum
#
# PURPOSE:
#     Decompose observed Spectra into local continuua and absorption
#     features, using Gaussian and Voigt profiles and returning the
#     fit parameters for each.
#
# CALLING SEQUENCE:
#     from edibles_spectrum import dibsFit
#
# INPUTS:
#     wave:         he wavelength in angstrom, In the desired rest-frame
#     flux:         The flux intensity. Spectrum should be normalized to 1,
#                   not nesecerily the local normalization
#     central_wave: The list containing the central wavelengths of the pick
#                   positions. If the list contain only one value, programm
#                   take it as central wavelength and return +-10 A around it
#                   but if it contain two value program take those as min and
#                   max values
#
# OUTPUT:
#     initial_val:  The initial values for fitting
#     params:       The fitted parameters
#     lambda_fit:   The used region for fitting
#
# AUTHOR:
#     Amin Farhang
#     University of Western Ontario
#
# DATE:
#     V1.1
#     March 7 2018
#
# +

import os, glob, sys
import mpfit
import ref_index
from scipy.special import wofz
from scipy.stats import chi2
import pyfits
import numpy as np
import matplotlib.pylab as plt
import warnings
warnings.simplefilter('ignore')




# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%  Spectrum  %%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class getSpectrum:

    EDIBLES_DR3 = '/Users/amin/DR3/'

    # =========================================================================
    #   initialization -- Define the required initial arguments for using
    #                     in the dibsFit class. The vel_wave is the central
    #                     wavelength for converting the angstrom into km/s
    #                     central_wave should be a list of data
    # =========================================================================
    def __init__(self, target_name, central_wave, obs_date=None, origin=None, vel_ref=None, normal=None, vel_wave=0):
        self.target_name = target_name
        self.central_wave = central_wave
        self.obs_date = obs_date
        self.origin = origin
        self.vel_ref = vel_ref
        self.normal = normal
        self.vel_wave = vel_wave

        # Initialize the basic properties use print
        # objectName.wave_unit to print this out
        self.wave_unit = 'unknown'
        self.flux_unit = 'unknonw'
        self.reference_frame = 'geocentric'

        # output data
        self.wave = []
        self.flux = []


    # =========================================================================
    #   findSpectrum -- Find the name and folder of spectrum if origin is not
    #                   defined. This Method is valid only for EDIBLES survey
    #                   for other fits file determine the origin of file
    # =========================================================================
    def findSpectrum(self):

        # check wavelength range is valid
        if self.central_wave == None:
            print 'central wavelength is not defined!'
            sys.exit()

        # central wavelength
        if len(self.central_wave) == 1:
            xmin = self.central_wave[0] - 10
            xmax = self.central_wave[0] + 10

        # wavelength range
        if len(self.central_wave) == 2:
            xmin = self.central_wave[0]
            xmax = self.central_wave[1]

        # select arm
        arms = ['BLUE_346', 'BLUE_437', 'REDL_564', 'REDU_564', 'REDL_860', 'REDU_860']
        if xmin <= 3876: l1 = 0
        if xmin <= 4990 and xmin >= 3754: l1 = 1
        if xmin <= 5667.9 and xmin >= 4616: l1 = 2
        if xmin <= 6693.9 and xmin >= 5668: l1 = 3
        if xmin <= 8649.9 and xmin >= 6694: l1 = 4
        if xmin >= 8650: l1 = 5

        if xmax <= 3876: l2 = 0
        if xmax <= 4990 and xmax >= 3754: l2 = 1
        if xmax <= 5667.9 and xmax >= 4616: l2 = 2
        if xmax <= 6693.9 and xmax >= 5668: l2 = 3
        if xmax <= 8649.9 and xmax >= 6694: l2 = 4
        if xmax >= 8650: l2 = 5

        if l1 == l2: arm_name = [arms[l1]]
        if l1 != l2: arm_name = arms[l1:l2+1]

        # choose the highest S/N
        source = []
        fits_name = []
        for it in range(len(arm_name)):

            # determine the folder location
            loc = getSpectrum.EDIBLES_DR3 + self.target_name + '/' + arm_name[it]
            if arm_name[it] == 'REDL_564' or arm_name[it] == 'REDU_564':
                loc = getSpectrum.EDIBLES_DR3 + self.target_name + '/RED_564'
            if arm_name[it] == 'REDL_860' or arm_name[it] == 'REDU_860':
                loc = getSpectrum.EDIBLES_DR3 + self.target_name + '/RED_860'

            # determine the arm
            if arm_name[it] == 'BLUE_346' or arm_name[it] == 'BLUE_437': arm_letter = 'B'
            if arm_name[it] == 'REDL_564' or arm_name[it] == 'REDL_860': arm_letter = 'L'
            if arm_name[it] == 'REDU_564' or arm_name[it] == 'REDU_860': arm_letter = 'U'
            globName = loc + '/' + '%s*_%s.fits'%(self.target_name, arm_letter)

            # extract the fits name
            if os.path.isdir(loc) == True:

                # extract fits name based on defined observation date
                if self.obs_date is not None:
                    if arm_name[it] == 'BLUE_346': w_ar = 'w346'
                    if arm_name[it] == 'BLUE_437': w_ar = 'w437'
                    if arm_name[it] == 'REDL_564' or arm_name[it] == 'REDU_564': w_ar = 'w564'
                    if arm_name[it] == 'REDL_860' or arm_name[it] == 'REDU_860': w_ar = 'w860'

                    glob_list = glob.glob(loc + '/' + self.target_name + '_' + w_ar + \
                                        '*' + self.obs_date + '_' + arm_letter + '.fits')
                    if len(glob_list) > 0:
                        fits_name.append(glob_list[0])
                        source.append(loc)

                # extract fits name based on the higest S/N
                else:
                    source.append(loc)
                    n_science = 0
                    file_count = 0
                    for file in glob.glob(globName):
                        file_count = file_count + 1
                        n_new = float(file.split('_')[3][1:])
                        if n_new > n_science:
                            n_science = n_new
                            temp_name = file
                        if file_count == len(glob.glob(globName)) - 1:
                            fits_name.append(temp_name)


        # check target be available
        if len(source) == 0:
            print 'This object have not been observed yet!'
            sys.exit()

        return [fits_name, source, xmin, xmax]




    # =========================================================================
    #   LoadSpectrum -- Load the spectrum fits file and apply the barycentric
    #                   correction from the header
    # =========================================================================
    def LoadSpectrum(self):

        # fits information -> for EDIBLES
        if self.origin == None:
            target_info = self.findSpectrum()
            fits_name = target_info[0]
            location = target_info[1]
            xmin = target_info[2]
            xmax = target_info[3]
        # fits information -> for common files
        else:
            if self.origin[-1] != '/':
                location = [self.origin + '/']
            else:
                location = [self.origin]
            fits_name = [location[0] + self.target_name]
            # wavelength range
            if len(self.central_wave) == 1:
                xmin = self.central_wave[0] - 10
                xmax = self.central_wave[0] + 10
            if len(self.central_wave) == 2:
                xmin = self.central_wave[0]
                xmax = self.central_wave[1]
            # check valid file
            if os.path.isfile(fits_name[0]) != True:
                print 'This object have not been observed yet!'
                sys.exit()


        # read the wavelength solution from header
        for it in range(len(location)):
            data = pyfits.open(fits_name[it])
            flux = data[0].data
            hdr = pyfits.getheader(fits_name[it]).copy()
            CRVAL1 = hdr['CRVAL1']
            CDELT1 = hdr['CDELT1']
            wave_air = CRVAL1 + CDELT1 * np.arange(len(flux))
            self.wave = np.concatenate((self.wave ,wave_air))
            self.flux = np.concatenate((self.flux ,flux))
            self.wave_unit = 'angstrom'


        # wavelength range
        idx = [i for i, e in enumerate(self.wave) if e > xmin and e < xmax]
        w_new = []
        f_new = []
        for it in idx:
            w_new.append(self.wave[it])
            f_new.append(self.flux[it])
        self.wave = np.array(w_new)
        self.flux = np.array(f_new)


        # apply barycentric velocity correction
        if self.vel_ref == 'hel':
            barycentric_corr = hdr['HIERARCH ESO QC VRAD BARYCOR']
            self.wave = self.wave * (1.0 + barycentric_corr/299792.458)
            self.reference_frame = 'barycentric'
            print 'wavelength is barycentric corrected'


        # median normalization
        if self.normal == 'median':
            self.flux = self.flux / np.median(self.flux)


        # convert angstrom to velocity
        if self.vel_wave > 0:
            self.wave = ((self.wave - self.vel_wave) / self.vel_wave) * 299792.458






# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%  fitting  %%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class fitSpectrum:

    # =========================================================================
    #   initialization -- Define the required initial arguments for using
    #                     in the fitSpectrum class.
    # =========================================================================
    def __init__(self, frq_os, gamma, fit_flag=None):
        self.frq_os = frq_os
        self.gamma = gamma
        self.fit_flag = fit_flag



    # =========================================================================
    #   Gaussian model -- The absorption Gaussian function
    # =========================================================================
    def gaussian_function(self, x, cw, amp, sig):
        gauss = amp * np.exp(-0.5*((x - cw)/sig)**2)
        return gauss



    # =========================================================================
    #   multi Gaussian model -- The multi absorption Gaussian function
    # =========================================================================
    def multi_gaussian_function(self, p, x=None, y=None, err=None, fjac=None):
        model = np.zeros_like(x)
        y_lin = p[0] + p[1]*(x-np.median(x)) + p[2]*(x-np.median(x))**2
        parm = p[3:]
        for i in range(0, len(parm), 3):
            cw = parm[i]
            amp = parm[i+1]
            sig = parm[i+2]
            model = model + self.gaussian_function(x, cw, amp, sig)
        model = y_lin - model
        self.gaussian_model = model
        status = 0
        return [status, (y-model)/err]



    # =========================================================================
    #   Voigt model -- the Mathemathical version of Voigt function
    # =========================================================================
    def math_voigt(self, x, cw, amp, sig, gam):
        sigma = sig / np.sqrt(2 * np.log(2))
        tau = amp * (np.real( wofz( (x - cw + 1j*gam) / (sigma/np.sqrt(2)) ) ) / (sigma/np.sqrt(2*np.pi)))
        y = np.exp(-tau) - 1.0
        return y



    # =========================================================================
    #   multi Voigt model -- the astronomical version of multi Voigt function
    # =========================================================================
    def multi_math_voigt(self, p, x=None, y=None, err=None, fjac=None):
        model = np.zeros_like(x)
        parm = p[3:]
        midle = []
        for i in range(0, len(parm), 4):
            cw = parm[i]
            amp = parm[i+1]
            sig = parm[i+2]
            gam = parm[i+3]
            model = model + self.math_voigt(x, cw, amp, sig, gam)
            midle.append(cw)
        y_lin = p[0] + p[1] * (x - np.mean(midle)) + p[2] * (x - np.mean(midle))**2
        model = y_lin + model
        self.math_voigt_model = model
        status = 0
        return [status, (y-model)/err]



    # =========================================================================
    #   Voigt model -- the astronomical version of Voigt function
    # =========================================================================
    def astro_voigt(self, x, cw, z, amp, b_eff, log_N, v_frq, v_gamma):
        c  = 2.99792458e10
        sigma0 = 0.0263
        x = np.array(x)
        nu = c / (x * 1.e-8)
        nu0 = c / (cw * 1.e-8)
        delta_nu = nu - nu0 / (1.0 + z)
        delta_nu_D = (b_eff*1.e5) * nu / c
        prf = amp / ((np.pi**0.5) * delta_nu_D)
        sig = delta_nu / delta_nu_D
        gam = v_gamma / (4 * np.pi * delta_nu_D)
        voigt = prf * wofz(sig + 1j*gam).real
        tau = (10**log_N) * sigma0 * v_frq * voigt
        voigt_model = np.exp(-tau) - 1

        return voigt_model



    # =========================================================================
    #   multi Voigt model -- the astronomical version of multi Voigt function
    # =========================================================================
    def multi_astro_voigt(self, p, x=None, y=None, err=None, fjac=None):
        model = np.zeros_like(x)
        y_lin = p[0] + p[1]*(x-np.median(x)) + p[2]*(x-np.median(x))**2
        parm = p[3:]
        peak_num = 0
        for i in range(0, len(parm), 5):
            cw = parm[i]
            z = parm[i+1]
            amp = parm[i+2]
            b_eff = parm[i+3]
            log_N = parm[i+4]
            model = model + self.astro_voigt(x, cw, z, amp, b_eff, log_N, self.frq_os[peak_num], self.gamma[peak_num])
            peak_num = peak_num + 1
        model = y_lin + model
        self.astro_voigt_model = model
        status = 0
        return [status, (y-model)/err]



    # =========================================================================
    #   multi double Voigt model -- the astronomical version of double multi
    #                               Voigt function, e.g. used for Na I line
    # =========================================================================
    def multi_double_astro_voigt(self, p, x=None, y=None, err=None, fjac=None):
        model = np.zeros_like(x)
        y_lin = p[0] + p[1]*(x-np.median(x)) + p[2]*(x-np.median(x))**2
        parm = p[3:]
        for i in range(0, len(parm), 9):
            cw1 = parm[i]
            cw2 = parm[i+1]
            z1 = parm[i+2]
            z2 = parm[i+3]
            amp1 = parm[i+4]
            amp2 = parm[i+5]
            b_eff1 = parm[i+6]
            b_eff2 = parm[i+7]
            log_N = parm[i+8]
            model = model + self.astro_voigt(x, cw1, z1, amp1, b_eff1, log_N, self.frq_os[0], self.gamma[0]) + \
                            self.astro_voigt(x, cw2, z2, amp2, b_eff2, log_N, self.frq_os[1], self.gamma[1])
        model = y_lin + model
        self.double_astro_voigt_model = model
        status = 0
        return [status, (y-model)/err]





    # =========================================================================
    #   initial value -- the initial value for all fits including the
    #                    constraints and limits for all methods
    # =========================================================================
    def initial_value(self, cw, continuum, min_continuum, max_continuum):

        # 3rd order polynominal as continuum
        initial_value = np.array([continuum, 0.005, 0.0], dtype='double')

        # gaussian
        if self.fit_flag == 'gaussian':
            for itr in range(len(cw)):
                initial_value = np.append(initial_value, cw[itr])
                initial_value = np.append(initial_value, 0.5)
                initial_value = np.append(initial_value, 1.0)

        # math voigt
        if self.fit_flag == 'math_voigt':
            for itr in range(len(cw)):
                initial_value = np.append(initial_value, cw[itr])
                initial_value = np.append(initial_value, 0.1)
                initial_value = np.append(initial_value, 0.1)
                initial_value = np.append(initial_value, 0.1)

        # voigt
        if self.fit_flag == 'astro_voigt':
            for itr in range(len(cw)):
                initial_value = np.append(initial_value, cw[itr])
                initial_value = np.append(initial_value, 0.0)
                initial_value = np.append(initial_value, 1.0)
                initial_value = np.append(initial_value, 3.0)
                initial_value = np.append(initial_value, 12.0)

        # double voigt
        if self.fit_flag == 'double_voigt':
            for itr in range(len(cw)):
                initial_value = np.append(initial_value, cw[itr][0])
                initial_value = np.append(initial_value, cw[itr][1])
                initial_value = np.append(initial_value, 0.0)
                initial_value = np.append(initial_value, 0.0)
                initial_value = np.append(initial_value, 1.0)
                initial_value = np.append(initial_value, 1.0)
                initial_value = np.append(initial_value, 2.0)
                initial_value = np.append(initial_value, 2.0)
                initial_value = np.append(initial_value, 12.0)


        # the constraints for the initial values
        p0 = initial_value
        parinfo = [{'value':0., 'fixed':0, 'limited':[0,0], 'limits':[0.,0.]} for i in range(len(p0))]
        parinfo[0]['limited'][0] = 1
        parinfo[0]['limits'][0] = min_continuum
        parinfo[0]['limited'][1] = 1
        parinfo[0]['limits'][1] = max_continuum
        # parinfo[2]['fixed'] = 1
        for itr in range(len(cw)):

            # ----------
            # gaussian -
            # ----------
            if self.fit_flag == 'gaussian':
                # CW
                parinfo[3+itr*3]['limited'][0] = 1
                parinfo[3+itr*3]['limits'][0] = cw[itr] - 0.3
                parinfo[3+itr*3]['limited'][1] = 1
                parinfo[3+itr*3]['limits'][1] = cw[itr] + 0.3
                # amp
                parinfo[4+itr*3]['limited'][0] = 1
                parinfo[4+itr*3]['limits'][0] = 0
                # sigma
                parinfo[5+itr*3]['limited'][0] = 1
                parinfo[5+itr*3]['limits'][0] = 0
                parinfo[5+itr*3]['limited'][1] = 1
                parinfo[5+itr*3]['limits'][1] = 10

            # ------------
            # math voigt -
            # ------------
            if self.fit_flag == 'math_voigt':
                # CW
                parinfo[3+itr*4]['limited'][0] = 1
                parinfo[3+itr*4]['limits'][0] = cw[itr] - 0.05
                parinfo[3+itr*4]['limited'][1] = 1
                parinfo[3+itr*4]['limits'][1] = cw[itr] + 0.05
                # sigma
                parinfo[5+itr*4]['limited'][0] = 1
                parinfo[5+itr*4]['limits'][0] = 0
                # gamma
                parinfo[6+itr*4]['limited'][0] = 1
                parinfo[6+itr*4]['limits'][0] = 0

            # -------------
            # astro voigt -
            # -------------
            if self.fit_flag == 'astro_voigt':
                # CW
                parinfo[3+itr*5]['limited'][0] = 1
                parinfo[3+itr*5]['limits'][0] = cw[itr] - 0.03
                parinfo[3+itr*5]['limited'][1] = 1
                parinfo[3+itr*5]['limits'][1] = cw[itr] + 0.03
                # redshift
                parinfo[4+itr*5]['fixed'] = 1
                # amp
                parinfo[5+itr*5]['limited'][0] = 1
                parinfo[5+itr*5]['limits'][0] = 0
                # b_eff
                parinfo[6+itr*5]['limited'][0] = 1
                parinfo[6+itr*5]['limits'][0] = 1
                parinfo[6+itr*5]['limited'][1] = 1
                parinfo[6+itr*5]['limits'][1] = 10
                # logN
                parinfo[7+itr*5]['limited'][0] = 1
                parinfo[7+itr*5]['limits'][0] = 7
                parinfo[7+itr*5]['limited'][1] = 1
                parinfo[7+itr*5]['limits'][1] = 20


            # --------------------
            # double astro voigt -
            # --------------------
            if self.fit_flag == 'double_voigt':
                # CW
                parinfo[3+itr*9]['limited'][0] = 1
                parinfo[3+itr*9]['limits'][0] = cw[itr][0] - 0.1
                parinfo[3+itr*9]['limited'][1] = 1
                parinfo[3+itr*9]['limits'][1] = cw[itr][0] + 0.1
                parinfo[4+itr*9]['limited'][0] = 1
                parinfo[4+itr*9]['limits'][0] = cw[itr][1] - 0.1
                parinfo[4+itr*9]['limited'][1] = 1
                parinfo[4+itr*9]['limits'][1] = cw[itr][1] + 0.1
                # redshift
                parinfo[5+itr*9]['fixed'] = 1
                parinfo[6+itr*9]['fixed'] = 1
                # amp
                parinfo[7+itr*9]['limited'][0] = 1
                parinfo[7+itr*9]['limits'][0] = 0
                parinfo[8+itr*9]['limited'][0] = 1
                parinfo[8+itr*9]['limits'][0] = 0
                # b_eff
                parinfo[9+itr*9]['limited'][0] = 1
                parinfo[9+itr*9]['limits'][0] = 1
                parinfo[9+itr*9]['limited'][1] = 1
                parinfo[9+itr*9]['limits'][1] = 12
                parinfo[10+itr*9]['limited'][0] = 1
                parinfo[10+itr*9]['limits'][0] = 1
                parinfo[10+itr*9]['limited'][1] = 1
                parinfo[10+itr*9]['limits'][1] = 12
                # logN
                parinfo[11+itr*9]['limited'][0] = 1
                parinfo[11+itr*9]['limits'][0] = 7
                parinfo[11+itr*9]['limited'][1] = 1
                parinfo[11+itr*9]['limits'][1] = 20

        for itr in range(len(p0)): parinfo[itr]['value'] = p0[itr]
        self.parinfo = parinfo
        self.p0 = p0



    # =========================================================================
    #   fitting -- all fitting accomplishing applying Levenberg-Marquardt
    #              fitting procedure using mpfit object
    # =========================================================================
    def afit(self, wave, spec, cw, measurement=None):

        # left and right point of peak
        if len(cw) == 1:
            xcl = cw[0]
            xcr = cw[0]
        if len(cw) > 1:
            xcl = cw[0]
            xcr = cw[-1]
        continuum_lim = np.mean(spec) + (max(spec) - np.mean(spec))/2.0

        # -------------------------------------------
        # iterate over +-8A around the peak to 0.3A
        # to select the best region around the peak.
        # -------------------------------------------
        temp_arr = np.zeros((2,26))
        for it in np.arange(8, 0.3, -0.3):
            idx = [i for i, e in enumerate(wave) if e > xcl-it and e < xcr+it]
            w_new = []
            f_new = []
            for i in idx:
                w_new.append(wave[i])
                f_new.append(spec[i])
            xn = np.array(w_new)
            yn = np.array(f_new)
            err = 0.05 * yn


            self.initial_value(cw, continuum_lim, continuum_lim-0.4, continuum_lim+0.4)
            fa = {'x':xn, 'y':yn, 'err':err}
            # gaussian
            if self.fit_flag == 'gaussian':
                m = mpfit.mpfit(self.multi_gaussian_function, self.p0, functkw=fa, quiet=1, parinfo=self.parinfo)
                chisq = (self.multi_gaussian_function(m.params, x=xn, y=yn, err=err)[1]**2).sum()
            # math voigt
            if self.fit_flag == 'math_voigt':
                m = mpfit.mpfit(self.multi_math_voigt, self.p0, functkw=fa, quiet=1, parinfo=self.parinfo)
                chisq = (self.multi_math_voigt(m.params, x=xn, y=yn, err=err)[1]**2).sum()
            # voigt
            if self.fit_flag == 'astro_voigt':
                m = mpfit.mpfit(self.multi_astro_voigt, self.p0, functkw=fa, quiet=1, parinfo=self.parinfo)
                chisq = (self.multi_astro_voigt(m.params, x=xn, y=yn, err=err)[1]**2).sum()

            # select the best regions based on the fit quality
            DOF = len(xn) - len(m.params)
            if chisq < chi2.isf(0.05, DOF):
                temp_arr[0,it] = chisq
                temp_arr[1,it] = it

        # remove zero iterations
        idx = np.nonzero(temp_arr[0,:] != 0)
        if idx[0].size > 0:
            col_chi = temp_arr[0,idx][0]
            col_lim = temp_arr[1,idx][0]
            lim_value = col_lim[np.nonzero(col_chi == min(col_chi))]
            if (lim_value == col_lim[0]) & (len(col_chi) > 1):
                lim_value = col_lim[np.argsort(col_chi)[-1]]
        else:
            lim_value = 0.45
            print 'CAUTION!!'
            print 'The procedure did not find any suitable region around'
            print 'the peak for fitting. This may be caused by very weak'
            print 'features or by deformed profiles, which caused a poor'
            print 'fitting. Maybe by changing the initial values, you could'
            print 'get a better fitting'


        # apply the determined wavelength limit
        idx = [i for i, e in enumerate(wave) if e > xcl-lim_value and e < xcr+lim_value]
        w_new = []
        f_new = []
        for i in idx:
            w_new.append(wave[i])
            f_new.append(spec[i])
        xn = np.array(w_new)
        yn = np.array(f_new)
        err = 0.05 * yn


        # -------------------------------------------
        # determine continuum level -- by iterate
        # over +-0.4A around the guessed continuum
        # -------------------------------------------
        cnt_value = np.zeros((2,10))
        m1 = continuum_lim - 0.4
        m2 = continuum_lim + 0.4
        cnt_temp = np.arange(m1, m2, 0.08)
        for it in range(10):
            self.initial_value(cw, cnt_temp[it], cnt_temp[it]-0.15, cnt_temp[it]+0.15)
            fa = {'x':xn, 'y':yn, 'err':err}
            # gaussian
            if self.fit_flag == 'gaussian':
                m = mpfit.mpfit(self.multi_gaussian_function, self.p0, functkw=fa, quiet=1, parinfo=self.parinfo)
                cnt_value[0,it] = (self.multi_gaussian_function(m.params, x=xn, y=yn, err=err)[1]**2).sum()
            # math voigt
            if self.fit_flag == 'math_voigt':
                m = mpfit.mpfit(self.multi_math_voigt, self.p0, functkw=fa, quiet=1, parinfo=self.parinfo)
                cnt_value[0,it] = (self.multi_math_voigt(m.params, x=xn, y=yn, err=err)[1]**2).sum()
            # voigt
            if self.fit_flag == 'astro_voigt':
                m = mpfit.mpfit(self.multi_astro_voigt, self.p0, functkw=fa, quiet=1, parinfo=self.parinfo)
                cnt_value[0,it] = (self.multi_astro_voigt(m.params, x=xn, y=yn, err=err)[1]**2).sum()

            cnt_value[1,it] = cnt_temp[it]

        max_value_f = cnt_value[1,np.nonzero(cnt_value[0,] == min(cnt_value[0,]))]
        max_value_f = max_value_f[0][0]
        p = m.params



        # -----------
        # final fit -
        # -----------
        if measurement is None:
            fa = {'x':xn, 'y':yn, 'err':err}
            self.initial_value(cw, max_value_f, max_value_f-0.15, max_value_f+0.15)
            # gaussian
            if self.fit_flag == 'gaussian':
                m = mpfit.mpfit(self.multi_gaussian_function, self.p0, functkw=fa, quiet=1, parinfo=self.parinfo)
            # math voigt
            if self.fit_flag == 'math_voigt':
                m = mpfit.mpfit(self.multi_math_voigt, self.p0, functkw=fa, quiet=1, parinfo=self.parinfo)
            # voigt
            if self.fit_flag == 'astro_voigt':
                m = mpfit.mpfit(self.multi_astro_voigt, self.p0, functkw=fa, quiet=1, parinfo=self.parinfo)

            self.fit_params = m.params

            # plot
            p = m.params
            temp_y = np.full((1, len(xn)), 1)
            temp_err = 0.02 * temp_y
            if self.fit_flag == 'gaussian':
                self.multi_gaussian_function(p, xn, temp_y, temp_err)
                yfit = self.gaussian_model
            if self.fit_flag == 'math_voigt':
                self.multi_math_voigt(p, xn, temp_y, temp_err)
                yfit = self.math_voigt_model
            if self.fit_flag == 'astro_voigt':
                self.multi_astro_voigt(p, xn, temp_y, temp_err)
                yfit = self.astro_voigt_model

            plt.step(wave,flux,'gray')
            plt.plot(xn, yfit, 'magenta')
            plt.show()



        # --------------
        #  measurement -
        # --------------
        if measurement is not None:

            # sigma value of the continuum region
            wave = np.array(wave)
            idx = np.nonzero((wave >= xcl-lim_value-2.0) & (wave <= xcl-lim_value))
            region_left = spec[idx]
            idx = np.nonzero((wave >= xcr+lim_value) & (wave <= xcr+lim_value+2))
            region_right = spec[idx]
            sigma = min([np.std(region_left), np.std(region_right)])


            y = yn
            if self.fit_flag == 'gaussian': parameters = np.zeros(shape=(3+len(cw)*3,100))
            if self.fit_flag == 'math_voigt': parameters = np.zeros(shape=(3+len(cw)*4,100))
            if self.fit_flag == 'astro_voigt': parameters = np.zeros(shape=(3+len(cw)*5,100))
            self.initial_value(cw, p[0], p[0]-0.3, p[0]+0.3)
            for it in range(100):

                # add noise to spectrum
                noise = np.random.normal(0, sigma/2.0, len(y))
                yn = y + noise
                err = 0.02 * yn

                # fit
                fa = {'x':xn, 'y':yn, 'err':err}

                if self.fit_flag == 'gaussian':
                    m = mpfit.mpfit(self.multi_gaussian_function, self.p0, functkw=fa, quiet=1, parinfo=self.parinfo)
                    parameters[0,it] = m.params[0]
                    parameters[1,it] = m.params[1]
                    parameters[2,it] = m.params[2]
                    for itr in range(len(cw)):
                        parameters[3+itr*3,it] = m.params[3+itr*3]
                        parameters[4+itr*3,it] = m.params[4+itr*3]
                        parameters[5+itr*3,it] = m.params[5+itr*3]


                if self.fit_flag == 'math_voigt':
                    m = mpfit.mpfit(self.multi_math_voigt, self.p0, functkw=fa, quiet=1, parinfo=self.parinfo)
                    parameters[0,it] = m.params[0]
                    parameters[1,it] = m.params[1]
                    parameters[2,it] = m.params[2]
                    for itr in range(len(cw)):
                        parameters[3+itr*4,it] = m.params[3+itr*4]
                        parameters[4+itr*4,it] = m.params[4+itr*4]
                        parameters[5+itr*4,it] = m.params[5+itr*4]
                        parameters[6+itr*4,it] = m.params[6+itr*4]


                if self.fit_flag == 'astro_voigt':
                    m = mpfit.mpfit(self.multi_astro_voigt, self.p0, functkw=fa, quiet=1, parinfo=self.parinfo)
                    # DOF = len(xn)-len(m.params)
                    # print (self.multi_astro_voigt(m.params, x=xn, y=yn, err=err)[1]**2).sum(), chi2.isf(0.05,DOF)
                    parameters[0,it] = m.params[0]
                    parameters[1,it] = m.params[1]
                    parameters[2,it] = m.params[2]
                    for itr in range(len(cw)):
                        parameters[3+itr*5,it] = m.params[3+itr*5]
                        parameters[4+itr*5,it] = m.params[4+itr*5]
                        parameters[5+itr*5,it] = m.params[5+itr*5]
                        parameters[6+itr*5,it] = m.params[6+itr*5]
                        parameters[7+itr*5,it] = m.params[7+itr*5]



            # median all parameters
            if self.fit_flag == 'gaussian': fited_params = np.zeros(shape=(3+len(cw)*3,2))
            if self.fit_flag == 'math_voigt': fited_params = np.zeros(shape=(3+len(cw)*4,2))
            if self.fit_flag == 'astro_voigt': fited_params = np.zeros(shape=(3+len(cw)*5,2))

            for it in range(len(parameters)):
                fited_params[it, 0] = np.median(parameters[it,:])
                fited_params[it, 1] = np.std(parameters[it,:])/np.sqrt(len(parameters[it,:]))

            self.fit_params = fited_params


            # plotting
            if self.fit_flag == 'gaussian':
                parms = np.zeros(shape=(6, len(cw)))
                for itr in range(len(cw)):
                    parms[0, itr] = fited_params[0, 0]
                    parms[1, itr] = fited_params[1, 0]
                    parms[2, itr] = fited_params[2, 0]
                    parms[3, itr] = fited_params[3+itr*3, 0]
                    parms[4, itr] = fited_params[4+itr*3, 0]
                    parms[5, itr] = fited_params[5+itr*3, 0]

            if self.fit_flag == 'math_voigt':
                parms = np.zeros(shape=(7, len(cw)))
                for itr in range(len(cw)):
                    parms[0, itr] = fited_params[0, 0]
                    parms[1, itr] = fited_params[1, 0]
                    parms[2, itr] = fited_params[2, 0]
                    parms[3, itr] = fited_params[3+itr*4, 0]
                    parms[4, itr] = fited_params[4+itr*4, 0]
                    parms[5, itr] = fited_params[5+itr*4, 0]
                    parms[6, itr] = fited_params[6+itr*4, 0]

            if self.fit_flag == 'astro_voigt':
                parms = np.zeros(shape=(8, len(cw)))
                for itr in range(len(cw)):
                    parms[0, itr] = fited_params[0, 0]
                    parms[1, itr] = fited_params[1, 0]
                    parms[2, itr] = fited_params[2, 0]
                    parms[3, itr] = fited_params[3+itr*5, 0]
                    parms[4, itr] = fited_params[4+itr*5, 0]
                    parms[5, itr] = fited_params[5+itr*5, 0]
                    parms[6, itr] = fited_params[6+itr*5, 0]
                    parms[7, itr] = fited_params[7+itr*5, 0]


            plt.plot(wave, flux, 'gray')

            temp_y = np.full((1, len(xn)), 1)
            temp_err = 0.02*temp_y
            if self.fit_flag == 'gaussian':
                self.multi_gaussian_function(fited_params[:,0], xn, temp_y, temp_err)
                yfit = self.gaussian_model
            if self.fit_flag == 'math_voigt':
                self.multi_math_voigt(fited_params[:,0], xn, temp_y, temp_err)
                yfit = self.math_voigt_model
            if self.fit_flag == 'astro_voigt':
                self.multi_astro_voigt(fited_params[:,0], xn, temp_y, temp_err)
                yfit = self.astro_voigt_model

            plt.plot(xn, yfit, color='blue')
            for itr in range(len(cw)):
                if self.fit_flag == 'gaussian':
                    self.multi_gaussian_function(parms[:,itr], xn, temp_y, temp_err)
                    yfit = self.gaussian_model
                if self.fit_flag == 'math_voigt':
                    self.multi_math_voigt(parms[:,itr], xn, temp_y, temp_err)
                    yfit = self.math_voigt_model
                if self.fit_flag == 'astro_voigt':
                    self.multi_astro_voigt(parms[:,itr], xn, temp_y, temp_err)
                    yfit = self.astro_voigt_model
                # sub-features
                plt.plot(xn, yfit, color='#c277dd')
            # continuum line
            plt.plot(xn,parms[0,0]+parms[1,0]*(xn-np.median(xn))+parms[2,0]*(xn-np.median(xn))**2,'#e035cf',linestyle='--')

            plt.show()
