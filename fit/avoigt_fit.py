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
import edibles.fit.avoigt as avoigt
import edibles.fit.ref_index as ref_index
warnings.simplefilter('ignore')


class fitSpectrum:

    # =========================================================================
    #   initialization -- Define the required initial arguments for using
    #                     in the fitSpectrum class.
    # =========================================================================
    def __init__(self):
        pass


    # =========================================================================
    #   multi voigt transision -- fit multiple transition simultaneously
    # =========================================================================
    def voigt_func(self, p, x=None, y=None, err=None, fjac=None):

        # define the overall params
        model = np.zeros_like(x)
        flg = self.flg

        # --------------------
        # simple voigt fitting
        # --------------------
        if flg == 'voigt_fit':
            # continuum
            if self.cheb_order == -1:
                cont_model = np.full(len(x), 1.0)
            else:
                cont_model = p[0] + p[1] * (x - np.mean(x)) + p[2] * (2*(x - np.mean(x))**2 - 1.0)
            parm = p[3:]
            peak_num = 0
            for loop_p in range(0, len(parm), 5):
                vel_cloud   = parm[loop_p]
                lambda_peak = parm[loop_p+1]
                gamma       = parm[loop_p+2]
                b_eff       = parm[loop_p+3]
                log_N       = parm[loop_p+4]
                model = model + avoigt.voigt(x, lambda_peak=lambda_peak, b_eff=b_eff,  # vel_cloud=vel_cloud,
                log_N=log_N, gamma=gamma, osc_freq=self.frq_os[peak_num], resolving_power=self.resolving_power)
                peak_num = peak_num + 1



        # ----------------
        # multi components
        # ----------------
        if flg == 'multi_voigt':
            # continuum
            if self.cheb_order == -1:
                cont_model = np.full(len(x), 1.0)
            else:
                cont_model = p[0] + p[1] * (x - np.mean(x)) + p[2] * (2*(x - np.mean(x))**2 - 1.0)
            parm = p[3:]
            Nc = len(parm)/float(self.Nd)
            for loop_m in range(self.Nd):
                sub_parm = parm[int(loop_m*Nc):int((loop_m+1)*Nc)]
                atom_count = 0
                for loop_s in range(0, len(sub_parm)-3, 2):
                    lambda_peak   = sub_parm[loop_s+0]
                    gamma         = sub_parm[loop_s+1]
                    vel_cloud     = sub_parm[-3]
                    b_eff         = sub_parm[-2]
                    log_N         = sub_parm[-1]
                    model = model + avoigt.voigt(x, lambda_peak=lambda_peak, b_eff=b_eff,  # vel_cloud=vel_cloud,
                    log_N=log_N, gamma=gamma, osc_freq=self.frq_os[atom_count], resolving_power=self.resolving_power)
                    atom_count += 1


        model = cont_model + model
        self.voigt_model = model
        status = 0
        return [status, (y-model)/err]



    # =========================================================================
    #   fitting -- all fitting accomplishing applying Levenberg-Marquardt
    #              fitting procedure using mpfit object
    # =========================================================================
    def afit(self, wave, spec, cw, yerr=None, lines=None, cheb_order=None, resolving_power=None):

        # convert to numpy array
        wave = np.array(wave)
        spec = np.array(spec)

        # estimate uncertanties
        error = err.errorEst(wave, spec)

        # determine the chebyshev order
        if cheb_order is None: self.cheb_order = -1
        else: self.cheb_order = cheb_order

        # check the instrumental resolving power
        if resolving_power is None: self.resolving_power = 80000
        else: self.resolving_power = resolving_power

        # determine the components of cloud
        if len(cw) == 1: self.flg = 'voigt_fit'
        else:
            self.flg = 'multi_voigt'
            try:
                if len(cw[0]) >= 1: self.Nd = len(cw)
            except TypeError: self.Nd = 1

        # estimate position of continuum
        cont_guess = cont_est.cont_est(wave, spec, cw)

        # initial guess
        if lines is not None:
            l0, self.frq_os, damp_coef = line_properties.line_properties(cw, lines)
            l0 = [x / 10.0 for x in l0]
            # convert vac to air
            line_lam_0 = []
            for loop_line in range(len(l0)):
                line_lam_0.append(10.0*ref_index.vac2air(l0[loop_line], t=0))
            parinfo, p0 = ini_guess.initial_value(cw, cont_guess, damp_coef=damp_coef, lambda_zero=line_lam_0)
        else:
            parinfo, p0 = ini_guess.initial_value(cw, cont_guess)

        # fitting
        fa = {'x':wave, 'y':spec, 'err':error}
        m = mpfit.mpfit(self.voigt_func, p0, functkw=fa, quiet=1, parinfo=parinfo)

        # return the fitting values
        self.fit_parm = m.params
        self.fit_err  = m.perror
        self.fit_norm = m.fnorm
        self.yfit = self.voigt_model
