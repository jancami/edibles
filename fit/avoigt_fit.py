import os, glob, sys
import numpy as np
import matplotlib.pyplot as plt
import mpfit
import ref_index
from scipy.special import wofz
from scipy.stats import chi2
from scipy import interpolate, signal, stats
import warnings
import edibles.fit.initial_suggest as inig
import edibles.fit.err_est as err
import edibles.fit.mpfit as mpfit
import edibles.fit.cont_est as cont_est
import edibles.fit.line_properties as line_properties
import edibles.fit.avoigt as avoigt
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

        # continuum
        # if flg == 'voigt_fit':
        #     if self.cheb_order > 0:
        #         # chebyshev order 2
        #         cont_model = p[0] + p[1] * (x - np.mean(x)) + p[2] * (2*(x - np.mean(x))**2 - 1.0)
        #         parm = p[3:]
        #     else:
        #         cont_model = np.full(len(x), 1.0)
        #         parm = p[3:]
        # else:
        #     cont_model = np.full(len(x), 1.0)
        #     parm = p[1:]


        # print flg
        # plt.plot(x,y)
        # plt.plot(x,cont_model,'magenta')
        # plt.show()
        # print cont_model


        # simple voigt fitting
        if flg == 'voigt_fit':
            # continuum
            if self.cheb_order == -1: cont_model = np.full(len(x), 1.0)
            else: cont_model = p[0] + p[1] * (x - np.mean(x)) + p[2] * (2*(x - np.mean(x))**2 - 1.0)
            parm = p[3:]
            peak_num = 0
            for loop_p in range(0, len(parm), 5):
                vel_cloud = parm[loop_p]
                lambda0   = parm[loop_p+1]
                gamma     = parm[loop_p+2]
                b_eff     = parm[loop_p+3]
                log_N     = parm[loop_p+4]
                model = model + avoigt.voigt(x, vel_cloud=vel_cloud,
                                lambda0=lambda0, b_eff=b_eff,
                                log_N=log_N, gamma=gamma,
                                osc_freq=self.frq_os[peak_num])
                peak_num = peak_num + 1


        # multi transition for a single component
        if flg == 'single_cloud':
            atom_count = 0
            for loop_s in range(0, len(parm)-2, 3):
                vel_cloud = parm[loop_s]
                lambda0   = parm[loop_s+1]
                gamma     = parm[loop_s+2]
                b_eff     = parm[-2]
                log_N     = parm[-1]
                model = model + avoigt.voigt(x, vel_cloud=vel_cloud, lambda0=lambda0, b_eff=b_eff,
                                log_N=log_N, gamma=gamma, osc_freq=self.frq_os[atom_count])
                atom_count += 1


        # multi components
        if flg == 'multi_cloud':
            # continuum
            if self.cheb_order == -1: cont_model = np.full(len(x), 1.0)
            else: cont_model = p[0] + p[1] * (x - np.mean(x)) + p[2] * (2*(x - np.mean(x))**2 - 1.0)
            parm = p[3:]
            Nc = len(parm)/float(self.Nd)
            for loop_m in range(self.Nd):
                sub_parm = parm[int(loop_m*Nc):int((loop_m+1)*Nc)]
                atom_count = 0
                for loop_s in range(0, len(sub_parm)-2, 2):
                    lambda0   = sub_parm[loop_s+0]
                    gamma     = sub_parm[loop_s+1]
                    vel_cloud = sub_parm[-3]
                    b_eff     = sub_parm[-2]
                    log_N     = sub_parm[-1]
                    model = model + avoigt.voigt(x, vel_cloud=vel_cloud, lambda0=lambda0, b_eff=b_eff,
                                    log_N=log_N, gamma=gamma, osc_freq=self.frq_os[atom_count])
                    atom_count += 1

        model = cont_model + model
        self.voigt_model = model
        status = 0
        return [status, (y-model)/err]






    # =========================================================================
    #   fitting -- all fitting accomplishing applying Levenberg-Marquardt
    #              fitting procedure using mpfit object
    # =========================================================================
    def afit(self, wave, spec, cw, yerr=None, lines=None, cheb_order=None):

        # convert to numpy array
        wave = np.array(wave)
        spec = np.array(spec)

        # estimate uncertanties
        error = err.errorEst(wave, spec)

        # determine the chebyshev order
        if cheb_order is None: self.cheb_order = -1
        else: self.cheb_order = cheb_order


        # determine the components of cloud
        if len(cw) == 1: self.flg = 'voigt_fit'
        else:
            try:
                if len(cw[0]) >= 1:
                    self.flg = 'multi_cloud'
                    self.Nd = len(cw)
            except TypeError:
                self.flg = 'multi_cloud'#'single_cloud'
                # cw_temp = []
                # cw_temp.append(cw)
                # cw = cw_temp
                self.Nd = 1



        # estimate position of continuum
        cont_guess = cont_est.cont_est(wave, spec, cw)
        # print cont_guess


        # initial guess
        if lines is not None:
            l0, self.frq_os, damp_coef = line_properties.line_properties(cw, lines)
            parinfo, p0 = inig.initial_value(cw, cont_guess, damp_coef)
        else:
            parinfo, p0 = inig.initial_value(cw, cont_guess)


        # print p0
        # print parinfo
        # fitting
        fa = {'x':wave, 'y':spec, 'err':error}
        m = mpfit.mpfit(self.voigt_func, p0, functkw=fa, quiet=1, parinfo=parinfo)

        self.fit_parm = m.params
        self.fit_err  = m.perror

        # print m.params
        #
        self.yfit = self.voigt_model
        # plt.step(wave, spec, 'b.')
        # plt.plot(wave, yfit, 'magenta')
        # plt.show()
        # print m.perror
