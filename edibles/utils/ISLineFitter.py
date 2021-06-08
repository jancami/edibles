import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import inspect
import collections
from math import floor
from lmfit import Model
from lmfit.models import update_param_vals
from edibles.utils.voigt_profile import voigt_absorption_line

class ISLineFitter():
    def __init__(self, wave, flux):
        # wave and flux from edibles spectrum or elsewhere
        self.wave = wave
        self.flux = flux

    def fit(self):
        # this will do the fitting
        pass


class ISLineModel(Model):
    def __init__(self, n_components,
                 lam_0 = [3302.369, 3302.978],
                 fjj = [8.26e-03, 4.06e-03],
                 gamma = [6.280e7, 6.280e7],
                 Nmag = 14,
                 v_res = 3.0,
                 independent_vars=["x"],
                 prefix="",
                 nan_policy="raise",
                 n_step=25,
                 **kwargs):
        """
        :param n_components: int, number of velocity components
        :param lam_0: list, air wavelength of target line. lam_0, fjj and gamma should have same length
        :param fjj: list, oscillator strengths
        :param gamma: list, Gamma parameter related to broadening.
        :param Nmag: int, magnitude of column density
        :param v_res: float, resolution in km/s
        :param independent_vars: from lmfit and Klay's code
        :param prefix: from lmfit and Klay's code
        :param nan_policy: from lmfit and Klay's code
        :param n_step: int, no. of points in 1*FWHM during calculation. Under-sample losses information
        but over-sample losses efficiency.
        :param kwargs: ???
        """
        self.n_components, self.lam_0, self.fjj, self.gamma, self.Nmag, self.n_setp = \
            self.__inputCheck(n_components, lam_0, fjj, gamma, Nmag, n_step)
        self.v_res = v_res

        kwargs.update({"prefix": prefix,
                       "nan_policy": nan_policy,
                       "independent_vars": independent_vars})

        # creat b, N, and V according to no. of components
        self.b_names = ["b_Cloud%i" % (i) for i in range(n_components)]
        self.N_names = ["N_Cloud%i" % (i) for i in range(n_components)]
        self.V_names = ["V_off_Cloud%i" % (i) for i in range(n_components)]
        kwargs["param_names"] = self.b_names + self.N_names + self.V_names
        params = {}
        for name in kwargs["param_names"]:
            if name[0] == "b":
                params[name] = 1.0
            if name[0] == "N":
                params[name] = 1.0
            if name[0] == "V":
                params[name] = 0.0

        # Other than the default parameters from Cloud0, other parameters are passed in kwargs
        def calcISLineModel(x, b_Cloud0=1.0, N_Cloud0=1.0, V_off_Cloud0=0.0, **kwargs):
            lambda0 = self.lam_0 * self.n_components
            f = self.fjj * self.n_components
            gamma = self.gamma * self.n_components
            N_mag = self.Nmag
            v_resolution = self.v_res

            # parse parameters
            bs = [b_Cloud0] * len(self.lam_0)
            Ns = [N_Cloud0 * 10 ** N_mag] * len(self.lam_0)
            V_offs = [V_off_Cloud0] * len(self.lam_0)

            # just something to print so we know it's working...
            print("Velocity of Cloud_0:", np.unique(V_offs))

            for name in kwargs.keys():
                if name[0] == "b":
                    bs = bs + [kwargs[name]] * len(self.lam_0)
                if name[0] == "N":
                    Ns = Ns + [kwargs[name] * 10 ** N_mag] * len(self.lam_0)
                if name[0] == "V":
                    V_offs = V_offs + [kwargs[name]] * len(self.lam_0)

            # no problem for convolution since we only call voigt_absorption_line once.
            flux = voigt_absorption_line(
                x,
                lambda0=lambda0,
                b=bs,
                N=Ns,
                f=f,
                gamma=gamma,
                v_rad=V_offs,
                v_resolution=v_resolution,
                n_step=self.n_setp
            )

            return flux

        # play with the signature of the calculation function
        # I don't really understand what is happening...
        sig = inspect.signature(calcISLineModel)
        base_b = inspect.signature(calcISLineModel).parameters["b_Cloud0"]
        base_N = inspect.signature(calcISLineModel).parameters["N_Cloud0"]
        base_V_off = inspect.signature(calcISLineModel).parameters["V_off_Cloud0"]

        d = {'x': sig.parameters['x']}
        for i in range(n_components):
            b_key = "b_Cloud" + str(i)
            b_val = base_b.replace(name=b_key)
            d[b_key] = b_val

            N_key = "N_Cloud" + str(i)
            N_val = base_N.replace(name=N_key)
            d[N_key] = N_val

            V_off_key = "V_off_Cloud" + str(i)
            V_off_val = base_V_off.replace(name=V_off_key)
            d[V_off_key] = V_off_val

        d = collections.OrderedDict(d)
        calcISLineModel.__signature__ = sig.replace(parameters=tuple(d.values()))

        super().__init__(calcISLineModel, **kwargs)

    def guess(self, V_off=[0.0], **kwargs):
        # For now just type in V_off but we can include v_correlation in the future

        assert len(V_off) == self.n_components, "Number of components do not match."

        pars = self.make_params()
        for i, v in enumerate(V_off):
            pars["%sb_Cloud%i" % (self.prefix, i)].set(value=1.0, min=0, max=10)
            pars["%sN_Cloud%i" % (self.prefix, i)].set(value=1.0, min=0, max=1000)
            pars["%sV_off_Cloud%i" % (self.prefix, i)].set(value=v, min=v-10, max=v+10)
            # if we have better estimate on v, consider using v pm 15 as min and max

        return update_param_vals(pars, self.prefix, **kwargs)

    def __inputCheck(self, n_components, lam_0, fjj, gamma, Nmag, n_step):
        # n_components should be int
        if not isinstance(n_components, int):
            raise TypeError("n_components (%.1f) must be an integer!" % n_components)

        # lam_0, fjj, gamma should be LIST have same length
        if not isinstance(lam_0, list) or not isinstance(fjj, list) or not isinstance(gamma, list):
            raise TypeError("lam_0, fjj, gamma should be list")
        len_array = [len(lam_0), len(fjj), len(gamma)]
        if not np.max(len_array) == np.max(len_array):
            raise TypeError("lam_0, fjj, gamma should have the same length")

        # Nmag and n_step should be int but just in case it's a float
        Nmag = floor(Nmag)
        n_step = floor(n_step)
        return n_components, lam_0, fjj, gamma, Nmag, n_step



if __name__ == "__main__":
    print("Hello Word!")
