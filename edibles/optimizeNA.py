import numpy as np
import astropy.constants as cst
from lmfit import Model
from lmfit.models import update_param_vals
from edibles.utils.voigt_profile import voigt_absorption_line
from scipy.interpolate import interp1d
from scipy.stats import pearsonr
import inspect
import collections

class NAModelOptimize(Model):
    """A model of the astronomical Voigt function."""

    def __init__(self, n_components, independent_vars=['x'], prefix='', nan_policy='raise',
                 **kwargs):
        if not isinstance(n_components, int):
            raise TypeError("n_components (%.1f) must be an integer!" % n_components)
        self.n_components = n_components
        kwargs.update({'prefix': prefix, 'nan_policy': nan_policy,
                       'independent_vars': independent_vars})

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

        def calcMultiNa(x, b_Cloud0=1.0, N_Cloud0=1.0, V_off_Cloud0=0.0, **kwargs):
            lambda0 = [3302.369, 3302.978, 5889.950, 5895.924] * self.n_components
            f = [8.26e-03, 4.06e-03, 6.49E-01, 3.24E-01] * self.n_components
            gamma = [6.280e7, 6.280e7, 6.280e7, 6.280e7] * self.n_components
            N_mag = 14.0
            v_resolution = 3.0

            # parse parameters
            bs = [b_Cloud0]*4
            Ns = [N_Cloud0 * 10**N_mag]*4
            V_offs = [V_off_Cloud0]*4
            # just something to print so we know it's working...
            print("Cloud velocities at:", np.unique(V_offs))

            for name in kwargs.keys():
                if name[0] == "b":
                    print(name) # will be popped up as Cloud idx??
                    bs = bs + [kwargs[name]]*4
                if name[0] == "N":
                    Ns = Ns + [kwargs[name] * 10**N_mag]*4
                if name[0] == "V":
                    V_offs = V_offs + [kwargs[name]]*4

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
            )

            return flux

        sig = inspect.signature(calcMultiNa)
        base_b = inspect.signature(calcMultiNa).parameters["b_Cloud0"]
        base_N = inspect.signature(calcMultiNa).parameters["N_Cloud0"]
        base_V_off = inspect.signature(calcMultiNa).parameters["V_off_Cloud0"]

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
        calcMultiNa.__signature__ = sig.replace(parameters=tuple(d.values()))

        super().__init__(calcMultiNa, **kwargs)

    def guess(self, data,  wavegrid=None, V_off=[0.0], **kwargs):
        V_off[-1] = self.calc_vrad(wavegrid, data)

        pars = self.make_params()
        for i, v in enumerate(V_off):
            pars["%sb_Cloud%i" % (self.prefix, i)].set(value=1.0, min=0, max=10)
            pars["%sN_Cloud%i" % (self.prefix, i)].set(value=1.0, min=0, max=500)
            pars["%sV_off_Cloud%i" % (self.prefix, i)].set(value=v, min=-40, max=40)
            # if we have better estimate on v, consider using v pm 15 as min and max

        return update_param_vals(pars, self.prefix, **kwargs)

    def calc_vrad(self,wave,flux):

        lambda0 = [3302.369, 3302.978, 5889.950, 5895.924]
        f = [8.26e-03, 4.06e-03, 6.49E-01,3.24E-01]
        gamma = [6.280e7, 6.280e7,6.280e7,6.280e7]
        v_resolution = 5.75
        b = 1.0
        N = 1.8e12
        v_rad = 0
        AbsorptionLine = voigt_absorption_line(
            wave,
            lambda0=lambda0,
            b=b,
            N=N,
            f=f,
            gamma=gamma,
            v_rad=v_rad,
            v_resolution=v_resolution,
        )
        bestV = None
        bestLine = None
        #1. calculate model v_rad = 0 (wave, flux)
        #2. calculate shifted model (wavev2, flux)
        # (wave2 -wave)/wave = v_rad/c    [-10,-9.5,-9,.....]  c = constant(look in astropy)
        # wave1 *(1+vrad/c)
        #3 for each vrad calculate correlation between model and data
        v = -40
        while v <= 40:
            mult =  1.0 + v/cst.c.to("km/s").value
            testW = (np.array(wave) * mult)

            interpolationfunction = interp1d(
                        testW, AbsorptionLine, kind="cubic", fill_value="extrapolate")
            interpolatedModel = interpolationfunction(wave)
            new, _ = pearsonr(flux, interpolatedModel)
            if not bestV:
                bestV = v
                bestLine = new
            elif new > bestLine:
                bestV = v
                bestLine = new
            v += 0.5

        return bestV
