# Check how to access individual parameter and components in a CompositeModel of LMFIT

import numpy as np
import inspect
import collections
from scipy.interpolate import CubicSpline
from lmfit import Model
from lmfit.models import update_param_vals
from edibles.utils.voigt_profile import voigt_absorption_line
from edibles.models import ContinuumModel


class multiNaModel(Model):
    def __init__(self, n_components, independent_vars=["x"],
                 prefix="",
                 nan_policy="raise",
                 n_step=25,
                 **kwargs):

        if not isinstance(n_components, int):
            raise TypeError("n_components (%.1f) must be an integer!" % n_components)
        self.n_components = n_components
        self.n_setp = n_step

        kwargs.update({"prefix": prefix,
                       "nan_policy": nan_policy,
                       "independent_vars": independent_vars})

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
                n_step=self.n_setp
            )

            return flux

        # play with the signature of the calculation function
        # I don't really understand what is happening...
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

    def guess(self, V_off=[0.0], **kwargs):
        # For now just type in V_off but we can include v_correlation in the future

        assert len(V_off) == self.n_components, "Number of components do not match."

        pars = self.make_params()
        for i, v in enumerate(V_off):
            pars["%sb_Cloud%i" % (self.prefix, i)].set(value=1.0, min=0, max=10)
            pars["%sN_Cloud%i" % (self.prefix, i)].set(value=1.0, min=0, max=1000)
            pars["%sV_off_Cloud%i" % (self.prefix, i)].set(value=v, min=-50, max=50)
            # if we have better estimate on v, consider using v pm 15 as min and max

        return update_param_vals(pars, self.prefix, **kwargs)


if __name__ == "__main__":
    from edibles.utils.edibles_oracle import EdiblesOracle
    from edibles.utils.edibles_spectrum import EdiblesSpectrum
    import matplotlib.pyplot as plt
    import time

    pythia = EdiblesOracle()

    # data for 3300 segment
    List = pythia.getFilteredObsList(object=["HD 183143"], OrdersOnly=True, Wave=3302.0)
    test = List.values.tolist()
    filename = test[1]
    wrange = [3301.5, 3304.0]
    sp = EdiblesSpectrum(filename)
    wave, flux = sp.bary_wave, sp.flux
    idx = np.where((wave > wrange[0]) & (wave < wrange[1]))
    wave, flux = wave[idx], flux[idx]

    print("=" * 20)
    print("Fitting....")

    continuum_model = ContinuumModel(n_anchors=5)
    pars_guess = continuum_model.guess(flux, wave)

    line_model = multiNaModel(n_components=2, n_step=12)
    pars_guess.update(line_model.guess(V_off=[-12, 0]))

    model2fit = continuum_model * line_model

    result = model2fit.fit(data=flux, params=pars_guess, x=wave)
    print("="*30)
    # play with parameters
    print(result.params)
    print(result.params["V_off_Cloud1"].value)
    print(result.params["V_off_Cloud1"].value+1.0)

    plt.plot(wave, flux)
    plt.plot(wave, result.best_fit)   # use stored result
    plt.plot(wave, model2fit.eval(params=result.params, x=wave))  # or to calculate it
    comps = result.eval_components(x=wave)    # get the components in the model
    plt.plot(wave, comps["cont"])
    plt.show()

    plt.plot(wave, comps["calcMultiNa"])
    plt.show()

