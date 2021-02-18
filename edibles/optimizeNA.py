import numpy as np
from lmfit import Model
from lmfit.models import update_param_vals
from edibles.utils.voigt_profile import voigt_absorption_line

def calc_vrad(wave,flux):

    lambda0 = [3302.369, 3302.978, 5889.950, 5895.924]
    f = [8.26e-03, 4.06e-03, 6.49E-01,3.24E-01]
    gamma = [6.280e7, 6.280e7,6.280e7,6.280e7]
    v_resolution = 5.75
    b = 1.1
    N = 18.9e13

    bestV = None
    bestLine = None
    v = -13
    while v <= -10:
        AbsorptionLine = voigt_absorption_line(
        wave,
        lambda0=lambda0,
        b=b,
        N=N,
        f=f,
        gamma=gamma,
        v_rad=v,
        v_resolution=v_resolution,
        )

        new = np.correlate(AbsorptionLine,flux)
        if not bestV:
            bestV = v
            bestLine = new
        elif new > bestLine:
            bestV = v
            bestLine = new
        v += 0.5

    return bestV


def guess_NA(model, data, v_rad, b, n, x):
    if data is None:
        raise ValueError('y does not exist')

    y = np.asarray(data)
    v_rad[-1] = calc_vrad(x,y)
    print ("vrad:", v_rad)


    voigt_pars = model.make_params(v_rad = v_rad, b = b, N = n)

    return voigt_pars


class NAModelOptimize(Model):
    """A model of the astronomical Voigt function."""


    def __init__(self, independent_vars=['wavegrid','lambda0','f','gamma','v_resolution'], prefix='', nan_policy='raise',
                 **kwargs):
        kwargs.update({'prefix': prefix, 'nan_policy': nan_policy,
                       'independent_vars': independent_vars})
        super().__init__(voigt_absorption_line, **kwargs)
        self._set_paramhints_prefix()

    def guess(self, data, vrad, b, n, wavegrid=None, **kwargs):
        pars = guess_NA(self, data, vrad, b, n, wavegrid)

        return update_param_vals(pars, self.prefix, **kwargs)

