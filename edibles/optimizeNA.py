import numpy as np
import astropy.constants as cst
from lmfit import Model
from lmfit.models import update_param_vals
from edibles.utils.voigt_profile import voigt_absorption_line
from scipy.interpolate import interp1d
from scipy.stats import pearsonr

def calc_vrad(wave,flux):

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
    v = -50
    while v <= 50:
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

def guess_NA(model, data, b, n, x):
    if data is None:
        raise ValueError('y does not exist')

    y = np.asarray(data)
    v_rad = calc_vrad(x,y)
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

    def _set_paramhints_prefix(self):
        self.set_param_hint('v_rad', min=-35, max=35)
        self.set_param_hint('b', min=0, max=10)
        self.set_param_hint('N', min=0, max=18.9e15)

    def guess(self, data,  b, n, wavegrid=None, **kwargs):
        pars = guess_NA(self, data,  b, n, wavegrid)

        return update_param_vals(pars, self.prefix, **kwargs)

