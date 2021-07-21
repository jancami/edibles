import numpy as np
import inspect
import collections
from scipy.interpolate import CubicSpline
from lmfit import Model
from lmfit.models import update_param_vals

from edibles.utils.voigt import voigtAbsorptionLine


def guess_voigt(model, data, x):
    """Estimate parameters of voigtAbsorptionLine, create params

    Args:
        model: VoigtModel
        data: flux data
        x: wavelength data
    Returns:
        lmfit.Parameters: Guessed parameters of voigt line

    TODO: This function could be extended in the future to guess peak
        parameters of multiple types of models.
    """

    if x is None:
        raise ValueError('x does not exist')
    if data is None:
        raise ValueError('y does not exist')

    x = np.asarray(x)
    y = np.asarray(data)

    lam_0 = x[np.argmin(y)]
    b = 2
    d = 0.001
    tau_0 = 0.1

    voigt_pars = model.make_params(lam_0=lam_0, b=b, d=d, tau_0=tau_0)

    return voigt_pars


class VoigtModel(Model):
    """A model of the astronomical Voigt function."""


    def __init__(self, independent_vars=['x'], prefix='', nan_policy='raise',
                 **kwargs):
        kwargs.update({'prefix': prefix, 'nan_policy': nan_policy,
                       'independent_vars': independent_vars})
        super().__init__(voigtAbsorptionLine, **kwargs)
        self._set_paramhints_prefix()


    def _set_paramhints_prefix(self):
        self.set_param_hint('lam_0', min=0, max=12000)
        self.set_param_hint('b', min=0, max=30)
        self.set_param_hint('d', min=0, max=10)
        self.set_param_hint('tau_0', min=0, max=5)


    def guess(self, data, x=None, **kwargs):
        """Estimate initial model parameter values from data.

        Args:
            data (array_like): y data points
            x (array_like): x data points

        Returns:
            lmfit.parameter.Parameters: Guessed parameters

        """
        pars = guess_voigt(self, data, x)

        return update_param_vals(pars, self.prefix, **kwargs)


class ContinuumModel(Model):
    """A model that puts a cubic spline through a small number (max 10) of evenly spaced
    anchor points, specified by ``n_anchors``. Only the y value of the anchor points is fit.

    Args:
        n_anchors (int): number of anchor points to fit spline through.
        independent_vars : ['x'] Arguments to func that are independent variables.
        prefix (str): optional, String to prepend to parameter names, needed to add two Models
            that have parameter names in common.
        nan_policy (str): optional, How to handle NaN and missing values in data.
            Must be one of: 'raise' (default), 'propagate', or 'omit'. See Notes below.
        **kwargs : optional,  Keyword arguments to pass to :class:`Model`.

    Notes:
        1. nan_policy sets what to do when a NaN or missing value is seen in the
        data. Should be one of:
        - 'raise' : Raise a ValueError (default)
        - 'propagate' : do nothing
        - 'omit' : drop missing data

    """


    def __init__(self, n_anchors, independent_vars=["x"], prefix="", nan_policy="raise", verbose=0, **kwargs):

        self.ANCHORS_ERR = "n_anchors (%.1f) must be an integer less than or equal to 10"

        kwargs.update({"prefix": prefix, "nan_policy": nan_policy,
                       "independent_vars": independent_vars})

        if not isinstance(n_anchors, int):
            raise TypeError(self.ANCHORS_ERR % n_anchors)

        self.ynames = ["y_%i" % (i) for i in range(n_anchors)]
        self.xnames = ["x_%i" % (i) for i in range(n_anchors)]
        kwargs["param_names"] = self.xnames + self.ynames

        params = {}
        for name in kwargs['param_names']:
            if name[0] == 'y':
                params[name] = 1
            if name[0] == 'x':
                params[name] = -999

        self.n_anchors = n_anchors
        self.verbose = verbose
        # if verbose >=2, print x and y anchors each time


        def cont(x, x_0=-999, y_0=1, **kwargs):

            spacing = np.linspace(np.min(x), np.max(x), self.n_anchors)
            x = np.asarray(x)
            spacing_idx = [np.argmin(np.abs(x - space)) for space in spacing]

            x_anchors = [x_0]
            y_anchors = [y_0]

            for arg in kwargs:
                if arg[0] == 'x':
                    x_anchors.append(kwargs[arg])
                elif arg[0] == 'y':
                    y_anchors.append(kwargs[arg])

            if all(anchor == -999 for anchor in x_anchors):
                x_anchors = [x[i] for i in spacing_idx]

            spline = CubicSpline(x_anchors, y_anchors)
            if self.verbose >= 2:
                print("====== Spline Continuum ======")
                x_anchors_p = ["%.2f" % item for item in x_anchors]
                print("Xs: ", x_anchors_p)

                y_anchors_p = ["%.2f" % item for item in y_anchors]
                print("Ys: ", y_anchors_p)

            return spline(x)


        sig = inspect.signature(cont)

        base_par = inspect.signature(cont).parameters['x_0']

        d = {'x': sig.parameters['x']}

        for i in range(n_anchors):

            y_key = 'y_' + str(i)
            x_key = 'x_' + str(i)

            y_val = base_par.replace(name=y_key)
            x_val = base_par.replace(name=x_key)

            d[y_key] = y_val
            d[x_key] = x_val

        d = collections.OrderedDict(d)

        cont.__signature__ = sig.replace(parameters=tuple(d.values()))

        super().__init__(cont, **kwargs)


    def guess(self, data, x=None, **kwargs):
        """
        Estimate initial anchor points through a dataset for the fitting of a cubic spline.

        Args:
            data (array_like): y data points
            x (array_like): x data points

        Returns:
            lmfit.parameter.Parameters: Guessed parameters
        """

        pars = self.make_params()

        spacing = np.linspace(np.min(x), np.max(x), self.n_anchors)
        spacing_idx = [(np.abs(x - space)).argmin() for space in spacing]

        data = np.asarray(data)

        x_anchors = [x[i] for i in spacing_idx]
        for i, coef in enumerate(x_anchors[::1]):
            pars['%sx_%i' % (self.prefix, i)].set(value=coef, vary=False)

        y_anchors = []
        for i in spacing_idx:
            if i == 0:
                ymin = i
                ymax = i + int(len(x) / self.n_anchors)
                y_anchors.append(np.median(data[ymin:ymax]))
            elif i == len(x) - 1:
                ymin = i - int(len(x) / self.n_anchors)
                ymax = i
                y_anchors.append(np.median(data[ymin:ymax]))
            else:
                ymin = i - int(len(x) / self.n_anchors)
                ymax = i + int(len(x) / self.n_anchors)
                y_anchors.append(np.median(data[ymin:ymax]))

        for i, coef in enumerate(y_anchors[::1]):
            pars['%sy_%i' % (self.prefix, i)].set(value=coef, min=0, max=2 * np.max(data))

        return update_param_vals(pars, self.prefix, **kwargs)


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    from edibles.utils.edibles_spectrum import EdiblesSpectrum

    filename = "/HD170740/RED_860/HD170740_w860_redl_20140915_O12.fits"
    sp = EdiblesSpectrum(filename)
    print(sp.target)
    sp.getSpectrum(xmin=7661, xmax=7670)
    # sp.flux = sp.flux / np.median(sp.flux)

    cont_model = ContinuumModel(n_anchors=4)
    cont_pars = cont_model.guess(sp.flux, x=sp.wave)

    result = cont_model.fit(data=sp.flux, params=cont_pars, x=sp.wave)
    result.plot_fit()
    plt.show()

    voigt_model = VoigtModel()
    voigt_pars = voigt_model.guess(sp.flux, x=sp.wave)

    model = cont_model * voigt_model
    pars = cont_pars + voigt_pars

    result = model.fit(data=sp.flux, params=pars, x=sp.wave)

    result.plot_fit()
    print(result.fit_report())
    plt.show()
