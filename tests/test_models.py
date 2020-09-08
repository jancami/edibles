import pytest

from edibles.utils.edibles_spectrum import EdiblesSpectrum
from edibles.models import ContinuumModel, VoigtModel


def testModels(filename="tests/HD170740_w860_redl_20140915_O12.fits"):


    method = 'least_squares'

    sp = EdiblesSpectrum(filename, noDATADIR=True)
    sp.getSpectrum(xmin=7661, xmax=7670)

    n_anchors = 4

    cont_model = ContinuumModel(n_anchors=n_anchors)
    assert len(cont_model.param_names) == n_anchors * 2
    assert cont_model.n_anchors == n_anchors

    with pytest.raises(TypeError):
        assert ContinuumModel()

    cont_pars = cont_model.guess(sp.flux, x=sp.wave)
    assert len(cont_pars) == len(cont_model.param_names)

    for name in cont_model.param_names:
        assert cont_pars[name].value is not None

    voigt = VoigtModel(prefix='voigt_')
    assert voigt.prefix == 'voigt_'
    assert len(voigt.param_names) == 4

    voigt_pars = voigt.guess(sp.flux, x=sp.wave)
    assert len(voigt_pars) == len(voigt.param_names)

    result = voigt.fit(data=sp.flux, params=voigt_pars, x=sp.wave, method=method)
    assert len(result.params) == n_anchors

    out = voigt.eval(data=sp.flux, params=result.params, x=sp.wave)
    assert len(out) == len(sp.flux)



if __name__ == "__main__":

    filename = "HD170740_w860_redl_20140915_O12.fits"
    testModels(filename=filename)
