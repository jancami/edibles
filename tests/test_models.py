import pytest

from edibles.utils.edibles_spectrum import EdiblesSpectrum
from edibles.models import ContinuumModel, VoigtModel


def testModels(filename="tests/HD170740_w860_redl_20140915_O12.fits"):


    method = 'least_squares'

    sp = EdiblesSpectrum(filename, noDATADIR=True)
    subset = sp.getSpectrum(xmin=7661, xmax=7670)

    n_anchors = 4

    cont_model = ContinuumModel(n_anchors=n_anchors)
    assert len(cont_model.param_names) == n_anchors
    assert cont_model.n_anchors == n_anchors

    with pytest.raises(TypeError):
        assert ContinuumModel()
    with pytest.raises(TypeError):
        assert ContinuumModel(n_anchors=11)

    cont_pars = cont_model.guess(subset.flux, x=subset.wave)
    assert len(cont_pars) == len(cont_model.param_names)

    for name in cont_model.param_names:
        assert cont_pars[name].value is not None



    voigt = VoigtModel(prefix='voigt_')
    assert voigt.prefix == 'voigt_'
    assert len(voigt.param_names) == 4

    voigt_pars = voigt.guess(subset.flux, x=subset.wave)
    assert len(voigt_pars) == len(voigt.param_names)

    result = voigt.fit(data=subset.flux, params=voigt_pars, x=subset.wave, method=method)
    assert len(result.params) == n_anchors

    out = voigt.eval(data=subset.flux, params=result.params, x=subset.wave)
    assert len(out) == len(subset.flux)



if __name__ == "__main__":

    filename = "HD170740_w860_redl_20140915_O12.fits"
    testModels(filename=filename)
