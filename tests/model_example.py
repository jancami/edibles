from edibles.utils.edibles_spectrum import EdiblesSpectrum
from edibles.models import ContinuumModel, VoigtModel


def modelExample(filename="tests/HD170740_w860_redl_20140915_O12.fits"):


    method = 'least_squares'

    sp = EdiblesSpectrum(filename, noDATADIR=True)
    print(sp.target)
    subset = sp.getSpectrum(xmin=7661, xmax=7670)

    def new_line(y, x, line_num):

        voigt = VoigtModel(prefix='voigt' + str(line_num) + '_')
        voigt_pars = voigt.guess(y, x=x)

        return voigt, voigt_pars


    cont_model = ContinuumModel(n_anchors=4)
    cont_pars = cont_model.guess(subset.flux, x=subset.wave)
    model = cont_model
    pars = cont_pars

    result = model.fit(data=subset.flux, params=pars, x=subset.wave, method=method)
    out = cont_model.eval(data=subset.flux, params=result.params, x=subset.wave)
    resid = subset.flux - out

    # result.plot_fit()
    # plt.show()

    # #################################################################################

    voigt1, voigt1_pars = new_line(resid, x=subset.wave, line_num=1)
    model = model * voigt1
    pars = result.params + voigt1_pars

    result = model.fit(data=subset.flux, params=pars, x=subset.wave, method=method)
    out = model.eval(data=subset.flux, params=result.params, x=subset.wave)
    resid = subset.flux - out

    # result.plot_fit()
    # plt.show()

    # #################################################################################

    voigt2, voigt2_pars = new_line(resid, x=subset.wave, line_num=2)
    model = model * voigt2
    pars = result.params + voigt2_pars

    result = model.fit(data=subset.flux, params=pars, x=subset.wave, method=method)
    out = model.eval(data=subset.flux, params=result.params, x=subset.wave)
    resid = subset.flux - out

    # result.plot_fit()
    # plt.show()

    # #################################################################################

    voigt3, voigt3_pars = new_line(resid, x=subset.wave, line_num=3)
    model = model * voigt3
    pars = result.params + voigt3_pars

    result = model.fit(data=subset.flux, params=pars, x=subset.wave, method=method)
    out = model.eval(data=subset.flux, params=result.params, x=subset.wave)
    resid = subset.flux - out

    # result.plot_fit()
    # plt.show()

    # #################################################################################

    voigt4, voigt4_pars = new_line(resid, x=subset.wave, line_num=4)
    model = model * voigt4
    pars = result.params + voigt4_pars

    result = model.fit(data=subset.flux, params=pars, x=subset.wave, method=method)
    out = model.eval(data=subset.flux, params=result.params, x=subset.wave)
    resid = subset.flux - out


    # #################################################################################

    # print(result.fit_report())

    # plt.plot(subset.wave, subset.flux)
    # plt.plot(subset.wave, out)
    # plt.plot(subset.wave, resid)
    # plt.show()



if __name__ == "__main__":

    filename = "HD170740_w860_redl_20140915_O12.fits"
    modelExample(filename=filename)
