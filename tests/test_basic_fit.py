def testBasicFit():

    from edibles.edibles.utils.edibles_spectrum import EdiblesSpectrum
    from edibles.edibles.models.create_model import createCont
    from edibles.edibles.models.model import Sightline
    from edibles.edibles.fitter import fit

    file1 = "tests/testdata/HD170740_w860_redl_20140915_O12.fits"
    sp1 = EdiblesSpectrum(file1, noDATADIR=True)
    xmin = 7661.5
    xmax = 7669.0

    subset = sp1.getSpectrum(xmin=xmin, xmax=xmax)
    data = (subset["wave"], subset["flux"])

    cont = createCont(data, n_points=4)
    sightline = Sightline(star_name=sp1.target, cont=cont)
    sightline.addSource(source_name="Telluric", b=0.001, d=0.05)
    sightline.addLine(name="tell1", lam_0=7664.5, tau_0=0.6)
    sightline.addLine(name="tell2", lam_0=7666, tau_0=0.6)

    sightline.addSource(source_name="Interstellar", b=0.01, d=0.02)
    sightline.addLine(name="int1", lam_0=7665.2, tau_0=0.2)
    sightline.addLine(name="int2", lam_0=7665.3, tau_0=0.01)

    fit_model = fit(sp1.target, data, sightline.model, breakdown=False, silent=True)

    return fit_model


if __name__ == "__main__":
    fit_model = testBasicFit()
    print(fit_model)
