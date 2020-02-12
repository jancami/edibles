from __future__ import print_function

from edibles.edibles.fit.models.create_model import *
from edibles.edibles.functions.edibles_spectrum import EdiblesSpectrum
from edibles.edibles.fit.models.models import Sightline
from edibles.edibles.fit.fit import fit
import matplotlib.pyplot as plt


def testAdvancedFit():
    star_name = "HD170740"
    file = "/HD170740/RED_860/HD170740_w860_redl_20140916_O12.fits"
    xmin = 7661.0
    xmax = 7670.0
    sp = EdiblesSpectrum(file)
    data = sp.getSpectrum(xmin, xmax)
    wave, flux = data

    # Cont parameters
    n_points = 4
    cont = createCont(data, n_points)

    slightline = Sightline(star_name=star_name, cont=cont)

    slightline.addSource(source_name="Telluric", b=1.07, d=0.046)
    slightline.addLine(name="tell_1", lam_0=7664.8, tau_0=0.75)
    slightline.addLine(name="tell_2", lam_0=7666, tau_0=0.75)

    slightline.addSource(source_name="Interstellar", b=1.0, d=0.001)
    slightline.addLine(name="KI_1", lam_0=7665.3, tau_0=0.1)
    slightline.addLine(name="KI_2", lam_0=7665.35, tau_0=0.05)

    fit_model = fit(star_name, data, slightline.model)

    file = "/HD148937/BLUE_346/HD148937_w346_blue_20150817_O11.fits"
    star = "HD148937"
    wave, flux = EdiblesSpectrum(file).getSpectrum(3301.5, 3304)
    plt.plot(wave, flux)
    plt.title(file)
    plt.show()

    cont = createCont((wave, flux), n_points=3)
    slightline = Sightline(star_name=star, cont=cont)
    slightline.addSource("Source 1", 1.42631e-07, 0.036356)
    slightline.addLine("NaI_1", lam_0=3302.46, tau_0=0.06)
    slightline.addLine("NaI_2", lam_0=3303.1, tau_0=0.03)

    slightline.dupSource("Source 1", "Source 2", 1.00005)
    slightline.dupSource("Source 1", "Source 3", 0.99995)

    fit_model = fit(star, (wave, flux), slightline.model, breakdown=True)

    return


if __name__ == "__main__":
    testAdvancedFit()
