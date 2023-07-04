import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

from edibles.models import ContinuumModel, VoigtModel
from edibles.utils.edibles_spectrum import EdiblesSpectrum

def fit_NaI_Lines(target, date):
    """A function to fit Na I 2S-2P doublet lines - still very basic


    NaI1 = 3302.368 -> from edibles_linelist_atoms
    NaI2 = 3302.978 -> from edibles_linelist_atoms
    NaI_blue_diff = 0.61 -> lab data from https://www.pa.uky.edu/~peter/newpage/

    NaI3 = 5889.951 -> from edibles_linelist_atoms
    NaI4 = 5895.924 -> from edibles_linelist_atoms
    NaI_red_diff = 5.975 -> lab data from https://www.pa.uky.edu/~peter/newpage/

    blue_red_diff = 2587.583 -> NaI3 - NaI1

    """

    wavelength = 3300

    pythia = EdiblesOracle()
    bluelist = pythia.GetObsListByWavelength(wavelength, OrdersOnly=True)

    files = []
    for filename in bluelist:
        if target in filename:
            if date in filename:
                files.append(filename)

    print(files)

    sp = EdiblesSpectrum(files[0])
    print(sp.target)
    sp.getSpectrum(xmin=3300, xmax=3305)

    sigma = np.std(sp.flux)
    prominence = sigma
    peaks, _ = find_peaks(-sp.flux, prominence=prominence)
    peak_wavelengths = [sp.wave[i] for i in peaks]

    # #########################################################################

    cont = ContinuumModel(n_anchors=4)
    cont_pars = cont.guess(sp.flux, x=sp.wave)

    # #########################################################################

    voigt1 = VoigtModel(prefix='v1_')
    # voigt1_pars = voigt1.make_params(lam_0=3302, b=1, d=0.001, tau_0=0.01)
    voigt1_pars = voigt1.guess(sp.flux, x=sp.wave)

    # #########################################################################

    voigt2 = VoigtModel(prefix='v2_')
    voigt2_pars = voigt2.make_params(lam_0=peak_wavelengths[1], b=1, d=0.001, tau_0=0.01)

    voigt2_pars['v2_lam_0'].set(expr='v1_lam_0 + 0.61')
    voigt2_pars['v2_b'].set(expr='v1_b')

    # #########################################################################

    model = cont * voigt1 * voigt2
    pars = cont_pars + voigt1_pars + voigt2_pars

    result = model.fit(data=sp.flux, params=pars, x=sp.wave)

    # #########################################################################

    result.params.pretty_print()
    print('Ratio: ', result.params['v1_tau_0'] / result.params['v2_tau_0'])
    print(result.fit_report())

    f, [ax1, ax2] = plt.subplots(ncols=2)
    result.plot_fit(ax1)
    ax1.set_title('Na I 3303')

    print()
    print()

    # #########################################################################
    # #########################################################################

    wavelength = 5890

    pythia = EdiblesOracle()
    redlist = pythia.GetObsListByWavelength(wavelength, OrdersOnly=True)

    files = []
    for filename in redlist:
        if target in filename:
            if date in filename:
                files.append(filename)

    print(files)

    sp = EdiblesSpectrum(files[1])
    print(sp.target)
    sp.getSpectrum(xmin=5885, xmax=5900)

    # #########################################################################

    cont = ContinuumModel(n_anchors=4)
    cont_pars = cont.guess(sp.flux, x=sp.wave)

    # #########################################################################

    prominence = (np.max(sp.flux) - np.min(sp.flux)) * 0.5
    peaks, _ = find_peaks(-sp.flux, prominence=prominence)
    peak_wavelengths = [sp.wave[i] for i in peaks]

    voigt3 = VoigtModel(prefix='v3_')
    voigt3_pars = voigt3.make_params(lam_0=peak_wavelengths[0], b=1, d=0.001, tau_0=0.4)
    # voigt3_pars = voigt3.guess(sp.flux, x=sp.wave)

    # #########################################################################

    voigt4 = VoigtModel(prefix='v4_')
    voigt4_pars = voigt4.make_params(lam_0=5896, b=1, d=0.001, tau_0=0.1)

    voigt4_pars['v4_lam_0'].set(expr='v3_lam_0 + 5.975')
    voigt4_pars['v4_b'].set(expr='v3_b')

    # #########################################################################

    model = cont * voigt3 * voigt4
    pars = cont_pars + voigt3_pars + voigt4_pars

    result = model.fit(data=sp.flux, params=pars, x=sp.wave)

    # #########################################################################

    result.params.pretty_print()
    print('Ratio: ', result.params['v3_tau_0'] / result.params['v4_tau_0'])
    print(result.fit_report())

    result.plot_fit(ax=ax2)
    ax2.set_title('Na I 5890')
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":

    from edibles.utils.edibles_oracle import EdiblesOracle

    target = 'HD170740'
    date = '20140916'

    dates = ['20140916', '20150424', '20160505', '20160612', '20170701']

    for date in dates:
        fit_NaI_Lines(target='HD170740', date=date)
