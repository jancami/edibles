import numpy as np


def print_results(wave, fit_parm, fit_err, fit_norm):

    # compute scaled uncertainties
    DOF = len(wave) - len(fit_parm)
    PCERROR = fit_err * np.sqrt(fit_norm/DOF)

    # decompose the results
    continuum = fit_parm[:3]
    fit_parms = fit_parm[3:]
