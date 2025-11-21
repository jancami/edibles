#############
# Test the Bayes, and maybe F-test method on stopping adding v-components
# Target Na 3300 for HD183143. We know it has at least two major components
# So let's say what does the Bayesian info attribute has to say
#############

# basic
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import copy
import math
from scipy.stats import f as scipy_f

# lmfit related
from lmfit import Model
from lmfit.models import update_param_vals

# model building related
import inspect
import collections

# edibles repo -- models and data IO
from edibles.models import ContinuumModel
from edibles.utils.ISLineFitter import ISLineModel
from edibles.utils.edibles_oracle import EdiblesOracle
from edibles.utils.edibles_spectrum import EdiblesSpectrum


####################
# function posterior_calc
####################
# Notes
# 1. This function is based on Klay's master project,
#    except getting ln of posterior.
# 2. !! The ranges of parameters will affect "range_term" as well as the
#    determinant of covariance matrix. But they would not cancel out,
#    although the difference is small in ln unit
# 3. !! The calculated posterior is highly determined by chi-sqr, by a few
#    magnitudes in ln unit. So the punishment of adding new parameters
#    cannot overtake the benefit of reduced chi-sqr. LMFIT assume noise=1
#    and that could be the problem.
def posterior_calc(result):
    range_term = 0

    com_idx = 0
    while True:
        if "V_off_Cloud%i" % com_idx in result.params.keys():
            for var in ["V_off", "b", "N"]:
                range = result.params[var + "_Cloud%i" % com_idx].max - result.params[var + "_Cloud%i" % com_idx].min
                range_term = range_term + np.log(range)

            com_idx = com_idx + 1
        else:
            break

    anchor_idx = 0
    while True:
        if "y_%i" % anchor_idx in result.params.keys():
            range = result.params["y_%i" % anchor_idx].max - result.params["y_%i" % anchor_idx].min
            range_term = range_term + np.log(range)

            anchor_idx = anchor_idx + 1
        else:
            break

    try:
        det = 1 / np.linalg.det(result.covar)
        det = np.log(det)
    except np.linalg.LinAlgError:
        print("Something wrong with Covar matrix...")
        if result.covar is None:
            print("It's None?")
    # Do we really need this part of debugging?

    a = math.factorial(com_idx) * \
        (2 * np.pi) ** ((3 * com_idx + anchor_idx) / 2)
    a = np.log(a)
    # b = range_term * (np.sqrt(det))
    b = range_term + 2 * det
    # print("range_term", range_term)
    # print("det", det)
    # print("b", b)
    c = -result.chisqr / 2
    posterior = a - b + c
    # posterior = a / b * c
    # print("a",a)
    # print("b",b)
    # print("c",c)

    return posterior


####################
# function for F test
####################
# Notes
# 1. This function is directly taken from Klay's code, only changed
#    variable names.
# 2. !! No. of parameters in the spline continuum needs more work,
#    e.g. x_anchor are always fixed, and in case of "no-continuum
#    model" all parameters in continuum would be fixed.
# 3. !! The output p_value seems tricky. In the following example,
#    adding the first three components yields p > 0.999, not
#    something close to 0.
# 4. !! Also, Klay require p > alpha = 0.05, and is strange. A
#    reasonable amend would be 1 - P < alpha = 0.05 I suppose.
def f_test(x_Grid, result_old, result_new, alpha=0.05):
    # print(len(result_new.params), len(result_old.params))
    df1 = len(result_new.params) - len(result_old.params)
    df2 = len(x_Grid) - len(result_new.params) + 4

    num = (result_old.chisqr - result_new.chisqr) / df1
    denom = result_new.chisqr / df2
    f = num / denom

    p_value = scipy_f.cdf(f, df1, df2)
    return p_value


####################
# Update to make for ISLineFitter
####################
# V  1. Plot individual components of line model. Say, a three panel
#    plot for original data, normalized data with components, and
#    residuals?
# 2. Better estimate on V_offs when the major components have been
#    addressed. Algorithm needed.
# 3. Output result in the text form, and creating/appending data
#    file.
# 4. Benchmark for measurements, check atomic data, and better
#    estimate
# 5. Constrain the range on N_tot, e.g. use log unit or provides
#    known magnitude (I prefer the latter).
# 6. Estimate SNR of the data, especially if we want posterior.
# 7. For the future, use Bayesian method to add stellar/telluric
#    line in the model for fitting.

####################
# Load data for 3300 segment
####################
pythia = EdiblesOracle()
List = pythia.getFilteredObsList(object=["HD 183143"], OrdersOnly=True, Wave=3302.0)
filename = List.values.tolist()[1]

wrange = [3301.0, 3304.0]
print(filename)
sp = EdiblesSpectrum(filename)
wave, flux = sp.bary_wave, sp.flux
idx = np.where((wave > wrange[0]) & (wave < wrange[1]))
wave2fit, flux2fit = wave[idx], flux[idx]

####################
# Continuum Model
####################
continuum_model = ContinuumModel(n_anchors=4)
cont_pars = continuum_model.guess(flux2fit, wave2fit)

####################
# Line Model
####################
model_all, result_all = [], []
v_off_known = [-11.1, 4.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
aics, bics, klays, chisqrs, fs = [], [], [], [], []
fig = plt.figure(figsize=[12, 6], dpi=180)
for i in range(6):
    line_model = ISLineModel(i)
    model2fit = continuum_model * line_model
    model_all.append(model2fit)

    pars_guess = copy.deepcopy(cont_pars)
    print(pars_guess)
    if i == 0:
        pars_guess.update(line_model.guess(V_off=[]))
    else:
        pars_guess.update(line_model.guess(V_off=v_off_known[0:i]))

    print(i)
    result = model2fit.fit(data=flux2fit, params=pars_guess, x=wave2fit)
    result_all.append(result)

    comps = result.eval_components(x=wave2fit)
    continuum = comps["cont"]

    ax = fig.add_subplot(2, 3, i + 1)
    ax.plot(wave2fit, flux2fit / continuum)
    if i >= 1:
        flux_comps = line_model.calcIndividualComponent(result.params, wave2fit)
        for com_idx, flux_single in enumerate(flux_comps):
            ax.plot(wave2fit, flux_single, label="Comp %i" % com_idx)
        ax.legend()
    ax.set_xlabel("n_components = %i" % i)

    aics.append(result.aic)
    bics.append(result.bic)
    chisqrs.append(result.chisqr)
    klays.append(posterior_calc(result))
    if i >= 1:
        fs.append(f_test(wave2fit, result_all[-2], result_all[-1]))

    print("=" * 40)
    print(result.params["x_0"].vary)

print("AIC: ", aics)
print("BIC: ", bics)
print("ChiSqr: ", chisqrs)
print("Posterior", klays)
print("F test", fs)

plt.tight_layout()
plt.show()
