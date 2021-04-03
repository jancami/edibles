import numpy as np
from scipy.stats import f as scipy_f
import statistics

def compareModels(old, new, iter, values):
    total_f = 0

    df1 = 4
    df2 = values - (len(new.params))

    num = (old.chisqr - new.chisqr) / df1
    denom = new.chisqr / df2

    f = num / denom
    print("f: ", f)

    alpha = 0.05  # Or whatever you want your alpha to be.
    p_value = scipy_f.cdf(f, df1, df2)
    print("p: ", p_value)
    if p_value > alpha:
        # Reject the null hypothesis that Var(X) == Var(Y)
        total_f += 1

    num_to_stop = np.ceil(iter / 2.0)

    if total_f >= 1:
        print(
            str(iter)
            + ' lines are better than '
            + str(iter - 1)
        )
        return False
    else:
        return True

def CompareModelBeyes(old, new, iter, values):
    total_bayes = 0

    """ posterior_calc(sightlines)

    if iter > 0:
        if np.isinf(sightline.posteriors[-1]):
            pass
        elif sightline.posteriors[-1] == 0:
            total_bayes += 1
        elif sightline.posteriors[-1] / old.posteriors > 1:
            total_bayes += 1
        else:
            pass

    else:
        total_bayes += 1

    num_to_stop = np.ceil(iter / 2.0)

    if total_bayes >= 1:
        print(
            str(iter)
            + ' lines are better than '
            + str(iter - 1)
        )
        return False
    else:
        return True"""
