import numpy as np

def CompareModelBeyes(old, new, iter1,iter2):
    total_bayes = 0

    oldPost = posterior_calc(old, iter1)
    newPost = posterior_calc(new, iter2)

    if np.isinf(newPost):
        pass
    elif newPost == 0:
        total_bayes += 1
    elif newPost / oldPost > 1:
        total_bayes += 1

    print ("New Post:", newPost, "Old Post:", oldPost, "Total Beyes:", total_bayes)
    if total_bayes >= 1:
        print(
            str(iter2)
            + ' lines are better than '
            + str(iter1)
        )
        return False

    else:
        print("Stop with", str(iter1))
        return True

def posterior_calc(sightline, lines):

    range_terms = []
    for i in range(lines):
        for name in sightline.params:
            if 'V_off_Cloud' + str(i) in name:
                v_off_name =name
        b_name = v_off_name.replace('V_off', 'b')
        N_name = v_off_name.replace('V_off', 'N')

        b_range = sightline.result.params[b_name].max - sightline.result.params[b_name].min
        v_off_range = sightline.result.params[v_off_name].max - sightline.result.params[v_off_name].min
        N_range = (
                sightline.result.params[N_name].max - sightline.result.params[N_name].min
        )
        range_term = b_range * v_off_range * N_range
        range_terms.append(range_term)

    denominator_term = 1


    for term in range_terms:
        denominator_term = denominator_term * term

    try:
        det = 1 / np.linalg.det(sightline.result.covar)
    except Exception as e:
        print(e)
        det = 0

    a = np.math.factorial(lines)*(2 * np.pi)**((4 * lines) / 2)
    b = denominator_term * (np.sqrt(det))
    c = np.exp((-sightline.chisqr *1000)/ 2)

    print('a: ', a)
    print('b: ', b)
    print('c: ', c)

    posterior = a  / b * c
    return posterior
