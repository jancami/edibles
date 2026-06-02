"""
Modelling a rotational transition for a linear molecule.
"""
import matplotlib.pyplot as plt
import numpy as np
from edibles.projects.edr5_integration import transformations
import pandas as pd
from scipy.stats import norm

kb = 0.6950356  # cm-1/K
my_k = 1  # electronic angular momentum
my_dk = 1  # change of electronic angular momentum


# Functions to calculate wave number of transitions
def get_coeffs(jp, jpp):
    cbp = jp * (jp + 1)
    cdp = -jp ** 2 * (jp + 1) ** 2
    cbpp = -jpp * (jpp + 1)
    cdpp = jpp ** 2 * (jpp + 1) ** 2

    return np.array([cbp, cdp, cbpp, cdpp])


def nu_rot(bp, dp, bpp, dpp, jp, jpp):
    coeffs = get_coeffs(jp, jpp)
    return np.dot(np.array([bp, dp, bpp, dpp]), coeffs)


def honl_london_factor(j, jp, k, dk):
    # calculate Honl-London factors
    if jp == j + 1:  # R branch
        if dk == 1:
            hlf = (j + 2 + k) * (j + 1 + k) / ((j + 1) * (2 * j + 1))
        elif dk == -1:
            hlf = (j + 2 - k) * (j + 1 - k) / ((j + 1) * (2 * j + 1))
        elif dk == 0:
            hlf = (j + 1 + k) * (j + 1 - k) / ((j + 1) * (2 * j + 1))

    if jp == j:  # Q branch
        if dk == 1:
            hlf = (j + 1 + k) * (j - k) / (j * (j + 1))
        elif dk == -1:
            hlf = (j + 1 - k) * (j + k) / (j * (j + 1))
        elif dk == 0:
            hlf = (2 * j + 1) * k ** 2 / (j * (j + 1))

    if jp == j - 1:  # P branch
        if dk == 1:
            hlf = (j - 1 - k) * (j - k) / (j * (2 * j + 1))
        elif dk == -1:
            hlf = (j - 1 + k) * (j + k) / (j * (2 * j + 1))
        elif dk == 0:
            hlf = (j + k) * (j - k) / (j * (2 * j + 1))

    return hlf


# Calculate transition intensities from J'' to J'
def line_intensity(bpp, dpp, jp, jpp, t, k=my_k, dk=my_dk):
    gpp = 2 * jpp + 1  # Degeneracy of J'' state
    epp = bpp * jpp * (jpp + 1) - dpp * jpp ** 2 * (jpp + 1) ** 2  # rotational energy of J'' state
    # thermal level population at J''
    njpp = gpp * np.exp(-epp / kb / t)

    hlf = honl_london_factor(jpp, jp, k, dk)

    return njpp * hlf


def linear_stick_model(bp, dp, b, d, origin, j_min, j_max, t, k=my_k, dk=my_dk):
    my_line_list = pd.DataFrame()
    # R-branch
    j_start = np.max([1, j_min])
    for my_j in range(j_start, j_max):
        jpp = my_j - 1
        jp = my_j
        my_nu_rot = nu_rot(bp, dp, b, d, jp, jpp) + origin
        my_intensity = line_intensity(b, d, jp, jpp, t, k=k, dk=dk)
        my_line_list = pd.concat((my_line_list,
                                  pd.Series([jpp, my_nu_rot, my_intensity], index=["J'", 'nu_rot', 'intensity'],
                                            name=f'R({jpp})')), axis=1)

    # Q-branch
    if dk != 0:
        for my_j in range(j_min, j_max):
            jp = jpp = my_j
            my_nu_rot = nu_rot(bp, dp, b, d, jp, jpp) + origin
            my_intensity = line_intensity(b, d, jp, jpp, t, k=k, dk=dk)
            my_line_list = pd.concat((my_line_list,
                                      pd.Series([jpp, my_nu_rot, my_intensity], index=["J'", 'nu_rot', 'intensity'],
                                                name=f'Q({jpp})')), axis=1)

    # P-branch
    for my_j in range(j_min, j_max):
        jpp = my_j + 1
        jp = my_j
        my_nu_rot = nu_rot(bp, dp, b, d, jp, jpp) + origin
        my_intensity = line_intensity(b, d, jp, jpp, t, k=k, dk=dk)
        my_line_list = pd.concat((my_line_list,
                                  pd.Series([jpp, my_nu_rot, my_intensity], index=["J'", 'nu_rot', 'intensity'],
                                            name=f'P({jpp})')), axis=1)

    return my_line_list.T


def convolve_stick_model(my_line_list, profile=None, rv_width=None, special_gauss=None, res=1000, x_data=None):
    if special_gauss is None:
        special_gauss = {}

    if rv_width is not None:
        if x_data is None:
            mean_gauss = my_line_list['nu_rot'].mean() - transformations.rv_to_wavenumber(rv_width, my_line_list[
                'nu_rot'].mean())
            x_data = np.linspace(my_line_list['nu_rot'].min() - mean_gauss * 5,
                                 my_line_list['nu_rot'].max() + mean_gauss * 5, res)
        y_data = np.zeros(len(x_data))
        for name, row in my_line_list.iterrows():
            if name in special_gauss:
                wn_width = row.nu_rot - transformations.rv_to_wavenumber(special_gauss[name], row.nu_rot)
                profile = norm.pdf(x_data, row.nu_rot, wn_width) * row.intensity
                y_data += profile
            else:
                wn_width = row.nu_rot - transformations.rv_to_wavenumber(rv_width, row.nu_rot)
                profile = norm.pdf(x_data, row.nu_rot, wn_width) * row.intensity
                y_data += profile

    return np.array([x_data, y_data])


def main():
    fit_res = [1.826543168326506, -0.0005865002642848938, 2.163380314296096, 0.027666048496364372, 16128.387554692832]
    rot_coeffs = fit_res[:4]
    origin = fit_res[4]
    # line list for CN
    line_list = linear_stick_model(1.973, 0, 1.8997, 0, 25752.0, 0, 6, 10, k=0, dk=0)
    print(line_list)
    # print(line_list.loc['R(2)', 'intensity'] / line_list.loc['R(1)', 'intensity'])
    # convolved_model = convolve_stick_model(line_list, rv_width=1, special_gauss={'Q(2)': 0.6, 'Q(3)': 0.6})
    convolved_model = convolve_stick_model(line_list, rv_width=1)
    plt.figure(figsize=(15, 10))
    plt.plot(convolved_model[0] + origin, convolved_model[1], label='convolved model')
    plt.bar(line_list.loc[:, 'nu_rot'] + origin, line_list.loc[:, 'intensity'], width=0.1)
    for name, row in line_list.iterrows():
        plt.annotate(name, (row.nu_rot + origin, row.intensity))
    plt.xlabel(r'$\tilde\nu(J)$')
    plt.ylabel('Intensity')
    plt.show()


if __name__ == '__main__':
    main()
