import matplotlib.pyplot as plt
import numpy as np
from edibles.functions.voigtMathematical import voigt_math


def decompose(x, fit_param, fit_error, fit_norm, num_params=4, print_data=False, plot=False):
    '''
    INPUT:

    fit_param  [ndarray]: Of the form - [amp, cent, alpha, gamma, amp, cent, alpha, gamma, ...]
    fit_error  [ndarray]: Of the form - [amp, cent, alpha, gamma, amp, cent, alpha, gamma, ...]
    fit_norm   [ndarray]: Of the form - [amp, cent, alpha, gamma, amp, cent, alpha, gamma, ...]
    num_params [int]:     Number of parameters for each line, default=4
    print_data [bool]:    Optional - print the data into a table
    plot       [bool]:    Optional - plot the data (NOT FINISHED)

    OUTPUT:

    params     [ndarray]: Of the form - [[p_1_1, p_1_2, p_1_3, ...], [p_2_1, p_2_2, p_2_3, ...], [...]]
    err_params [ndarray]: Of the form - [[p_1_1, p_1_2, p_1_3, ...], [p_2_1, p_2_2, p_2_3, ...], [...]]
    '''


    num_lines = int(np.ceil(len(fit_param) / float(num_params)))


    DOF = len(x) - len(fit_param)
    PCERROR = fit_error * np.sqrt(fit_norm/DOF)


    params = []
    err_params = []
    for i in range(num_lines):
        params.append([])
        err_params.append([])

    index=0
    for i in range(num_lines):
        for j in range(num_params):
            params[i].append(fit_param[index])
            err_params[i].append(PCERROR[index])
            index += 1

            if index == len(fit_param):
                break


    # Optional print functionality
    if print_data is True:

        print('=========================================================')
        print('                    *** results ***')
        for i in range(num_lines):


            print('')
            for j in range(len(params[i])):
                print('  param_{:.0f}_{:.0f}         :          {:.5f} +- {:.7f}'.format(i, j, params[i][j], err_params[i][j]))

            print('')

        print('=========================================================')

    # LOOKS LIKE:
    # =========================================================
    #                     *** results ***

    #   param_0_0         :          0.99985 +- 0.0007536
    #   param_0_1         :          0.00031 +- 0.0003013
    #   param_0_2         :          0.00010 +- 0.0001002
    #   param_0_3         :          0.00000 +- 0.0000000


    #   param_1_0         :          5890.00071 +- 0.0003346
    #   param_1_1         :          60640000.00000 +- 0.0000000
    #   param_1_2         :          3.33599 +- 0.0627548
    #   param_1_3         :          12.92497 +- 0.0399954

    # =========================================================



    # Optional plot functionality - NOT COMPLETE
    if plot is True:

        for i in range(num_lines):
            amp, cent, alpha, gamma = params[i]
            y = voigt_math(x, params[i][1], params[i][2], params[i][3])

            plt.figure()
            plt.plot(x, y)

        plt.show()
    return params, err_params