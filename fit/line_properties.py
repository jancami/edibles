import numpy as np

def line_properties(cw, lines):

    # read atomic dataset
    atomfile = '../data/atomic_data.dat'
    lineList = np.loadtxt(atomfile, dtype=[
                                    ('trans', 'S13'),
                                    ('ion',    'S6'),
                                    ('l0',     'f4'),
                                    ('f',      'f4'),
                                    ('gam',    'f4'),
                                    ('mass',   'f4')  ])

    l0_temp = []
    f_temp = []
    gam_temp = []
    for loop_l in range(len(lines)):
        index = lineList['trans'].tolist().index(lines[loop_l])
        tag, ion, l0, f, gam, mass = lineList[index]
        l0_temp.append(l0)
        f_temp.append(f)
        gam_temp.append(gam)

    # repeat the pattern for all transitions
    l0_list = []
    f_list = []
    gam_list = []
    for loop_l in range(len(cw)):
        for loop_l2 in range(len(lines)):
            l0_list.append(l0_temp[loop_l2])
            f_list.append(f_temp[loop_l2])
            gam_list.append(gam_temp[loop_l2])



    return l0_list, f_list, gam_list
