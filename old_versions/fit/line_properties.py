import numpy as np
import edibles.functions.TextFileParser as tfp

def line_properties(cw, lines):

    # read atomic dataset
    atomfile = '../data/atomic_data.dat'

    trans = []
    ion = []
    l0 = []
    f = []
    gam = []
    mass = []

    lineList = tfp.parseTextFile(atomfile, delimiter='\t', header=3)
    for file_line in range(len(lineList)):
        try:

            if lineList[file_line][0][0] == '#':
                continue

            lineList[file_line][0] = lineList[file_line][0].split()
            trans.append(lineList[file_line][0][0])
            ion.append(lineList[file_line][0][1])
            l0.append(float(lineList[file_line][0][2]))
            f.append(float(lineList[file_line][0][3]))
            gam.append(float(lineList[file_line][0][4]))
            mass.append(float(lineList[file_line][0][5]))

        except IndexError:
            continue

    # turn lists into arrays & merge to one large array
    np.asarray(trans)
    np.asarray(ion)
    np.asarray(l0)
    np.asarray(f)
    np.asarray(gam)
    np.asarray(mass)
    lineList = np.column_stack((trans, ion, l0, f, gam, mass))



    l0_temp = []
    f_temp = []
    gam_temp = []

    for loop_l in range(len(lines)):
        index = lineList[:, 0].tolist().index(lines[loop_l])
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
            l0_list.append(float(l0_temp[loop_l2]))
            f_list.append(float(f_temp[loop_l2]))
            gam_list.append(float(gam_temp[loop_l2]))


    return l0_list, f_list, gam_list
