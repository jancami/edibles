from edibles.functions.continuum_guess import generate_continuum
from edibles.functions.find_f_known import find_F
from edibles.models import Cont1D, VoigtAbsorptionLine


def createLine(name, lam_0, b, d, tau_0):
    '''Creates an instance of a VoigtAbsorptionLine object

    INPUT:
        name:
        data:        [tuple]  In the form (wave, flux)
        line params: [floats]
            lam_0
            b
            d
            N
            f_known
    OUTPUT:
        line model [object instance]

    '''
    line = VoigtAbsorptionLine(name=name)
    line.lam_0          = lam_0
    line.b              = b
    line.d              = d
    line.tau_0          = tau_0

    return line


def createKnownLine(name, lam_0, b, d, N, f_known):
    '''Creates an instance of a VoigtAbsorptionLine object 
    where the oscillator strength is KNOWN

    INPUT:
        name:
        data:        [tuple]  In the form (wave, flux)

        line params: [floats]
            lam_0
            b
            d
            N
            f_known
    OUTPUT:
        line model [object instance]

    '''
    line = VoigtAbsorptionLine(name=name)

    line.N.hidden = False
    line.f.hidden = False
    line.tau_0.hidden = True

    line.N.frozen = False
    line.tau_0.frozen = True

    line.lam_0          = lam_0
    line.b              = b
    line.d              = d
    line.N              = N
    line.f              = f_known

    return line


def createKnownCloud(name, num_lines, lam_0, b, d, N, f_known):
    '''Creates multiple instances of a VoigtAbsorptionLine object 
    where the Oscillator strength is KNOWN

    INPUT:
        name:
        data:        [tuple]  In the form (wave, flux)
        num_lines:   [int]
        line params: [lists of floats]
            lam_0
            b
            d
            N
            f_known
    OUTPUT:
        line model [object instance]

    '''
    line0 = createKnownLine(name[0], lam_0[0], b[0], d[0], N[0], f_known[0])
    cloud = line0

    if num_lines > 1:
        for i in range(1, num_lines):
            line = createKnownLine(name[i], lam_0[i], b[i], d[i], N[i], f_known[i])
            line.N = line0.N
            line.b = line0.b

            cloud *= line
    return cloud


def createCont(data, n_points):
    '''Creates an instance of a Cont1D object

    INPUT:
        name:
        data:       [tuple]  In the form (wave, flux)
        n_points:

    OUTPUT:
        line model [object instance]

    '''
    n_piece = n_points - 1
    y_spline, y_points= generate_continuum(data, delta_v=1000, n_piece=n_piece)

    cont = Cont1D()
        # always at least 2 points / 1 piece
    if n_points >= 1:
        cont.y1            = y_points[0]
        cont.y1.frozen     = False
    if n_points >= 2:
        cont.y2            = y_points[1]
        cont.y2.frozen     = False
    if n_points >= 3:
        cont.y3            = y_points[2]
        cont.y3.frozen     = False
    if n_points >= 4:
        cont.y4            = y_points[3]
        cont.y4.frozen     = False
    if n_points >= 5:
        cont.y5            = y_points[4]
        cont.y5.frozen     = False
    if n_points >= 6:
        cont.y6            = y_points[5]
        cont.y6.frozen     = False
    if n_points >= 7:
        cont.y7            = y_points[6]
        cont.y7.frozen     = False
    if n_points >= 8:
        cont.y8            = y_points[7]
        cont.y8.frozen     = False

    return cont
