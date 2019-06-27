from edibles.functions.continuum_guess import generate_continuum
from edibles.models import Cont1D, VoigtAbsorptionLine


def create_line(name, lam_0, b, d, tau_0):
    '''Creates an instance of a VoigtAbsorptionLine object

    INPUT:
        data:        [tuple]  In the form (wave, flux)
        line params: lam_0, b, d, tau_0 [floats]

    OUTPUT:
        line model [object instance]

    '''

    line = VoigtAbsorptionLine(name=name)
    # line.name           = name
    line.lam_0          = lam_0
    line.b              = b
    line.d              = d
    line.tau_0          = tau_0

    return line


def create_cont(data, n_points):
    '''Creates an instance of a Cont1D object

    INPUT:
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







