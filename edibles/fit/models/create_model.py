from edibles.edibles.functions.continuum_guess import generate_continuum
from edibles.edibles.fit.models.models import (
    Cont1D,
    VoigtAbsorptionLine,
    KnownVelocityLine,
)


def createCont(data, n_points=4):
    """
    Creates an spline model through a spectrum.


    :param data: Spectrum data in the form (wave, flux)
    :type data: tuple
    :param n_points: Number of anchor points to use - must be between 0 and 8!
    :type n_points: int

    :return: spline model through anchor points
    :rtype: object
    """

    assert (0 <= n_points <= 8), 'n_points must be between 0 and 8!'

    n_piece = n_points - 1
    y_spline, y_points = generate_continuum(data, delta_v=1000, n_piece=n_piece)

    cont = Cont1D()
    # always at least 2 points / 1 piece
    if n_points >= 1:
        cont.y1 = y_points[0]
        cont.y1.frozen = False
    if n_points >= 2:
        cont.y2 = y_points[1]
        cont.y2.frozen = False
    if n_points >= 3:
        cont.y3 = y_points[2]
        cont.y3.frozen = False
    if n_points >= 4:
        cont.y4 = y_points[3]
        cont.y4.frozen = False
    if n_points >= 5:
        cont.y5 = y_points[4]
        cont.y5.frozen = False
    if n_points >= 6:
        cont.y6 = y_points[5]
        cont.y6.frozen = False
    if n_points >= 7:
        cont.y7 = y_points[6]
        cont.y7.frozen = False
    if n_points >= 8:
        cont.y8 = y_points[7]
        cont.y8.frozen = False

    return cont


###############################################################
# OLDER CODE BELOW - Newer method available (see Sightline oblect)

# def createLine(name, lam_0, b, d, tau_0):
#     """Creates an instance of a VoigtAbsorptionLine object

#     INPUT:
#         name:
#         data:        [tuple]  In the form (wave, flux)
#         line params: [floats]
#             lam_0
#             b
#             d
#             tau_0
#     OUTPUT:
#         line model [object instance]

#     """
#     line = VoigtAbsorptionLine(name=name)
#     line.lam_0 = lam_0
#     line.b = b
#     line.d = d
#     line.tau_0 = tau_0

#     return line


# def createKnownLine(name, lam_0, b, d, N, f_known):
#     """Creates an instance of a VoigtAbsorptionLine object 
#     where the oscillator strength is KNOWN

#     INPUT:
#         name:
#         data:        [tuple]  In the form (wave, flux)

#         line params: [floats]
#             lam_0
#             b
#             d
#             N
#             f_known
#     OUTPUT:
#         line model [object instance]

#     """
#     line = VoigtAbsorptionLine(name=name)

#     line.N.hidden = False
#     line.f.hidden = False
#     line.tau_0.hidden = True

#     line.N.frozen = False
#     line.tau_0.frozen = True

#     line.lam_0 = lam_0
#     line.b = b
#     line.d = d
#     line.N = N
#     line.f = f_known

#     return line


# def createKnownVelocityLine(name, v_cloud, b, d, N, f_known, lab_lam_0):
#     """Creates an instance of a KnownVelocityLine object 
#     where the oscillator strength is KNOWN and the x axis is in velocity space

#     INPUT:
#         name:
#         data:        [tuple]  In the form (wave, flux)
#         line params: [floats]
#             v_cloud
#             b
#             d
#             N
#             f_known
#             lab_lam_0

#     OUTPUT:
#         line model [object instance]
#     """

#     line = KnownVelocityLine(name=name)
#     line.v_cloud = v_cloud
#     line.b = b
#     line.d = d
#     line.N = N
#     line.f = f_known
#     line.lab_lam_0 = lab_lam_0

#     return line


# def createKnownCloud(name, num_lines, lam_0, b, d, N, f_known):
#     """Creates multiple instances of a VoigtAbsorptionLine object 
#     where the Oscillator strength is KNOWN

#     INPUT:
#         name:
#         data:        [tuple]  In the form (wave, flux)
#         num_lines:   [int]
#         line params: [lists of floats]
#             lam_0
#             b
#             d
#             N
#             f_known

#     OUTPUT:
#         line model [object instance]
#     """

#     line0 = createKnownLine(name[0], lam_0[0], b[0], d[0], N[0], f_known[0])
#     cloud = line0

#     if num_lines > 1:
#         for i in range(1, num_lines):
#             line = createKnownLine(name[i], lam_0[i], b[i], d[i], N[i], f_known[i])
#             line.N = line0.N
#             line.b = line0.b

#             cloud *= line

#     return cloud


# def createKnownVelocityCloud(name, num_lines, v_cloud, b, d, N, f_known, lab_lam_0):
#     """Creates multiple instances of a KnownVelocityLine object 
#     where the oscillator strength is KNOWN and the x axis is in velocity space

#     INPUT:
#         name:           [list of strings]
#         data:           [tuple]  In the form (wave, flux)
#         num_lines:      [int]
#         line params:
#             v_cloud     [float]
#             b           [list of floats]
#             d           [list of floats]
#             N           [float]
#             f_known     [list of floats]
#             lab_lam_0   [list of floats]

#     OUTPUT:
#         line model [object instance]
#     """

#     line0 = createKnownVelocityLine(name[0], v_cloud, b, d, N, f_known[0], lab_lam_0[0])
#     cloud = [line0]

#     if num_lines > 1:
#         for i in range(1, num_lines):
#             line = createKnownVelocityLine(
#                 name[i], v_cloud, b, d, N, f_known[i], lab_lam_0[i]
#             )
#             line.v_cloud = line0.v_cloud
#             line.N = line0.N
#             line.b = line0.b
#             line.d = line0.d

#             cloud.append(line)

#     return cloud
