from edibles.edibles.utils import generate_continuum
from edibles.edibles.models import Cont1D


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
