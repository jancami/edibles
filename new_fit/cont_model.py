import numpy as np

from sherpa.models.model import ArithmeticModel
from sherpa.models.parameter import Parameter

__all__ = ('Cont1D', )


from scipy.interpolate import CubicSpline


def _cont1d(pars, x):
    '''
    This function fits a continuum to data separated into n sections
    where the x and y-values are the median of each section using a cubic spline

    INPUT:
    x:            [ndarray]                   wavelength grid (angstroms)

    pars:         [list]                      (y_pts, delta_v, n_piece)
        y1 - y8:      [floats]                    input y_points to fit spline
        n_piece:      [int, default=2]            number of sections to split data into

    OUTPUT:
    y_spline:     [ndarray]                   continuum flux value array
    '''

    (y1, y2, y3, y4, y5, y6 ,y7, y8, n_piece) = pars

    if type(n_piece) is not int:
        n_piece = int(n_piece)

    # split x & y arrays into n_piece*2 sections
    x_sections = np.array_split(x, n_piece*2)
    # y_sections = np.array_split(flux, n_piece*2)

    # initialize list of points to spline fit
    x_points = [np.min(x)]
    # y_points = [np.median(y_sections[0])]


    # loop through every other section (1, 3, 5...)
    # make n_piece+1 points to fit a spline through
    # create a spline point on edge of each piece
    for i in range(1, len(x_sections), 2):
        # set x_point 
        x_point = np.max(x_sections[i])

        if i == range(len(x_sections))[-1]:
            # span = y_sections[i]
            # y_point = np.median(span)
            x_point = np.max(x)

        # else:
            # span = np.append(y_sections[i], y_sections[i+1])
            # y_point = np.median(span)

        x_points.append(x_point)
        # y_points.append(y_point)



    y_pts = [y1, y2, y3, y4, y5, y6, y7, y8]
    y_points = []
    for i in range(n_piece+1):
        if y_pts[i] is not None:
            y_points.append(y_pts[i])
        else:
            break


    spline = CubicSpline(x_points, y_points)
    y_spline = spline(x)



    return y_spline


class Cont1D(ArithmeticModel):
    """A one-dimensional continuum spline.

    The model parameters are:

    y_pts
        initial y_points to fit
    n_points
        number of segments to break spectrum into.
    

    """

    def __init__(self, name='cont1d'):

        # self.flux    = Parameter(name, 'flux', None, alwaysfrozen=True, hidden=True)
        # self.y_pts   = Parameter(name, 'y_pts', None)
        self.y1 = Parameter(name, 'y1', 1.0, frozen=True)
        self.y2 = Parameter(name, 'y2', 1.0, frozen=True)
        self.y3 = Parameter(name, 'y3', 1.0, frozen=True)
        self.y4 = Parameter(name, 'y4', None, frozen=True)
        self.y5 = Parameter(name, 'y5', None, frozen=True)
        self.y6 = Parameter(name, 'y6', None, frozen=True)
        self.y7 = Parameter(name, 'y7', None, frozen=True)
        self.y8 = Parameter(name, 'y8', None, frozen=True)

        self.n_piece = Parameter(name, 'n_piece', 2, min=2, max=8, frozen=True)

        ArithmeticModel.__init__(self, name,
            (self.y1, self.y2, self.y3, self.y4, self.y5, self.y6, 
                self.y7, self.y8, self.n_piece))

    def calc(self, pars, x, *args, **kwargs):
        """Evaluate the model"""

        # If given an integrated data set, use the center of the bin
        

        return _cont1d(pars, x)