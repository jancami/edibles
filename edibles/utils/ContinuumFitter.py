import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
import astropy.constants as cst
import math
from scipy.interpolate import CubicSpline

from lmfit import Model
from lmfit.models import update_param_vals
from edibles.utils.voigt_profile import voigt_absorption_line
from edibles.models import ContinuumModel

class ContinuumFitter():
    def __init__(self, wave, flux):
        assert len(wave) == len(flux), "Wave grid and flux must have the same number of elements"
        self.wave = wave
        self.flux = flux

    def SelectPoints(self, n=5, nearest=True, y_message=None, vetoTimeout = False):
        """
        Use interactive backend to select points
        :param n: int, max number of points to be selected
        :param nearest: bool, pick data points from spec that is closest to the selection
        :return: coordinate of selected points, np_array, [[x0, y0], [x1, y1], ... [xn, yn]]
        """

        # an interactive backends will be needed
        import matplotlib
        matplotlib.use('Qt5Agg', force=True)
        import matplotlib.pyplot as tmp_plt

        fig1, ax = tmp_plt.subplots(1, 1)
        ax.plot(self.wave, self.flux, marker=".", linestyle="--", linewidth=0.5)
        ax.grid()
        if y_message is not None:
            ax.set_ylabel(y_message)

        if vetoTimeout:
            timeout = 0
        else:
            timeout = np.median([30, 200, n*5])
        points = tmp_plt.ginput(n, timeout=timeout, mouse_add=1, mouse_pop=3, mouse_stop=2)
        tmp_plt.close()
        points = np.asarray(points)
        # points = [[x0, y0], [x1, y1], [x2, y2], ... [xn, yn]]

        if nearest:
            spec = np.asarray([self.wave, self.flux])
            spec = spec.T
            point_idx = nearest_point(points, spec, scale=True)
            points = spec[point_idx]
            # Make plot to check
            # plt.scatter(self.wave, self.flux, color="0.5", s=10)
            # plt.scatter(self.wave[point_idx], self.flux[point_idx], marker="x", color="r")
            # for point in points:
            #     plt.scatter(point[0], point[1])
            # plt.show()

        return points

    def SplineManualAnchor(self):
    #def SplineManualAnchor(self, n_anchors=5):
        """
        Fit spline continuum from anchor points.
        These anchor points are manually selected by user.
        Up to 99 anchor points
        #:param n_anchors: int, max number max anchor points, user can break in the middle
        :return: continuum, coordinates of anchor points
        """

        # Anchor points do not have to be from the spectrum
        anchor_points = self.SelectPoints(n=99, nearest=False,
                                          y_message="Please Select Anchor Points")
        anchor2use = anchor_points.T
        x_points = anchor2use[0]
        y_points = anchor2use[1]
        spline_continuum = CubicSpline(x_points, y_points)
        return spline_continuum, anchor_points

    def SplineManualRegion(self, n_anchors=5):
        """
        Fit spline continuum from regions selected by users
        Up to 99 regions
        #:param n_regions: int, max number of data regions, user can break in the middle
        :param n_anchors: int, number of anchor points
        :return: continuum, coordinates of anchor points
        """

        # Select continuum regions to fit anc create a "marker" array
        # boundaries have to be from spectrum
        boundary_points = self.SelectPoints(n=99*2, nearest=True,
                                            y_message="Please Select Continuum Regions")
        boundary_x = boundary_points.T[0]
        data2fit = np.zeros_like(self.wave)
        idx = 0
        while idx < len(boundary_x) - 1:
            data2fit[(self.wave >= boundary_x[idx]) & (self.wave <= boundary_x[idx+1])] = 1
            idx = idx + 2

        # tmp check if boundaries work
        # plt.scatter(self.wave, self.flux, s=10, color="0.5")
        # plt.scatter(boundary_points.T[0], boundary_points.T[1], marker="X", color="orange", s=100)
        # plt.scatter(self.wave[data2fit == 1], self.flux[data2fit == 1], s=10, color="r")
        # plt.show()

        # Fit spline continuum model
        wave2fit, flux2fit = self.wave[data2fit == 1], self.flux[data2fit == 1]
        continuum_model = ContinuumModel(n_anchors=n_anchors, verbose=0)
        pars_guess = continuum_model.guess(self.flux, x=self.wave)

        result = continuum_model.fit(data=flux2fit,
                                     params=pars_guess,
                                     x=wave2fit,
                                     weights=np.ones_like(wave2fit))

        # generate continuum and scipy.CubicSpline using the fitted parameters
        params2report = result.params
        x_points, y_points = [], []
        for i in range(n_anchors):
            x_points.append(params2report["x_%i" % (i)].value)
            y_points.append(params2report["y_%i" % (i)].value)

        spline_continuum = CubicSpline(np.asarray(x_points), np.asarray(y_points))
        anchor_points = np.asarray([x_points, y_points]).T

        return spline_continuum, anchor_points



def nearest_point(points, spec, scale=True):
    """
    Find points on the spectrum that is closest to input coordinates
    :param points: coordinates of input points, [[x0, y0], [x1, y1], ... [xn, yn]]
    :param spec: coordinates of spectral data to be compared with
    :param scale: bool, if true, coordinates will be scaled to [0, 1], default to true
    :type scal: bool

    :return: point_idx, the index of points in the spec that are most close to the input points
    :rtype: int array
    """
    from scipy.spatial.distance import cdist

    #points = np.asarray(points)
    spec = np.asarray(spec)
    assert len(points.shape) == 2, "Inputs must have 2 dimension"
    assert len(spec.shape) == 2, "Inputs must have 2 dimension"

    if scale:
        points, spec = points.T, spec.T
        xmax, xmin = np.max(spec[0]), np.min(spec[0])
        ymax, ymin = np.max(spec[1]), np.min(spec[1])
        points = points - np.array([[xmin] * len(points[0]), [ymin]*len(points[0])])
        points = points / np.array([[xmax - xmin] * len(points[0]), [ymax - ymin] * len(points[0])])
        spec = spec - np.array([[xmin] * len(spec[0]), [ymin]*len(spec[0])])
        spec = spec / np.array([[xmax - xmin] * len(spec[0]), [ymax - ymin] * len(spec[0])])
        points, spec = points.T, spec.T

    d_matrix = cdist(points, spec, 'euclidean')

    point_idx = []
    for i in range(points.shape[0]):
        d_matrix_point = d_matrix[i]
        point_idx.append(np.where(d_matrix_point == np.min(d_matrix_point))[0][0])

    return np.asarray(point_idx)


if __name__ == "__main__":
    from edibles.utils.edibles_oracle import EdiblesOracle
    from edibles.utils.edibles_spectrum import EdiblesSpectrum


    def make_test_plot(tester, continuum, anchor):
        fig = plt.figure(figsize=(10, 6.5))
        plt.gcf().subplots_adjust(hspace=0)
        spec = gridspec.GridSpec(ncols=1, nrows=2,
                                 height_ratios=[4, 4])

        # Top panel for raw data and overall fitting
        ax0 = fig.add_subplot(spec[0])
        plt.gca().xaxis.set_visible(False)
        plt.step(tester.wave, tester.flux, color="0.5")
        plt.scatter(anchor.T[0], anchor.T[1], marker="x", s=80, color="r")
        plt.plot(tester.wave, continuum(tester.wave), color="orange")
        plt.ylabel("Raw Data")

        # Lower panel for normalized data and multi components
        ax1 = fig.add_subplot(spec[1])
        plt.step(testTester.wave, testTester.flux / continuum(testTester.wave),color="0.5")
        plt.scatter(anchor.T[0], np.ones_like(anchor.T[1]), marker="x", s=80, color="r")
        plt.plot(testTester.wave, np.ones_like(testTester.wave), linestyle="--", color="orange")
        plt.ylabel("Normalized Data")
        plt.xlabel('Wavelenght $\AA$')
        plt.show()


    pythia = EdiblesOracle()
    List = pythia.getFilteredObsList(object=["HD 147889"], OrdersOnly=True, Wave=6707.0)
    file_all = List.values.tolist()
    for filename in file_all:
        sp = EdiblesSpectrum(filename)
        wave, flux = sp.bary_wave, sp.flux
        wave, flux = wave[(wave > 6706) & (wave < 6711)], flux[(wave > 6706) & (wave < 6711)]
        testTester = ContinuumFitter(wave, flux)

        cont, anchor = testTester.SplineManualAnchor()
        make_test_plot(testTester, cont, anchor)

        cont, anchor = testTester.SplineManualRegion(n_anchors=6)
        make_test_plot(testTester, cont, anchor)




