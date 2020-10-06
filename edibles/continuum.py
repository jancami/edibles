import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import csv

from edibles.utils.edibles_spectrum import EdiblesSpectrum
from edibles.models import ContinuumModel


class Continuum:
    """A class that has multiple methods for fitting different types of continuua.

    Args:
        Spectrum (EdiblesSpectrum): The input EiblesSpectrum data
        plot (bool): If true, plots the continuum fit
        verbose (int): If > 0, display more status messages

    """

    def __init__(self, Spectrum, method="spline", plot=False, verbose=0, *args, **kwargs):

        # check existing the available continuum csv files
        try:
            # if Spectrum.continuum_filename:
            num_saved_continuua = 0
            with open("/home/kulik/python/ediblesdr4/DR4/continuum/HD23466/BLUE_346/HD23466_w346_blue_20180731_O11.csv") as f:
                # reader = csv.reader(csvfile, delimiter=' ', quotechar='|')
                for row in f:
                    print(row)
                    if len(row) > 0:
                        if row[0] == '######':
                            num_saved_continuua += 1
        except AttributeError:
            print('No previously saved data')


        self.method = method
        self.Spectrum = Spectrum
        self.plot = plot
        self.verbose = verbose

        if method == "spline":
            self.spline(*args, **kwargs)
        elif method == "alphashape":
            self.alphashape(*args, **kwargs)
        elif method == "polynomial":
            self.polynomial(*args, **kwargs)

    def spline(self, *args, **kwargs):
        """A spline function through a set number of anchor points

        Args:
            n_anchors (int): The number of anchor points in the spline

        """

        if self.verbose > 0:
            print("method: ", self.method)
            print()

        n_anchors = kwargs["n_anchors"]

        self.model = ContinuumModel(n_anchors=n_anchors)
        cont_pars = self.model.guess(self.Spectrum.flux, x=self.Spectrum.wave)

        self.result = self.model.fit(
            data=self.Spectrum.flux, params=cont_pars, x=self.Spectrum.wave
        )

        if self.plot:
            self.result.plot_fit()
            plt.show()

    def alphashape(self):
        if self.verbose > 0:
            print("method: ", self.method)
            print()

        print("This method is not available yet.")

    def polynomial(self):
        if self.verbose > 0:
            print("method: ", self.method)
            print()
        print("This method is not available yet.")


    def add_to_csv(self, user, comments=False):

        # Tests not in testing folder beceause we dont want to write the testing data
        assert isinstance(cont.model, ContinuumModel)
        assert isinstance(user, str)
        assert len(user) > 0, "A name must be entered"
        assert isinstance(comments, str)
        assert isinstance(cont.model.n_anchors, int)
        assert isinstance(datetime.now(), datetime)


        csv_file = self.Spectrum.filename.replace(".fits", ".csv").replace(
            "/DR4/data/", "/DR4/continuum/"
        )

        line1 = "method=" + self.method + ", " + "n_anchors=" + str(self.model.n_anchors)
        line2 = (
            "datetime=" + str(datetime.now()) + ", "
            + "user=" + user + ", " + "Comments: " + comments
        )

        x_points = [self.result.params[xname].value for xname in self.model.xnames]
        y_points = [self.result.params[yname].value for yname in self.model.ynames]

        header = line1 + "\n" + line2
        with open(csv_file, mode="a") as f:
            np.savetxt(f, (x_points, y_points), delimiter=",", header=header, comments="# ")
            f.write("\n")

        if self.verbose > 0:
            print("Appended to file!")
        if self.verbose > 1:
            print("File appended to: " + csv_file)


if __name__ == "__main__":

    sp = EdiblesSpectrum("/HD23466/BLUE_346/HD23466_w346_blue_20180731_O11.fits")
    subset = sp.getSpectrum(xmin=3270, xmax=3305)

    cont = Continuum(sp, method="spline", n_anchors=5, plot=True, verbose=2)

    print("X names: ", cont.model.xnames)
    print("X values: ", [cont.result.params[param].value for param in cont.model.xnames])
    print("Y names: ", cont.model.ynames)
    print("Y values: ", [cont.result.params[param].value for param in cont.model.ynames])


    cont.add_to_csv(user="First Last", comments="These are test points and should not be used.")

