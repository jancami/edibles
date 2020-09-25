import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
from scipy.interpolate import CubicSpline
from lmfit.models import update_param_vals


from edibles.utils.edibles_spectrum import EdiblesSpectrum
from edibles.models import ContinuumModel


class Continuum:
    """A class that has multiple methods for fitting different types of continua.

    Args:
        Spectrum (EdiblesSpectrum): The input EiblesSpectrum data
        method (str): The method of fitting
        plot (bool): If true, plots the continuum fit
        verbose (int): If > 0, display more status messages

    """


    def __init__(self, Spectrum, method="None", plot=False, verbose=0, *args, **kwargs):

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
        self.cont_pars = self.model.guess(self.Spectrum.flux, x=self.Spectrum.wave)

        self.result = self.model.fit(
            data=self.Spectrum.flux, params=self.cont_pars, x=self.Spectrum.wave
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


    def prebuilt_model(self, chosen_save_num=0, plot=False, verbose=0):

        # read and parse file contents
        csv_file = self.Spectrum.continuum_filename

        saves_counter = 0
        saves_dict = {}
        with open(csv_file, mode="r") as f:
            for line in f:
                line = line.split("\n")[0]

                # initialize new save group
                if len(line) > 0:
                    if line == "######":
                        name = "save" + str(saves_counter)
                        saves_counter += 1
                        saves_dict[name] = {"x": None, "y": None}

                    # update dict
                    elif line[0:2] == "# ":
                        key, val = line.split("# ")[1].split("=")

                        if key == "n_anchors":
                            val = int(val)
                        if key == "datetime":
                            val = datetime.strptime(val, "%Y-%m-%d %H:%M:%S.%f")

                        saves_dict[name][key] = val

                    else:
                        if saves_dict[name]["x"] is None:
                            saves_dict[name]["x"] = [float(item) for item in line.split(',')]
                        else:
                            saves_dict[name]["y"] = [float(item) for item in line.split(',')]


        if verbose > 0:
            print("Number of saved continuum datasets: ", saves_counter)


        chosen_save = saves_dict["save" + str(chosen_save_num)]

        print(chosen_save)


        # Create spline
        spline = CubicSpline(chosen_save['x'], chosen_save['y'])
        out = spline(self.Spectrum.wave)

        # TODO: Implement this with a ContinuumModel - almost done, still buggy

        # cont = ContinuumModel(n_anchors=chosen_save['n_anchors'])

        # params = cont.make_params()
        # for i in range(cont.n_anchors):
        #     # params.add(name="x_" + str(i), )
        #     params['%sx_%i' % (cont.prefix, i)].set(value=chosen_save["x"][i], vary=False)
        #     params['%sy_%i' % (cont.prefix, i)].set(value=chosen_save["y"][i], vary=False)
        # update_param_vals(params, cont.prefix)

        # print(cont.param_names)
        # print([params["x_" + str(i)].value for i in range(cont.n_anchors)])
        # print([params["y_" + str(i)].value for i in range(cont.n_anchors)])

        # out = cont.eval(params=params, x=self.Spectrum.wave)

        if plot:
            plt.plot(self.Spectrum.wave, self.Spectrum.flux)
            plt.plot(self.Spectrum.wave, out)
            plt.scatter(chosen_save['x'], chosen_save['y'], marker='x', s=80, color='k')
            plt.show()


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

        x_points = [self.result.params[xname].value for xname in self.model.xnames]
        y_points = [self.result.params[yname].value for yname in self.model.ynames]

        with open(csv_file, mode="a") as f:
            f.write("######\n")
            f.write("# method=" + str(self.method) + "\n")
            f.write("# n_anchors=" + str(self.model.n_anchors) + "\n")
            f.write("# datetime=" + str(datetime.now()) + "\n")
            f.write("# user=" + str(user) + "\n")
            f.write("# comments=" + str(comments) + "\n")
            np.savetxt(f, (x_points, y_points), delimiter=",")
            f.write("\n")

        if self.verbose > 0:
            print("Appended to file!")
        if self.verbose > 1:
            print("File appended to: " + csv_file)


if __name__ == "__main__":

    sp = EdiblesSpectrum("/HD23466/BLUE_346/HD23466_w346_blue_20180731_O11.fits")
    sp.getSpectrum(xmin=3270, xmax=3305)

    # cont = Continuum(sp, method="spline", n_anchors=5, plot=False, verbose=2)

    # print("X names: ", cont.model.xnames)
    # print("X values: ", [cont.result.params[param].value for param in cont.model.xnames])
    # print("Y names: ", cont.model.ynames)
    # print("Y values: ", [cont.result.params[param].value for param in cont.model.ynames])


    # cont.add_to_csv(user="First Last", comments="These are test points and should not be used.")

    cont = Continuum(sp).prebuilt_model(chosen_save_num=1, plot=True, verbose=1)
