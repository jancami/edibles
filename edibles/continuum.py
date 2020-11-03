import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
from lmfit.models import update_param_vals
from pprint import pprint

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

        # check existing the available continuum csv files
        try:
            if Spectrum.continuum_filename:
                self.num_saved_continua = 0
                with open(Spectrum.continuum_filename) as f:
                    for row in f:
                        if "######" in row:
                            self.num_saved_continua += 1
        except AttributeError:
            print("No previously saved data")

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
        """A function that uses alphashape to find a continuum.

        Note:
            Currently not implemented

        """

        if self.verbose > 0:
            print("method: ", self.method)
            print()

        print("This method is not available yet.")

    def polynomial(self):
        """A function that uses a polynomial to find a continuum.

        Note:
            Currently not implemented

        """

        if self.verbose > 0:
            print("method: ", self.method)
            print()
        print("This method is not available yet.")

    def prebuilt_model(self, chosen_save_num=None, plot=False, verbose=0):
        """A function that generates continua based on data saved in csv files.


        Args:
            chosen_save_num (int): The 'save number' of the continuum data, default=None.
                If None, the function will create all saved models and (possibly) plot them.
            plot (bool): If True, plot the model(s) once it is created
            verbose (int): If > 0, print more information about the data

        """

        # assert self.num_saved_continua is not 0, otherwise there is no known continuum point
        assert self.num_saved_continua > 0, "There is no saved continuum."

        # read and parse file contents
        saves_counter = 0
        saves_dict = {}
        with open(self.Spectrum.continuum_filename, mode="r") as f:
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
                            saves_dict[name]["x"] = [
                                float(item) for item in line.split(",")
                            ]
                        else:
                            saves_dict[name]["y"] = [
                                float(item) for item in line.split(",")
                            ]

        if chosen_save_num is not None:
            assert chosen_save_num < self.num_saved_continua, (
                "There are only " + str(self.num_saved_continua) + " saved continua."
            )
            chosen_save = saves_dict["save" + str(chosen_save_num)]

            if verbose > 0:
                print("Number of saved continuum datasets: ", saves_counter)
                print("Save chosen: save" + str(chosen_save_num))
                pprint(chosen_save)

            cont_model = ContinuumModel(n_anchors=chosen_save['n_anchors'])
            params = cont_model.make_params()
            for i in range(cont_model.n_anchors):
                params['%sx_%i' % (cont_model.prefix, i)].set(value=chosen_save["x"][i], vary=False)
                params['%sy_%i' % (cont_model.prefix, i)].set(value=chosen_save["y"][i], vary=False)

            params = update_param_vals(params, cont_model.prefix)

            out = cont_model.eval(params=params, x=self.Spectrum.wave)

            if plot:
                plt.plot(self.Spectrum.wave, self.Spectrum.flux)
                plt.plot(self.Spectrum.wave, out)
                plt.scatter(chosen_save["x"], chosen_save["y"], marker="x", s=80, color="k")
                plt.show()

            return out
        else:
            for j in range(saves_counter):
                chosen_save = saves_dict["save" + str(j)]
                try:
                    method = str(chosen_save["method"])
                except:
                    method = "unknown"

                try:
                    n_anchors = str(chosen_save["n_anchors"])
                except:
                    n_anchors = "unknown"

                try:
                    date_time = str(chosen_save["datetime"])
                except:
                    date_time = "unknown"

                try:
                    user = str(chosen_save["user"])
                except:
                    user = "unknown"

                try:
                    comments = str(chosen_save["comments"])
                except:
                    comments = "None"

                print("Continuum group {a}/{b}".format(a = j, b=saves_counter-1))
                print("Method: " + method + ", n_anchors: " + n_anchors + ";")
                print("Added by: " + user + ", at time:" + date_time + ";")
                print("Comments: " + comments)

                if plot:
                    cont_model = ContinuumModel(n_anchors=chosen_save['n_anchors'])
                    params = cont_model.make_params()
                    for i in range(cont_model.n_anchors):
                        params['%sx_%i' % (cont_model.prefix, i)].set(value=chosen_save["x"][i], vary=False)
                        params['%sy_%i' % (cont_model.prefix, i)].set(value=chosen_save["y"][i], vary=False)

                    params = update_param_vals(params, cont_model.prefix)

                    out = cont_model.eval(params=params, x=self.Spectrum.wave)

                    plt.plot(self.Spectrum.wave, self.Spectrum.flux)
                    plt.plot(self.Spectrum.wave, out)
                    plt.scatter(chosen_save["x"], chosen_save["y"], marker="x", s=80, color="k")
                    plt.show()
                print("Please make your selection back in the script.")

    def add_to_csv(self, user, comments):
        """A function that saves the continuum model parameters to a csv file.

        Each save appears as follows:

        | ######
        | # method=spline
        | # n_anchors=4
        | # datetime=2020-10-06 10:56:02.192865
        | # user=First Last
        | # comments=This is a comment.
        | x1, x2, x3, x4
        | y1, y2, y3, y4


        Args:
            user (str): The name of the person adding the data
            comments (str): Any comments the user wishes to make about the data to be saved

        Note:
            The data is saved to the ediblesdr4 github repository,
            with the same filepath as the original FITS file.



        """

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

    # build a 4 anchor points spline
    cont = Continuum(sp, method="spline", n_anchors=4, plot=False, verbose=2)
    # Guess the model parameters
    params = cont.model.guess(sp.flux, x=sp.wave)
    # Fit the model
    result = cont.model.fit(data=sp.flux, params=params, x=sp.wave)
    # Get the output of the fit model
    out = result.eval(params=result.params, x=sp.wave)
    # Print the result parameters
    print(result.params)
    # Plot
    plt.plot(sp.wave, sp.flux)
    plt.plot(sp.wave, out)
    plt.show()
    cont.add_to_csv(
        user="Mario", comments="Test of 4 anchor points spline"
    )


    # build a 8 anchor points spline
    cont = Continuum(sp, method="spline", n_anchors=8, plot=False, verbose=2)
    # Guess the model parameters
    params = cont.model.guess(sp.flux, x=sp.wave)
    # Fit the model
    result = cont.model.fit(data=sp.flux, params=params, x=sp.wave)
    # Get the output of the fit model
    out = result.eval(params=result.params, x=sp.wave)
    # Plot
    plt.plot(sp.wave, sp.flux)
    plt.plot(sp.wave, out)
    plt.show()
    cont.add_to_csv(
        user="Luigi", comments="Test of 8 anchor points spline"
    )

    # reinitialize the edibles spectrum class, to get the continuum_filename in it
    # sp = EdiblesSpectrum("/HD23466/BLUE_346/HD23466_w346_blue_20180731_O11.fits")
    # sp.getSpectrum(xmin=3270, xmax=3305)
    cont2 = Continuum(sp)
    cont2.prebuilt_model(chosen_save_num=None, plot=True, verbose=1)
