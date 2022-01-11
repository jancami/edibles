import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
import astropy.constants as cst
from scipy.interpolate import interp1d
from scipy.stats import pearsonr, f
import math

import inspect
import collections
from math import floor
from lmfit import Model
from lmfit.models import update_param_vals
from lmfit import Parameter

from edibles.utils.voigt_profile import voigt_absorption_line
from edibles.models import ContinuumModel

from pathlib import Path
from edibles import DATADIR
from edibles import PYTHONDIR

# Known Issues:
# 1. Some Gamma values are missing in the csv table, for now I replace them with 0.0
# 2. When adding a known species, the name match is based on string match and thus not smart enough


# To Do:
# 1. Change the way of indexing clouds, by name (str) only
# 2. Add something to store the info on wavelengths segments, regarding boundaries and resolution
# 3. Function to compare wave_grid to compare to boundaries (so we can decided the resolution to use)

class ROIs():
    def __init__(self):
        self.all_regions = []

    def addROI(self, roi_range, resolution=2.8):
        new_ROI={"Left": np.min(roi_range),
                 "Right": np.max(roi_range),
                 "Resolution": resolution}
        self.all_regions.append(new_ROI)

    def getResolution(self, wave_input):
        grid_left = np.min(wave_input)
        grid_right = np.max(wave_input)
        matching_score = 0.0
        resolution = 2.8

        for region in self.all_regions:
            if grid_left >= region["Right"] or grid_right <= region["Left"]:
                continue

            score_left = np.sort([grid_left, grid_right, region["Right"], region["Left"]])[1]
            score_right = np.sort([grid_left, grid_right, region["Right"], region["Left"]])[2]
            score = (score_right - score_left) / (grid_right - grid_left)
            if score > matching_score:
                resolution = region["Resolution"]
                matching_score = score

        return resolution


class Species():
    # Data structure to hold information of one species. This is used as child data struct of
    # the VelocityComponent class.
    # b, N, and V_off are optional, and is usually passed by VelocityComponent in the **kwargs
    # Species name is mandatory.
    # Code will look up the name in edibles_linelist_atoms.csv and load data
    # User can define his/her target species and line list, by passing name, lambda_0, fjj, and

    def __init__(self, name, lam0=None, fjj=None, Gamma=None, **kwargs):
        # inputs:
        # name: str, name of the species
        # lam0, fjj, Gamma: float, or list/np.array of same length
        #                   wavelengths, fjj, and Gamma of species,
        #                   required when using different data than
        #                   edibles_linelist_atoms.csv
        # kwargs: information on N, b, V_off, fit_status, from Velocity Component

        self.name = name
        self.b = Parameter("b")
        self.N = Parameter("N")

        if None in [lam0, fjj, Gamma]:
            self.__LoadKnownSpecies()
        else:
            self.__LoadNewSpecies(lam0, fjj, Gamma)
        # We have all info for one species, yet we do not know if we need all lines
        self.species_df["In_Range"] = [0] * len(self.species_df)

        self.__parse_b(**kwargs)
        self.__parse_N(**kwargs)

    def __LoadKnownSpecies(self):
        # read in atomic line data frame
        folder = Path(PYTHONDIR + "/data")
        filename = folder / "auxiliary_data" / "line_catalogs" / "edibles_linelist_atoms.csv"
        self.species_df = pd.read_csv(filename)

        # check name is known
        known_species = self.species_df["Species"].to_list()
        # print(np.unique(known_species))
        assert self.name in known_species, "Unknown Species: %s." % self.name

        # load data when species name in known_species
        bool_species_matches = self.species_df["Species"].str.contains(self.name)
        self.species_df = self.species_df.loc[bool_species_matches].reset_index()

        # Also change items with Gamma = NaN to 0.0
        self.species_df["Gamma"] = self.species_df["Gamma"].fillna(0.0)

    def __LoadNewSpecies(self, lam0, fjj, Gamma):
        # Convert parameters to list and assert they have the same length
        para_all = [lam0, fjj, Gamma]
        para_len = []
        for i, item in enumerate(para_all):
            if not (isinstance(item, list) or isinstance(item, np.ndarray)):
                para_all[i] = [item]
            para_len.append(len(item))

        assert_str = "Invalid Input for Species %s. \n"\
                     "lam0, fjj, and Gamma should have same length" % self.name
        assert np.min(para_len) == np.max(para_len), assert_str

        data = np.asarray(para_all).T
        columns = ["WavelengthAir", "OscillatorStrength", "Gamma"]
        self.species_df = pd.DataFrame(data=data, columns=columns)
        self.species_df["Species"] = self.name

    def __parse_b(self, **kwargs):
        if "known_b" in kwargs.keys():
            self.b.set(value=kwargs["known_b"])
            self.b.set(vary=False)
        else:
            if "b" in kwargs.keys():
                self.b.set(value=kwargs["b"])
            else:
                self.b.set(value=1.0)

            if "b_min" in kwargs.keys():
                self.b.set(min=kwargs["b_min"])
            else:
                self.b.set(min=0.1)

            if "b_max" in kwargs.keys():
                self.b.set(max=kwargs["b_max"])
            else:
                self.b.set(max=10.0)

    def __parse_N(self, **kwargs):
        if "known_N" in kwargs.keys():
            self.N.set(value=kwargs["known_N"])
            self.N.set(vary=False)
        else:
            if "N" in kwargs.keys():
                self.N.set(value=kwargs["N"])
            else:
                # set N to -99 and estimate N for each line (to have 3% depth)
                # will update when the lines are filtered
                self.N.set(value=-99)
                N_guess = []
                lambda_0, fjj, Gamma, N, b = self.GetData()
                for lam0, f, gamma in zip(lambda_0, fjj, Gamma):
                    N_guess.append(GuessN(lam0, f, gamma))
                self.species_df["N"] = N_guess

            if "N_min" in kwargs.keys():
                self.N.set(min=kwargs["N_min"])
            else:
                self.N.set(min=0)

            if "N_max" in kwargs.keys():
                self.N.set(max=kwargs["N_max"])
            else:
                self.N.set(max=np.inf) # or do we want a real upper limit?

    def FilterLines(self, wave_regions):
        wave_min, wave_max = np.min(wave_regions), np.max(wave_regions)

        bool_min = self.species_df.WavelengthAir > wave_min
        bool_max = self.species_df.WavelengthAir < wave_max

        self.species_df._set_value(bool_min&bool_max, "In_Range", 1)
        return True

    def GetData(self, ROI=None, margin=5):
        # margin is used to expand the region a little so
        # if a line is just off the edge it will still be
        # included for its possible impact

        # position: lam_0 (sp_df) and v_off (cl)
        # width: Gamma (sp_df) and b (sp_par)
        # strength: fjj (sp_df) and N (sp_par)

        if ROI is None:
            df_out = self.species_df
        else:
            left, right = np.min(ROI) - margin, np.max(ROI) + margin
            df_out = self.species_df.loc[(self.species_df["WavelengthAir"] > left)
                                         & (self.species_df["WavelengthAir"] < right)]

        n_lines = len(df_out)

        return df_out["WavelengthAir"].to_list(), \
               df_out["OscillatorStrength"].to_list(), \
               df_out["Gamma"].to_list(), \
               [self.N.value] * n_lines, [self.b.value] * n_lines

    def __str__(self):
        N_str = "N = %.2E cm-2" % self.N.value
        if self.N.vary:
            N_str = N_str + " (%.2E to %.2E)" % (self.N.min, self.N.max)
        else:
            N_str = N_str + " (fixed)"

        b_str = "b = %.2f km/s" % self.b.value
        if self.b.vary:
            b_str = b_str + " (%.2f to %.2f)" % (self.b.min, self.b.max)
        else:
            b_str = b_str + " (fixed)"

        lambda_0, fjj, Gamma, N, b = self.GetData()
        line_str = "Lines at: " + ", ".join([str(item) for item in lambda_0])

        str_out = self.name.join([" "*4, ":"]) + "\n"
        str_out = str_out + N_str.join([" "*8, "\n"])
        str_out = str_out + b_str.join([" " * 8, "\n"])
        str_out = str_out + line_str.join([" " * 8, "\n"])

        return str_out


class VelocityComponent():

    def __init__(self, name, species=None, **kwargs):
        self.name = name

        self.v_cloud = Parameter("v_cloud")
        self.__ParseV_cloud(kwargs)

        self.species = {}
        if species is not None:
            self.AddSpecies(species)

    def __ParseV_cloud(self, kwargs):
        if "known_v_cloud" in kwargs.keys():
            self.v_cloud.set(value=kwargs["known_v_cloud"])
            self.v_cloud.set(vary=False)
        else:
            if "v_cloud" in kwargs.keys():
                self.v_cloud.set(value=kwargs["v_cloud"])
            else:
                self.v_cloud.set(value=0.0)

            if "v_cloud_min" in kwargs.keys():
                self.v_cloud.set(min=kwargs["v_cloud_min"])
            else:
                self.v_cloud.set(min=-100)

            if "v_cloud_max" in kwargs.keys():
                self.v_cloud.set(max=kwargs["v_cloud_max"])
            else:
                self.v_cloud.set(max=100)

    def AddSpecies(self, species, lam0=None, fjj=None, Gamma=None, **kwargs):
        if isinstance(species, str):
            new_species = Species(species, lam0=lam0, fjj=fjj, Gamma=Gamma, **kwargs)
            self.species[new_species.name] = new_species
        else:
            for item in species:
                self.AddSpecies(item, **kwargs)

    def RemoveSpecies(self, species2remove):
        if isinstance(species2remove, str):
            if species2remove in self.species.keys():
                self.species.pop(species2remove)
        else:
            for item in species2remove:
                self.RemoveSpecies(item)

    def GetData(self, ROI=None):
        # position: lam_0 (sp) and v_off (cl)
        # width: Gamma (sp) and b (sp)
        # strength: fjj (sp) and N (sp)
        lambda_0_cloud, fjj_cloud, Gamma_cloud, N_cloud, b_cloud = [], [], [], [], []
        for key in self.species.keys():
            lambda_0, fjj, Gamma, N, b = self.species[key].GetData(ROI=ROI)
            lambda_0_cloud = lambda_0_cloud + lambda_0
            fjj_cloud = fjj_cloud + fjj
            Gamma_cloud = Gamma_cloud + Gamma
            N_cloud = N_cloud + N
            b_cloud = b_cloud + b

        n_lines = len(lambda_0_cloud)

        return lambda_0_cloud, fjj_cloud, Gamma_cloud, \
               N_cloud, b_cloud, [self.v_cloud.value] * n_lines

    def FilterLines(self, wave_regions):
        for key in self.species.keys():
            print(self.species[key])
            self.species[key].FilterLines(wave_regions)

    def __str__(self):
        # title_str = cloudxxx at xxx km/s (xxx-xxx) or (fixed)
        v_off_str = " at %.2f km/s" % self.v_cloud.value
        if self.v_cloud.vary:
            v_off_str = v_off_str + " (%.2f - %.2f)" % (self.v_cloud.min, self.v_cloud.max)
        else:
            v_off_str = v_off_str + " (fixed)"
        title_str = self.name + v_off_str

        species_str = "%.i species: " % len(self.species.keys())
        species_str = species_str + ", ".join([str(self.species[key].name) for key in self.species.keys()])

        str_out = "\n".join([title_str, species_str] + [str(self.species[key]) for key in self.species.keys()])
        return str_out


class SightLine():
    def __init__(self, name=None):
        if name is None:
            self.name = "SightLine"
        else:
            self.name = name

        self.clouds = {}
        self.ROI = ROIs()

    def AddSpecies(self, species, cloud=None, lam0=None, fjj=None, Gamma=None, **kwargs):
        # use "cloud" to choose which (known) cloud to work with.
        # if not exist, create it

        next_cloud_name = "Cloud_%i" % (len(self.clouds.keys()))
        cloud2add = VelocityComponent(next_cloud_name, **kwargs)
        if cloud is not None:
            for key in self.clouds.keys():
                if self.clouds[key].name == cloud:
                    cloud2add = self.clouds[key]

        cloud2add.AddSpecies(species, lam0=lam0, fjj=fjj, Gamma=Gamma, **kwargs)
        self.clouds[cloud2add.name] = cloud2add

    def AddCloud(self, cloud_name=None, species=None, **kwargs):
        if cloud_name is None:
            cloud_name = "Cloud_%i" % (len(self.clouds.keys()))

        all_names = [self.clouds[key].name for key in self.clouds.keys()]
        assert cloud_name not in all_names, "Cloud name already taken, use another."

        cloud2add = VelocityComponent(cloud_name, species=species, **kwargs)
        self.clouds[cloud_name] = cloud2add

    def AddROI(self, roi_range, resolution=2.8):
        self.ROI.addROI(roi_range=roi_range, resolution=resolution)

    def FilterLines(self, wave_regions):
        for key in self.clouds.keys():
            self.clouds[key].FilterLines(wave_regions)

    def __str__(self):
        out_str = "%i Clouds in " % (len(self.clouds.keys())) + self.name + "\n"
        for key in self.clouds.keys():
            out_str = out_str + str(self.clouds[key])
        return out_str

    def __call__(self, input_wave, plot=True):
        input_wave = np.asarray(input_wave)

        grids = GetWaveGridSegment(input_wave)
        if plot:
            figure = plt.figure(figsize=[4, len(grids)*2], dpi=200)

        for i, grid in enumerate(grids):
            res2use = self.ROI.getResolution(grid)
            lambda_0_all, fjj_all, Gamma_all, N_all, b_all, v_all = [], [], [], [], [], []
            for key in self.clouds.keys():
                pars = self.clouds[key].GetData(ROI=grid)
                lambda_0_all = lambda_0_all + pars[0]
                fjj_all = fjj_all + pars[1]
                Gamma_all = Gamma_all + pars[2]
                N_all = N_all + pars[3]
                b_all = b_all + pars[4]
                v_all = v_all + pars[5]

            flux_grid = voigt_absorption_line(grid,
                                              lambda0=lambda_0_all,
                                              b=b_all,
                                              N=N_all,
                                              f=fjj_all,
                                              v_rad=v_all,
                                              gamma=Gamma_all,
                                              v_resolution=res2use,
                                              n_step=25)

            if plot:
                ax = figure.add_subplot(len(grids), 1, i+1)
                ax.plot(grid, flux_grid)
                ax.tick_params(axis="both", which="major", labelsize=6)
                ax.set_ylabel("Flux", fontsize=8)
                ax.grid()
                if i == len(grids) - 1:
                    ax.set_xlabel("Wavelength (AA)", fontsize=8)

            try:
                flux_out = np.append(flux_out, flux_grid)
            except:
                flux_out = flux_grid

        if plot:
            plt.tight_layout()
            plt.show()

        return(flux_out)


# auxiliary functions
def GuessN(lam0, fjj, Gamma):
    # return estimate N so the central depth of the
    # absorption is 5% of continuum
    FWHM = 2.8 / cst.c.to("km/s").value * lam0
    x_Grid = np.arange(start=-5 * FWHM, stop=5 * FWHM, step=0.005) + lam0
    y = voigt_absorption_line(x_Grid,
                              lambda0=lam0,
                              b=0.5,
                              N=1E10,
                              f=fjj,
                              gamma=Gamma,
                              v_rad=0.0,
                              v_resolution=2.8)

    scale = 0.05 / (1 - np.min(y))
    return scale * 1E10


def GetWaveGridSegment(wave_grid, n_tolerance=10):
    """
    Evaluate if input wave grid is made up of several segments. If so, identify
    these segments and break them up.
    Algorithm: Calculate the step of wave grid, break up when a given step larger
    than n_tolerance * median(step).
    !!! This algorithm does not work if the data is highly segmental!!!

    :param wave_grid: ndarray
           n_tolerance: int, n_tolerance * median(step) is used to identify
                        boundaries of segments. Step of EDIBLES data is usually
                        0.02AA so the default tolerance is 0.2AA for most case.
    :return: grid_segment: list of ndarray
    """

    steps = np.diff(wave_grid)
    tolerance = n_tolerance * np.median(steps)
    break_idx = np.where(steps >= tolerance)[0]
    if len(break_idx) == 0:
        # one piece
        return [wave_grid]
    else:
        return np.split(wave_grid, break_idx+1)



if __name__ == "__main__":
    print("hello sight line!")

    # Known Species (available in the csv file)
    # '6LiI' '7LiI' '85RbI' '87RbI' 'AlI' 'CaI' 'CaII' 'CrI' 'FeI' 'HeI*' 'KI' 'NaI' 'TiII'

    model_builder = SightLine("test_model")

    # Add known species NaI to cloud with v_off = 1.0
    # wavelength, fjj, and Gamma are loaded from csv file, b and N set to default value and range
    # If there is no cloud at 1.0 km/s, it will be created with name Cloud_i
    # v_off for this new cloud is default, with range between (-100, 100) km/s
    model_builder.AddSpecies(species="NaI", cloud="Cloud_0")

    # similarly, add more "known" species to one cloud
    model_builder.AddSpecies(species=["6LiI", "7LiI"], cloud="Cloud_0")  # or set cloud="Cloud_0"

    # add customized species, with name/data different from the table
    model_builder.AddSpecies("Vibranium", cloud="Cloud_0",
                             lam0=[1234.567, 7654.321],
                             fjj=[1.0, 10.0],
                             Gamma=[0, 0],
                             known_b=0.8,
                             N_max=10E13, N_min=10E12)
    # this will fix b_value to 0.8, and N set between 10e12 to 10E13

    # To constrain the v_off of a cloud, you need to add it via .AddCloud method:
    model_builder.AddCloud(cloud_name="Fix_V_Cloud", known_v_cloud=50.0)

    # Can still add multiple species in one time if data is available in csv table and do not constrain b or N
    model_builder.AddSpecies(species=["6LiI", "7LiI"], cloud="Fix_V_Cloud")

    # Can only add one customized species at a time though
    model_builder.AddSpecies(species="Vibranium", cloud="Fix_V_Cloud",
                             lam0=[1234.567, 7654.321],
                             fjj=[1.0, 10.0],
                             Gamma=[0, 0])
    model_builder.AddSpecies(species="NaI", cloud="Fix_V_Cloud", known_N=5 * 10e10)

    # see what we have
    print("=" * 40)
    print(model_builder)
    print("\n" * 5)

    # Add target wavelength range, it will mark lines within the range
    model_builder.FilterLines(np.arange(start=3300, stop=3500, step=0.02))

    # You can add multiple ranges
    model_builder.FilterLines([4000, 4500])
    # The filtering process only checks min and max, so the input here can be the x_grid or its boundary
    # but one range at a time.

    # The output will be slightly different then
    print("=" * 40)
    print(model_builder)


    ######################################
    # Test for Calculation
    ######################################
    print("\n"*10)

    model_builder = SightLine(name="CalcTest")
    model_builder.AddCloud(cloud_name="Cloud_0", v_cloud=-20)
    model_builder.clouds["Cloud_0"].AddSpecies("NaI", N=1E14)

    model_builder.AddSpecies(["NaI", "KI"], cloud="Cloud_1", known_v_cloud=10, N=3E13)

    print(model_builder)

    wave_grid1 = np.arange(start=3301.5, stop=3304, step=0.015)
    wave_grid2 = np.arange(start=5887, stop=5897.5, step=0.015)
    model_builder.AddROI([3300, 3305], resolution=2.8)
    model_builder.AddROI(wave_grid2, resolution=5.0)

    print("\n" * 10)
    flux_out = model_builder(np.append(wave_grid1, wave_grid2), plot=True)

    plt.plot(np.append(wave_grid1, wave_grid2), flux_out)
    plt.show()