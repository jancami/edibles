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

from edibles.utils.voigt_profile import voigt_absorption_line
from edibles.models import ContinuumModel

from pathlib import Path
from edibles import DATADIR
from edibles import PYTHONDIR

class Species():
    # Data structure to hold information of one species. This is used as child data struct of
    # the VelocityComponent class.
    # b, N, and V_off are optional, and is usually passed by VelocityComponent in the **kwargs
    # Species name is mandatory.
    # Code will look up the name in edibles_linelist_atoms.csv and load data
    # User can define his/her target species and line list, by passing name, lambda_0, fjj, and

    # Can be defined by name and the class will lookup in edibles_linelist_atoms.csv
    # Or the user can specify their line, with input lambda_0, fjj, Gamma
    # b, N, and V_off are always optiona
    def __init__(self, name, lam0=None, fjj=None, Gamma=None, **kwargs):
        # inputs:
        # name: str, name of the species
        # lam0, fjj, Gamma: float, or list/np.array of same length
        #                   wavelengths, fjj, and Gamma of species,
        #                   required when using different data than
        #                   edibles_linelist_atoms.csv
        # kwargs: information on N, b, V_off, fit_status, from Velocity Component
        self.name = name

        # read in atomic line data frame
        folder = Path(PYTHONDIR + "/data")
        filename = folder / "auxiliary_data" / "line_catalogs" / "edibles_linelist_atoms.csv"
        self.species_df = pd.read_csv(filename)
        known_species = self.species_df["Species"].to_list()

        if self.name in known_species:
            self.__LoadKnownSpecies()
        else:
            self.__LoadNewSpecies(lam0, fjj, Gamma)

    def __LoadKnownSpecies(self):
        # load data when species name in known_species
        bool_species_matches = self.species_df["Species"].str.contains(self.name)
        self.species_df = self.species_df.loc[bool_species_matches].reset_index()

    def __LoadNewSpecies(self, lam0, fjj, Gamma):
        success = False
        while success is False:
            if lam0 is None or fjj is None or Gamma is None:
                break

            para_all = [lam0, fjj, Gamma]
            for i, item in enumerate(para_all):
                if not (isinstance(item, list) or isinstance(item, np.ndarray)):
                    para_all[i] = [item]

            para_len = [len(item) for item in para_all]
            if not np.min(para_len) == np.max(para_len):
                break

            data = np.asarray(para_all).T
            columns = ["WavelengthAir", "OscillatorStrength", "Gamma"]
            self.species_df = pd.DataFrame(data=data, columns=columns)
            self.species_df["Species"] = self.name
            success = True

        assert success, "Invalid Input for Species %s. " \
                        "Make sure input lam0, fjj, and Gamma " \
                        "have the same length" % self.name

    def FilterLines(self, wave_regions):
        # trim self.species_df according to input wave_regions
        # only keep lines within the wave_regions
        return True

    def GetData(self):
        lambda_0 = self.species_df["WavelengthAir"].to_list()
        fjj = self.species_df["OscillatorStrength"].to_list()
        Gamma = self.species_df["Gamma"].to_list()

        return lambda_0, fjj, Gamma

    def __str__(self):
        # ideally, something like:
        # Na at 3301 and 3302, N=xxx ??

        return "Species %s" % self.name


class VelocityComponent():

    def __init__(self, name, V_off=0.0, species=None):
        self.name = name
        self.V_off = V_off
        self.species = {}
        if species is not None:
            self.AddSpecies(species)

    def AddSpecies(self, species, lam0=None, fjj=None, Gamma=None):
        if isinstance(species, str):
            new_species = Species(species, lam0=lam0, fjj=fjj, Gamma=Gamma)
            self.species[new_species.name] = new_species
        else:
            for item in species:
                self.AddSpecies(item)

    def RemoveSpecies(self, species2remove):
        if isinstance(species2remove, str):
            if species2remove in self.species.keys():
                self.species.pop(species2remove)
        else:
            for item in species2remove:
                self.RemoveSpecies(item)

    def GetData(self):
        return True

    def FilterLines(self):
        return True





if __name__ == "__main__":
    from edibles.utils.voigt_profile import voigt_optical_depth
    lam0 = 3241
    fjj = 0.232
    Gamma = 0
    b = 2.8
    FWHM = 2.8 / cst.c.to("km/s").value * lam0
    x_Grid = np.arange(start=-3*FWHM, stop=3*FWHM, step=0.005) + lam0

    tau_profile = voigt_optical_depth(x_Grid,
                                      lambda0=lam0,
                                      b=b,
                                      N=1E9,
                                      f=fjj,
                                      gamma=0.0, v_rad=0.0)
    flux_profile = np.exp(-1 * tau_profile)
    CD = 1 - np.min(flux_profile)
    scale = 0.03 / CD
    print(np.log10(scale * 1E9))



    print("hello sight line!")
