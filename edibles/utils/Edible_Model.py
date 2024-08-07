"""
This script is part of a summer project that focuses on modeling the absorption profiles
of molecular and atomic species using Voigt profiles. The class `EdiblesModeling` is designed
to fit the absorption features in spectra obtained from the EDIBLES database.

*** see the examples at the end ***
"""

import glob as g
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from edibles.utils.voigt_profile import *

import io
import contextlib
import sys
import os
from edibles.utils.edibles_oracle import EdiblesOracle
from edibles.utils.edibles_spectrum import EdiblesSpectrum



@contextlib.contextmanager
def suppress_stdout():
    """
    Context manager to suppress the standard output.
    Useful to silence print statements from external libraries.
    """
    with open(os.devnull, 'w') as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout



class EdiblesModeling:
    """
    Class to perform model fitting on Edibles spectra for any molecular or atomic species.
    """


    def __init__(self, obj_name, ext, wrange, isotope, f, gamma, b, N, v_rad, v_resolution, zoom=None):
        """
        Initialize the EdiblesModelFit class with necessary parameters.

        Args:
        obj_name (str): Name of the target star
        ext : give your file extension # e,g 14,15, # incase you don't know the ext then you can use "" 
        wrange (float64): The wavelength range of the observation (in Angstroms)
        isotope (float64): Central (rest) wavelength for the absorption line, in Angstrom
        f (float64): The oscillator strength (dimensionless)
        gamma (float64):  Lorentzian gamma (=HWHM) component
        b (float64): The b parameter (Gaussian width), in km/s.
        N (float64): The column density (in cm^{-2})
        v_rad (float64): Radial velocity of absorption line (in km/s)
        v_resolution (float64): Instrument resolution in velocity space (in km/s)
        zoom (float64): The zoom range of the observation (in Angstroms)
        """
        
        
        self.obj_name = obj_name
        self.ext = ext
        self.wrange = wrange
        self.isotope = isotope
        self.f = f
        self.gamma = gamma
        self.b = [float(b)]
        self.N = N
        self.v_rad = v_rad
        self.v_resolution = v_resolution
        self.zoom = zoom if zoom else wrange

    def model_fit(self):
        """
        Identify the files containing absorption features and fit the absorption profile to 
        a Voigt profile for any molecular or atomic species. Generates and displays plots of the model fit.
        """
        

        # Suppress output from external libraries
        with suppress_stdout():
            pythia = EdiblesOracle()
            obs_list_all_files= pythia.getFilteredObsList(
                object=[self.obj_name], MergedOnly=False, Wave=(self.wrange[0] + self.wrange[1]) / 2
            )
        obs_list_all_files = obs_list_all_files.values.tolist()

        with suppress_stdout():
            pythia = EdiblesOracle()
            obs_list_marged_only= pythia.getFilteredObsList(
                object=[self.obj_name], MergedOnly=True, Wave=(self.wrange[0] + self.wrange[1]) / 2
            )
        obs_list_marged_only = obs_list_marged_only.values.tolist()
        obs_list = [obs for obs in obs_list_all_files if obs not in obs_list_marged_only]

        # Filter files based on extensions
        files = obs_list
        files = [filename for filename in files if any(filename.endswith(f'{e}.fits') for e in self.ext)]
        

        print(f"{150 * '='}\n{f'{self.obj_name}'.center(100)}\n{150 * '='}")
        
        print(" ")
        
        print(f'The total number of searched files for {self.obj_name} is : {len(files)}')
        for i in files:
            print(i)
        print("\n ")

        # Loop through each file to generate plots
        for i in files:
            index = files.index(i)
            sp = EdiblesSpectrum(i)
            sp.getSpectrum(self.wrange[0], self.wrange[1])
            wave = sp.bary_wave
            flux = sp.bary_flux
            idx = np.where((wave > self.wrange[0]) & (wave < self.wrange[1]))
            wave = wave[idx]
            flux = flux[idx]
            flux = flux / np.median(flux)
            absorption_lines_dict = {}

            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 4), dpi=100)

            print(f"plot of {i} is done")

            # Calculate and plot individual absorption lines
            for j in range(len(self.isotope)):
                AbsorptionLin = voigt_absorption_line(
                    wave,
                    lambda0=self.isotope[j],
                    b=self.b,
                    N=self.N[j],
                    f=self.f[j],
                    gamma=self.gamma[j],
                    v_rad=self.v_rad,
                    v_resolution=self.v_resolution,
                )

                key_name = f'AbsorptionLin_{j+1}'
                absorption_lines_dict[key_name] = AbsorptionLin

                ax1.plot(wave, AbsorptionLin, marker="*", label=key_name)

            # Plot original flux and individual absorption lines
            ax1.plot(wave, flux, color="blue", label="Original flux")
            ax1.set_xticks(np.arange(self.wrange[0], self.wrange[1], 0.7))
            ax1.set_xlim(self.zoom[0], self.zoom[1])
            ax1.set_xlabel("Wavelength ($\AA$)")
            ax1.set_ylabel("Normalised flux")
            ax1.set_title("Individual Absorption Lines")
            ax1.legend()

            # Calculate and plot the overall absorption profile
            absorption_product = np.ones_like(wave)
            for key in absorption_lines_dict:
                absorption_product *= absorption_lines_dict[key]

            ax2.plot(wave, flux, color="blue", label="Original flux")
            ax2.plot(wave, absorption_product, color="red", marker="o", label="model")
            ax2.set_xticks(np.arange(self.wrange[0], self.wrange[1], 0.7))
            ax2.set_xlim(self.zoom[0], self.zoom[1])
            ax2.set_xlabel("Wavelength ($\AA$)")
            ax2.set_ylabel("Normalised flux")
            ax2.set_title("Overall Profile")
            ax2.legend()

            # Display the plots
            fig.suptitle(
                f'{i[1:9]} absorption lines --> file {index}',
                color='red',
                fontsize=12,
                fontweight='bold',
                y=1.05
            )
            plt.show()

            print("\n")




''''
In case you don't know whats is the file name ends with (for ext parameter)then you can just give : `ext=[""]` this 
will give you all the file for the target. Thus in the output you will be able to see all the file corresponding to your target.
After getting the file names you may get any error in the plotting section saying file is not found: this is because some files 
are merged files which are not there in DR4, so if you are interested in those files downloade those files from 
edibles merged data base.
'''




# Example usage:
# ============================================
# example1: Li (2 doublet- 7Li & 6Li)
# ============================================

# obj = "HD 147084"
# ext = [1]
# wrange = [6705, 6710]
# isotope = [[6707.761, 6707.912], [6707.921, 6708.072]]
# f = [[0.4982, 0.2491], [0.4982, 0.2491]]
# gamma = [[3.69e07, 3.69e07], [3.69e07, 3.69e07]]
# b = 3
# N = [[0.001e13], [0.001e13 / 5]]
# v_rad = 11.5
# v_resolution = 3
# zoom = [6705, 6710]

# model = EdiblesModeling(obj, ext, wrange, isotope, f, gamma, b, N, v_rad, v_resolution, zoom)
# model.model_fit()



# ============================================
# example2: CN (2 singlet & 1 doublet of 12CN)
# ============================================

# obj = "HD 149757"
# ext = [5]
# wrange = [3870, 3880]
# isotope = [[3873.9939], [3874.6024, 3874.6059], [3875.759]]
# f = [[0.022800], [0.034200, 0.039200], [0.011400]]
# gamma = [[1.610e07], [1.610e07, 1.610e07], [1.610e07]]
# b = 1.85
# N = [[0.075e13], [0.075e13], [0.075e13]]
# v_rad = -15.35
# v_resolution = 1.54748141144159
# zoom = [3873, 3876]

# model = EdiblesModeling(obj, ext, wrange, isotope, f, gamma, b, N, v_rad, v_resolution, zoom)
# model.model_fit()



# ============================================
# example3: CH+ (2 singlet of 12CH+)
# ============================================

# obj = "HD 149757"
# ext = [5]
# wrange = [4231.8, 4233]
# isotope = [[4232.288], [4232.548]]
# f = [[0.003], [0.0054]]
# gamma = [[1e8], [1e8]]
# b = 1
# N = [[0.47e14 / 80], [0.47e14]]
# v_rad = -15
# v_resolution = 5.75

# model = EdiblesModeling(obj, ext, wrange, isotope, f, gamma, b, N, v_rad, v_resolution)
# model.model_fit()





'''
#### To run all the sightlins togather

import os

parent_folder = "./DR4"  # change this to your data file location

subfolders = [f.name for f in os.scandir(parent_folder) if f.is_dir()]

SUBNAMES = [''.join(filter(str.isdigit, name)) for name in subfolders]
Target=[]
for i in range(len(SUBNAMES)):
    target = "HD " + SUBNAMES[i]
    Target.append(target)# CN lines

    

for i in Target:

    obj = i   
    ext= [5]
    wrange = [3870,3880] # this data are corresponding to CN but you can change it whatever element you like 
    isotope = [[3873.9939], [3874.6024, 3874.6059],[3875.759]]
    f = [[0.022800],[0.034200, 0.039200], [0.011400]]
    gamma = [[1.610e07],[1.610e07, 1.610e07], [1.610e07]]
    b=1.85
    N = [[0.075e13],[0.075e13],[0.075e13]]
    v_rad = -15.35
    v_resolution = 1.54748141144159
    zoom = [3873,3876]

    model = EdiblesModeling(obj, ext, wrange, isotope, f, gamma, b, N, v_rad, v_resolution, zoom)
    model.model_fit()





'''
