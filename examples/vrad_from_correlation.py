import numpy as np
import matplotlib.pyplot as plt
from edibles import PYTHONDIR
from edibles.utils.edibles_oracle import EdiblesOracle
from edibles.utils.edibles_spectrum import EdiblesSpectrum
from edibles.utils.voigt_profile import *
from pathlib import Path
import astropy.constants as cst
from scipy.interpolate import interp1d
from scipy.stats import pearsonr


def determine_vrad_from_correlation(wave, flux, model):
    """
    Function to calculate the correlation between an observed spectrum and a model as a function of
    radial velocity and return the radial velocity with the highest correlation coefficient. 

    Args:
        wave (float64): array of wavelengths
        flux (float64): Flux (observed)
        model(float64): model

    Returns:
        vrad_best: radial velocity corresponding to highest correlation. 

        """
    # Create the grid of velocities at which to calculate the correlation. 
    # Using a step size of 0.1 km/s here, and a range of -50 to 50; this should 
    # suffice for most sightlines. 
    v_rad_grid = np.arange(-50.,50.,.1) # in km / s
    all_corr = v_rad_grid * 0.
    #print(v_rad_grid)
    for loop in range(len(v_rad_grid)):
        v_rad = v_rad_grid[loop]
        Doppler_factor = 1. + v_rad / cst.c.to("km/s").value
        #print(Doppler_factor)
        new_wave = wave * Doppler_factor
        # Interpolate shifted model to original wavelength grid
        interpolationfunction = interp1d(
        new_wave, model, kind="cubic", fill_value="extrapolate")
        interpolatedModel = interpolationfunction(wave)
        # Calculate correlation coefficient
        this_c, _ = pearsonr(flux, interpolatedModel)
        all_corr[loop] = this_c
        #print(v_rad, this_c)
        #plt.plot(wave, flux)
        #plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
        #plt.plot(wave, interpolatedModel, color="orange", marker="*")
        #plt.show()
    # Return the radial velocity at the maximum correlation.  
    v_rad_best = v_rad_grid[np.argmax(all_corr)]
    print("Highest correlation for v_rad = ",v_rad_best)
    plt.plot(v_rad_grid, all_corr, marker='*')
    plt.xlabel("v_rad [km/s]")
    plt.ylabel("Correlation coefficient")
    plt.axvline(v_rad_best, color='red')
    plt.show()
    Doppler_factor = 1. + v_rad_best / cst.c.to("km/s").value
    new_wave = wave * Doppler_factor
    interpolationfunction = interp1d(
        new_wave, model, kind="cubic", fill_value="extrapolate")
    interpolatedModel = interpolationfunction(wave)
    plt.plot(wave, flux)
    plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
    plt.plot(wave, interpolatedModel, color="orange", marker="*")
    plt.show()

    return v_rad_best
   

if __name__ == "__main__":
   #############################################################
   #
   # EXAMPLE for HD 183143
   #
   #############################################################

   # Get one of the spectra -- this is just the 3302 region. 
   pythia = EdiblesOracle()
   List = pythia.getFilteredObsList(object=["HD 183143"], MergedOnly=True, Wave=3302.0)
   test = List.values.tolist()
   filename = test[0]
   print(filename)
   wrange = [3301.5, 3304.0]
   sp = EdiblesSpectrum(filename)
   wave = sp.wave
   flux = sp.flux
   idx = np.where((wave > wrange[0]) & (wave < wrange[1]))
   wave = wave[idx]
   flux = flux[idx]
   flux = flux / np.median(flux)
   
   # A key parameter is the velocity resolution. We will estimate that here from the data, 
   # Assuming an oversampling rate of 2. 
   v_resolution = 2. * (wave[1]-wave[0])/wave[0] * cst.c.to("km/s").value
   print("Velocity resolution:", v_resolution)
   
   # Now create a model for a single cloud, at v_rad = 0. Let's pick N and b such that they produce a decent absorption line. 
   lambda0 = [3302.369, 3302.978]
   f = [8.26e-03, 4.06e-03]
   gamma = [6.280e7, 6.280e7]
   b = [1.0]
   N = [5e13]
   v_rad = [0.0]

   ReferenceAbsorptionModel = voigt_absorption_line(
    wave,
    lambda0=lambda0,
    b=b,
    N=N,
    f=f,
    gamma=gamma,
    v_rad=v_rad,
    v_resolution=v_resolution,
    )
   #plt.plot(wave, flux)
   #plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
   #plt.plot(wave, AbsorptionLine, color="orange", marker="*")
   #plt.show()

   # Now get the guess for the v_rad of the first component.    
   v_rad_guess = determine_vrad_from_correlation(wave, flux, ReferenceAbsorptionModel)

   # Now let's look at step two. Make a model that more or less fits the strongest line, 
   # then calculate residuals and search again. 
   # These parameters produce an OK fit. 
   lambda0 = [3302.369, 3302.978]
   f = [9.21e-03, 4.60e-03]
   gamma = [6.280e7, 6.280e7]
   b = [2.0]
   N = [1e14]
   v_rad = [v_rad_guess]
   v_resolution = 4.00

   AbsorptionLine = voigt_absorption_line(
    wave,
    lambda0=lambda0,
    b=b,
    N=N,
    f=f,
    gamma=gamma,
    v_rad=v_rad,
    v_resolution=v_resolution,
    )
   # Calculate residuals
   ResidualFlux = flux-AbsorptionLine + 1

   plt.plot(wave, flux)
   plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
   plt.plot(wave, AbsorptionLine, color="orange", marker="*")
   plt.plot(wave,ResidualFlux, color="green")
   plt.show()

   # And run the radial velocity correlation thing again.... 
   v_rad_guess2 = determine_vrad_from_correlation(wave, ResidualFlux, ReferenceAbsorptionModel)

   # Now create a 2-line model   
   lambda0 = [3302.369, 3302.369, 3302.978, 3302.978]
   f = [9.21e-03, 9.21e-03, 4.60e-03, 4.60e-03]
   gamma = [6.280e7, 6.280e7, 6.280e7, 6.280e7]
   b = [2.0, 2.0, 2.0, 2.0]
   N = [1e14, 8e13, 1e14, 8e13]
   v_rad = [v_rad_guess, v_rad_guess2, v_rad_guess, v_rad_guess2]
   v_resolution = 4.00
   AbsorptionLine = voigt_absorption_line(
    wave,
    lambda0=lambda0,
    b=b,
    N=N,
    f=f,
    gamma=gamma,
    v_rad=v_rad,
    v_resolution=v_resolution, 
    )
   ResidualFlux = flux-AbsorptionLine + 1
   plt.plot(wave, flux)
   plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
   plt.plot(wave, AbsorptionLine, color="orange", marker="*")
   plt.plot(wave,ResidualFlux, color="green")
   plt.show()