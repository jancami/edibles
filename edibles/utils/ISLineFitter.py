import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from edibles.utils.edibles_oracle import EdiblesOracle
from edibles.utils.edibles_spectrum import EdiblesSpectrum
from edibles.utils.voigt_profile import voigt_absorption_line
import astropy.constants as cst
from scipy.interpolate import interp1d
from scipy.stats import pearsonr

class ISLineFitter():
    def __init__(self, wave, flux):
        # wave and flux from edibles spectrum or elsewhere
        self.wave = wave
        self.flux = flux

    def fit(self):
        # this will do the fitting
        pass
    
    def determine_vrad_from_correlation(self,wave, flux, model):
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

        return v_rad_best

    
    
if __name__ == "__main__":
    print("Hello Word!")      
    
    # Get one of the spectra -- this is just the 6708 region. 
    pythia = EdiblesOracle()
    List = pythia.getFilteredObsList(object=["HD 147889"], OrdersOnly=True, Wave=6708.0)
    test = List.values.tolist()
    filename = test[1]
    print(filename)
    wrange = [6707, 6709.5]
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
    lambda0 = [6707.761, 6707.912]
    f = [4.99e-01, 2.49e-01]
    gamma = [6e7,6e7]
    b = [1]
    N = [1.75e10]
    v_rad = 0

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
    # plt.plot(wave, flux)
    # plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
    # plt.plot(wave, ReferenceAbsorptionModel, color="orange", marker="*")
    # plt.show()


    # Now get the guess for the v_rad of the first component.    
    fit = ISLineFitter(wave,flux)
    v_rad_guess = fit.determine_vrad_from_correlation(wave, flux, ReferenceAbsorptionModel)
    print("Highest correlation for v_rad = ",v_rad_guess)
