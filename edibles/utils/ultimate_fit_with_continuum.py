#################################################################
# fitting with continuum 
#################################################################


from edibles.utils.voigt_profile_new_A import *

# useful packages to import
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import glob as g
from scipy.ndimage import gaussian_filter
from lmfit import Parameters, minimize,Model
from scipy.optimize import fmin

from scipy.interpolate import CubicSpline
from scipy.optimize import curve_fit
from scipy.special import wofz



# Wrapper function around fit_multi_voigt_absorptionlines. 
def multi_voigt_absorption_line( **params_list):
    """
    This function is essentially a wrapper around voigt_absorption_line in voigt_profile.py, that exists
    to allow it to create an instance of the Model class and create models with multiple Voigt commponents
    (both transitions and clouds). 
    This function will parse the list of parameters, and reformat them to call the voigt_absorption_line function. 
    It will then call that function, and return the resulting model. 

    Args:
        wavegrid: the wavelength grid on which to calculate the models. 
        n_trans: the number of unique transitions to calculate. 
                 If n_trans == 1, just one rest wavelength, f and gamma value is passed on. 
                 If n_trans == 2, we have a doublet and so on. 
        lambda0, lambda1, ...: rest wavelengths for each transition. 
        f0, f1, ...:           oscillator strengths for each transition. 
        gamma0, gamma1, ...:   Lorentz broadening parameter for each transition. 
        n_components: the number of different cloud components for which to calculate the profiles. 
        b0, b1, ....:  Doppler b-values for each cloud component. 
        N0, N1, ....:  column densities for each cloud component. 
        v_rad0, v_rad1, ...:  radial velocities for each cloud component/transition. 
        v_resolution: the velocity resolution (in km/s) of the desired final result. 
        n_step: the number of steps to sample the Voigt profile (default: 25). 
        debug:   Boolean: print debug info while running or not. 
    Returns:
    np.array
        Model array with normalized Voigt profiles corresponding to the specified parameters. 
    """

    # adding different parameters in the lmfit param_list 
    n_trans = params_list['n_trans']
    wavegrid = params_list['wavegrid']
    n_components = params_list['n_components'] # cloud components
    v_resolution = params_list['v_resolution']
    n_step = params_list['n_step']
    n_species = params_list['n_species']
    
    
    
    #----------------------------chemical properties of species----------------------------------------------------------------
    # Initialize the lambda,f,gamma arrays for voigt_absorption_line
    all_lambda = np.empty(n_trans) 
    all_f = np.empty(n_trans)
    all_gamma = np.empty(n_trans)
    map = np.empty(n_trans)

    for i in range(n_trans):
        all_lambda[i] = params_list[f'lambda{i}']
        all_f[i] = params_list[f'f{i}']
        all_gamma[i] = params_list[f'gamma{i}']
        map[i] = params_list[f'map{i}'] # mapping the species with lambdas
    #==============================================================================================================================
   

    
     # ---------------cloud component(right now we only assume single cloud)--------------------------------------------------------
    all_v_rad = np.empty(n_components)
    for i in range(n_components): # for single cloud the loop will run only once

        all_v_rad[i] = params_list[f'v_rad{i}']
    #==============================================================================================================================
    
    
   # ---------------species--------------------------------------------------------------------------------------------------
    all_b= np.empty(n_species)
    all_N= np.empty(n_species)
    
    for i in range(n_species):
        all_b[i] = params_list[f'b{i}']
        all_N[i] = params_list[f'N{i}'] 
    #==============================================================================================================================
   
   
   
   #**********************  continuum  **********************************************************************************************
    anchor_fluxes = np.interp(anchor_wavelengths, wavegrid, flux) # anchor flux corresponding to anchor wavelengths
    spline_fit = CubicSpline(anchor_wavelengths, anchor_fluxes) #using cubic spline to fit the continuum
    continuum_model = spline_fit(wavegrid) # continuum model
    #**********************************************************************************************************************
    
    
    # Normalize the spectrum by dividing by the continuum model
    normalized_flux = flux / continuum_model
     
     
     

    # Now call voigt_absorption_line with these parameters.... 
    model = voigt_absorption_line_new(wavegrid, lambda0=all_lambda, f=all_f, gamma=all_gamma, b=all_b, N=all_N, v_rad=all_v_rad, 
                                  v_resolution=v_resolution, n_step=n_step, map=map, debug=False)
    
    
    
    # adding continuum to the voigt model 
    combined_fit = model * continuum_model
    
    return combined_fit






#The main fitting function 
def fit_multi_voigt_absorptionlines(wavegrid=np.array, ydata=np.array, restwave=np.array, f=np.array, gamma=np.array, 
             b=np.array, N=np.array, v_rad=np.array, v_resolution=0., n_step=0, std_dev = 1,map =np.array,anchor_wavelengths=np.array([])):
    """
    This function will take an observed spectrum contained in (wavegrid, ydata) and fit a set of Voigt profiles to
    it. The transitions to consider are specified by restwave, f, and gamma, and can be single floats or numpy arrays 
    containing the information for all transitions that we want to consider. 
    There can then be 1 or more clound components to consider in the sightline. The components are specified by b, N and v_rad
    that can be again floats or numpy arrays. 

    This function essentially creates a Model instance from the multi_voigt_absorption_line function, and most of the work done
    here is to "translate" the parameters from arrays to unique parameters (using the Parameters class) that will then be parsed 
    by multi_voigt_absorption_line.
    
     Args:
        wavegrid (float64): Wavelength grid (in Angstrom) on which the final result is desired.
        lambda0 (float64): Central (rest) wavelength(s) for the absorption lines, in Angstrom (numpy array, p elements)
        f (float64): The oscillator strength (dimensionless) for the absorption lines (numpy array, p elements)
        gamma (float64): Lorentzian gamma (=HWFM) component for the absorption lines (numpy array, p elements)
        b (float64): The b parameter (Gaussian width), in km/s (numpy array, n elements)
        N (float64): The column density (in cm^{-2}) (numpy array, n elements)
        v_rad (float64): Radial velocity of cloud (in km/s) -- single value!! 
        v_resolution (float64): Instrument resolution in velocity space (in km/s)
        composition (string): n element string array, listing the n species corresponding to N (n elements, optional)
        map (int): p-element numpy array that for each line points to what species it corresponds to (index into composition)
        n_step (int): no. of point per FWHM length, governing sampling rate and efficiency
        debug (bool): If True, info on the calculation will be displayed
        
    Returns:
            Fitted plot for single/multiple lines with single cloud model.

    """
    
    # We should probably do lots of parameter checking first!!! To be done later.... 


#------------------------------------------------------------------------------  
    # How many transitions do we have? 
    restwave = np.asarray(restwave)  # reatwave is basically the lambdas
    n_trans = restwave.size # number of lines
    
    #sorting different parameters to use them later
    N = np.asarray(N)
    n_species = N.size # number of species
    map = np.asarray(map) 

    b = np.asarray(b)
    v_rad = np.asarray(v_rad)
    n_components = v_rad.size # cloud component
#==============================================================================




#----------------Now create the Parameters class----------------------------------------------
 
    params = Parameters()
    # Create the parameters for the transitions. Those should *not* be free parameters for lmfit. 
    if n_trans == 1:
        params.add('lambda0', value=restwave, vary=False)
        params.add('f0', value=f, vary=False)
        params.add('gamma0', value=gamma, vary=False)
        params.add('map0', value=map[0], vary=False)
    else: # for multiple lines 
        
        for i in range(n_trans):
             
            params.add(f'lambda{i}', value=restwave[i], vary=False)
            params.add(f'f{i}', value=f[i], vary=False) 
            params.add(f'gamma{i}', value=gamma[i], vary=False) 
            params.add(f'map{i}', value=map[i], vary=False)

    # Also create parameters for the other keywords -- these too should *not* be free parameters. 
    params.add('v_resolution', value=v_resolution, vary=False)
    params.add('n_step', value=n_step, vary=False)
    params.add('n_trans', value=n_trans, vary=False)
    params.add('n_components', value=n_components, vary=False)
    params.add('n_species', value=n_species, vary=False)
    
    
   
    # Create the parameters for the clouds.  These *should* be free parameters.
    if n_components == 1:
        params.add('v_rad0', value=v_rad)
        
    else:
        for i in range(n_components):

            params.add(f'v_rad{i}', value=v_rad[i])
            
    for i in range(n_species): # for Li the for loop will run twice
        params.add(f'b{i}', value=b[i],min=0) #b =[b0,b1]
        params.add(f'N{i}', value=N[i],min=0) # N =[N0,N1]
        
#=============================================================================================================  
       

    # Now create the Model instance. 
    voigtmod = Model(multi_voigt_absorption_line, independent_vars=['wavegrid', 'flux', 'anchor_wavelengths'])
    #print("Resolution: ", v_resolution)

    # and do the fitting with the parameters we have created. 
    result = voigtmod.fit(ydata, params, wavegrid=wavegrid, flux=ydata, anchor_wavelengths=anchor_wavelengths, weights=1/std_dev)
    return result
















'''
EXAMPLES
'''


# ========================================================================
# HERE ARE FEW EXAMPLES THAT CAN BE RUN TO TEST THE CODE
# ========================================================================
if __name__ == "__main__":
    
    from math import *

    show_example = 1  # <---------------- GIVE YOUR INPUT HERE ************************




    if show_example == 1:
        #############################################################
        #
        # EXAMPLE 1: Li (Li7 and Li6 ), with single cloud
        #
        #############################################################
        pythia = EdiblesOracle()
        List = pythia.getFilteredObsList(object=["HD 147084"], MergedOnly=False, Wave=6708.0)
        test = List.values.tolist()
        print(test)
        for i in test:
            print(i)
        filename = test[2]
        print("the file is ",filename)



        wrange = [6706.5, 6708.9]
        sp = EdiblesSpectrum(filename)
        sp.getSpectrum(xmin=wrange[0], xmax=wrange[1])
        wave = sp.bary_wave
        flux = sp.bary_flux

        # Filter the spectrum within the specified wavelength range
        idx = np.where((wave > wrange[0]) & (wave < wrange[1]))
        wave = wave[idx]
        flux = flux[idx]

        # Normalize the flux
        flux = flux / np.median(flux)

        restwave = [6707.761, 6707.912, 6707.921, 6708.072]
        f = [0.4982, 0.2491, 0.4982/18, 0.2491/18]
        gamma = [3.69e7, 3.69e7, 3.69e7, 3.69e7]

        #continuum fit
        anchor_wavelengths= [6706.95,6707.2,6707.3,6708.2,6708.8]
        anchor_fluxes = np.interp(anchor_wavelengths, wave, flux)


        fitresult = fit_multi_voigt_absorptionlines(wavegrid=wave, ydata= flux,f= [0.4982,0.2491,0.4982,0.2491], restwave=restwave, gamma=gamma, 
                            b=[2,2], N=[0.00073e13,0.00073e13/18], v_rad=-7.1, v_resolution=3, n_step=25, map =[0,0,1,1],std_dev= 0.0035554091884622933,anchor_wavelengths=anchor_wavelengths)
        fitresult.params.pretty_print()

        print("chi-square value ",fitresult.chisqr)
        print("reduced chi-square value ",fitresult.redchi)

        plt.scatter(anchor_wavelengths, anchor_fluxes, marker='D', c='darkorange', label='Anchor Points')
        plt.plot(wave, flux,color ='gray',label ='data')
        plt.plot(wave,fitresult.best_fit,color ='purple',label ="fit")
        plt.xlabel("Wavelength ($\AA$)")
        plt.ylabel("Normalised flux")
        plt.grid()
        plt.legend()
        plt.show()

    elif show_example == 2:
    
        #############################################################
        #
        # EXAMPLE 2: Na doublet: multiple lines with continuum fit , single cloud.
        #
        #############################################################


        lambda0 = [3302.369, 3302.978]
        f = [8.26e-03, 4.06e-03]
        gamma = [6.280e7, 6.280e7]
        b = [0.8]
        N = [1.8e13]
        map = [0, 0]
        v_rad = -10
        v_resolution = 5.75

        pythia = EdiblesOracle()
        List = pythia.getFilteredObsList(
            object=["HD 145502"], MergedOnly=False, Wave=3302.0
        )
        test = List.values.tolist()
        filename = test[0]
        print(filename)
        wrange = [3301.5, 3304.0]
        sp = EdiblesSpectrum(filename)
        sp.getSpectrum(wrange[0],wrange[1])
        wave = sp.bary_wave
        flux = sp.bary_flux
        idx = np.where((wave > wrange[0]) & (wave < wrange[1]))
        wave = wave[idx]
        flux = flux[idx]
        flux = flux / np.median(flux)

        # Filter the spectrum within the specified wavelength range
        idx = np.where((wave > wrange[0]) & (wave < wrange[1]))
        wave = wave[idx]
        flux = flux[idx]

        # Normalize the flux
        flux = flux / np.median(flux)

        #continuum fit
        anchor_wavelengths= [6706.95,6707.2,6707.3,6708.2,6708.8]
        anchor_fluxes = np.interp(anchor_wavelengths, wave, flux)



        fitresult = fit_multi_voigt_absorptionlines(wavegrid=wave, ydata= flux,f= f, restwave=lambda0, gamma=gamma, 
                               b=b, N=N, v_rad=v_rad, v_resolution=v_resolution, n_step=25, map =[0,0],std_dev= 0.0035554091884622933,anchor_wavelengths=anchor_wavelengths)
        fitresult.params.pretty_print()

        print("chi-square value ",fitresult.chisqr)
        print("reduced chi-square value ",fitresult.redchi)

        plt.plot(wave, flux,color ='gray',label ='data')
        plt.plot(wave,fitresult.best_fit,color ='purple',label ="fit")
        plt.xlabel("Wavelength ($\AA$)")
        plt.ylabel("Normalised flux")
        plt.grid()
        plt.legend()
        plt.show()

    elif show_example == 3:
        #############################################################
        #
        # EXAMPLE 3: 12CH+ and 13CH+ lines with continuum fit : multiple lines, single cloud.
        #
        #############################################################
        lambda0 = [4232.288, 4232.548]
        f = [0.005450, 0.005450]
        gamma = [1e8,1e8]
        b = [1, 1.5]
        N = [4.7e13/50, 4.7e13]
        map = [0, 1]
        v_rad = -15
        v_resolution = 3

        pythia = EdiblesOracle()
        List = pythia.getFilteredObsList(
            object=["HD 149757"], MergedOnly=False, Wave=4232.6
        )
        print(List)
        test = List.values.tolist()
        print(test)
        filename = test[1]
        print(filename)
        wrange = [4231, 4233]
        sp = EdiblesSpectrum(filename)
        sp.getSpectrum(wrange[0],wrange[1])
        wave = sp.bary_wave
        flux = sp.bary_flux
        idx = np.where((wave > wrange[0]) & (wave < wrange[1]))
        wave = wave[idx]
        flux = flux[idx]
        flux = flux / np.median(flux)

        # Filter the spectrum within the specified wavelength range
        idx = np.where((wave > wrange[0]) & (wave < wrange[1]))
        wave = wave[idx]
        flux = flux[idx]

        # Normalize the flux
        flux = flux / np.median(flux)



        #continuum fit
        anchor_wavelengths= [4231.2,4231.5,4231.75,4231.8,4231.95,4232.6,4232.7,4233]
        anchor_fluxes = np.interp(anchor_wavelengths, wave, flux)



        fitresult = fit_multi_voigt_absorptionlines(wavegrid=wave, ydata= flux,f= f, restwave=lambda0, gamma=gamma, 
                               b=b, N=N, v_rad=v_rad, v_resolution=v_resolution, n_step=25, map =[0,1],std_dev= 0.0015,anchor_wavelengths=anchor_wavelengths)
        fitresult.params.pretty_print()

        print("chi-square value ",fitresult.chisqr)
        print("reduced chi-square value ",fitresult.redchi)

        plt.scatter(anchor_wavelengths, anchor_fluxes, marker='D', c='darkorange', label='Anchor Points')
        plt.plot(wave, flux,color ='gray',label ='data')
        plt.plot(wave,fitresult.best_fit,color ='purple',label ="fit")
        plt.xlabel("Wavelength ($\AA$)")
        plt.ylabel("Normalised flux")
        plt.grid()
        plt.legend()
        plt.show()

    elif show_example == 4:
        #############################################################
        #
        # EXAMPLE 4: CN with continuum fit, single cloud 
        #
        #############################################################

        lambda0 = [ 3874.6024,3874.6059, 3874.783] 
        f = [0.0342,0.0392,0.0114]
        gamma = [1.610e7,1.610e7,1.610e7] 
        b = [2.5,2.5]
        N = [0.38e13,0.38e13/10]
        map = [0,0, 1]
        v_rad = -8
        v_resolution = 3
        std_dev = 0.0037752588672860937 

        pythia = EdiblesOracle()
        List = pythia.getFilteredObsList(
            object=["HD 147889"], MergedOnly=False, Wave=3874.65
        )
        print(List)
        test = List.values.tolist()
        print(test)
        filename = test[-1]
        print(filename)
        wrange = [3874.1, 3875]
        sp = EdiblesSpectrum(filename)
        sp.getSpectrum(wrange[0],wrange[1])
        wave = sp.bary_wave
        flux = sp.bary_flux
        idx = np.where((wave > wrange[0]) & (wave < wrange[1]))
        wave = wave[idx]
        flux = flux[idx]
        flux = flux / np.median(flux)

        # Filter the spectrum within the specified wavelength range
        idx = np.where((wave > wrange[0]) & (wave < wrange[1]))
        wave = wave[idx]
        flux = flux[idx]

        # Normalize the flux
        flux = flux / np.median(flux)

        #continuum fit
        anchor_wavelengths= [3874.15,3874.8,3875]
        anchor_fluxes = np.interp(anchor_wavelengths, wave, flux)


        fitresult = fit_multi_voigt_absorptionlines(wavegrid=wave, ydata= flux,f= f, restwave=lambda0, gamma=gamma, 
                               b=b, N=N, v_rad=v_rad, v_resolution=v_resolution, n_step=25, map =map ,std_dev= std_dev,anchor_wavelengths=anchor_wavelengths)
        fitresult.params.pretty_print()

        print("chi-square value ",fitresult.chisqr)
        print("reduced chi-square value ",fitresult.redchi)


        plt.scatter(anchor_wavelengths, anchor_fluxes, marker='D', c='darkorange', label='Anchor Points')
        plt.plot(wave, flux,label ='data')
        plt.plot(wave,fitresult.best_fit,color ='purple',label ="fit")
        plt.xlabel("Wavelength ($\AA$)")
        plt.ylabel("Normalised flux")
        plt.grid()
        plt.legend()
        plt.show()


