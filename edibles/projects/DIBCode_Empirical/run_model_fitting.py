from numpy.polynomial import chebyshev
import numpy as np
import pandas as pd
import os
from edibles.projects.DIBCode_Empirical.data_cleaning import data_cleaning
from edibles.projects.DIBCode_Empirical.format_params_grouped import format_params_grouped
from edibles.projects.DIBCode_Empirical.plot_final_with_components import plot_final_with_components
from edibles.projects.DIBCode_Empirical.iterative_fit import iterative_fit

from edibles.projects.DIBCode_Empirical.calculate_chisqr_manual import calculate_chisqr_manual
from edibles.projects.DIBCode_Empirical.remove_components import remove_components
from edibles.projects.DIBCode_Empirical.refit_components import refit_components
def run_model_fitting(minrange, maxrange, folder, csv_name, removed_components, recalibrated_wavelength, coadded_flux, coadd_SNR, figpath):
    '''
    Runs the model fitting for the target star
    Args:
        minrange: Lowest wavelength of the model plot.
        maxrange: Highest wavelength of the model plot
        folder: Location where the csv plot containing parameters is stored.
        csv_name: Location of the csv in the folder.
    Returns:
        Final parameter data in csv format inside the folder.
        - Plot of the components of the model for each target.
        - List of parameter data.

    '''
    #if __name__ == "__main__":
    if 1 == 1:

        x, data = data_cleaning(x=recalibrated_wavelength, data=coadded_flux, minrange=minrange, maxrange=maxrange)
        """                                                                                                                
        x: Array of wavelength values for the DIB of interest.                                                             
        data: Coadded y-values of the data.                                                                                
        Use the data_cleaning function if the data contains more data than the data for the fitting                        
        """

        noise_level = 1 / coadd_SNR
        uncertainties = np.full_like(data, noise_level)

        # Run the iterative fit
        result, best_params, history = iterative_fit(x, data, uncertainties, max_iterations=10, significance_level=0.01,
                                                     criterion="bic", verbose=True,
                                                     plot=True)  # shows plots of each iteration
        """                                                                                                                
        max_iterations=10,         # increase if you expect many features                                                  
        significance_level=0.01,   # controls when to stop adding components                                               
        criterion="bic",           # can use 'aic' or 'bic'                                                                
        plot=True)                  # shows plots of each iteration                                                        

        """

        # Create the folder first to avoid errors
        # folder = "C:/Users/brook/Desktop"
        if not os.path.exists(folder):
            os.makedirs(folder)

        # Drop whatever looks like continuum noise and refit
        best_params = remove_components(best_params, removed_components)
        result, best_params = refit_components(best_params, x, data, uncertainties)

        best_params_list1 = list(best_params.valuesdict())
        best_params_list2 = list(best_params.valuesdict().values())
        df1 = pd.DataFrame(best_params_list1)
        df2 = pd.DataFrame(best_params_list2)
        df = pd.concat([df1, df2], axis=1, ignore_index=True)

        # Save to the specific location
        df.to_csv(os.path.join(folder, csv_name), index=False)
        print("\nFinal Best-Fit Parameters:")
        print(format_params_grouped(best_params))
        plot_final_with_components(x, data, best_params, figpath)
        return
