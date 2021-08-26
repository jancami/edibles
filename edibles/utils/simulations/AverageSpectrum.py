"""This script computes the average spectrum for a given DIB and sightline.

For a given sightline, a profile is created with the average spectrum. The S/N ratio
is computed too. The resulting profile is for a given DIB, in a range of (DIB-4, DIB+4).
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from edibles.utils.simulations.SRC.Functions import Signal_Noise_Calculator, weighted_average
from edibles.utils.edibles_oracle import EdiblesOracle
from edibles.utils.edibles_spectrum import EdiblesSpectrum


def CreateAverageSpectrum(DIB, Target, save_to_file=False, save_figure=False, verbose=True):
    """Create a DIB profile with an average spectrum for a sightline.

    Args:
        DIB (int): Center wavelength of DIB.
        Target (str): Sightline of observation. Usually "HD {integer number}"
        save_to_file (bool, optional): If True, it saves the resulting profile.
            Defaults to False.
        save_figure (bool, optional): If True, it creates a plot of the resulting profile
            and saves it. This plot contains the SN ratio. Defaults to False.
        verbose (bool, optional): If True, print the target date and resulting profile.
            Defaults to True.

    Returns:
        2darray: Resulting average profile in the form of (wavelength, intensity).
    """
    # Get Edibles files of that sightline.
    oracle = EdiblesOracle()
    List = oracle.getFilteredObsList([Target], Wave=DIB, MergedOnly=True)

    # Dataframe to save the data.
    df = pd.DataFrame()

    # Offset for plotting
    offset = 0.05

    # Check if data is found.
    if len(List) == 0:
        print("Sightline files not found!")
        return np.array([[0], [0]])

    # Iterate over found datafiles.
    for file in List:

        # Get Edibles data.
        sp = EdiblesSpectrum(file)

        # Get target observation date and print it.
        target_date = str(sp.datetime.date()).replace('-', '_')
        if verbose:
            print(target_date)

        # Obtain data in wavelength range of interest.
        sp.getSpectrum(xmin=DIB-4, xmax=DIB+4)

        # Get wavelenth and flux
        DIB_wavelength = np.asarray(sp.grid.byteswap().newbyteorder(), dtype='float64')
        DIB_flux = np.asarray(sp.interp_bary_flux.byteswap().newbyteorder(), dtype='float64')
        DIB_flux = DIB_flux/np.max(DIB_flux)

        # Select datapoints to compute SN ratio.
        cont_x1 = np.array(DIB_wavelength[-20:])
        cont_y1 = np.array(DIB_flux[-20:])
        cont_x2 = np.array(DIB_wavelength[:15])
        cont_y2 = np.array(DIB_flux[:15])

        # Get SN ratio
        SN1, Fit1 = Signal_Noise_Calculator(cont_x1, cont_y1)
        SN2, Fit2 = Signal_Noise_Calculator(cont_x2, cont_y2)

        # Save data unvertainty
        Signal_Noise = 0.5*(SN1+SN2)
        New_Noise = DIB_flux/Signal_Noise
        uncertainty = np.full(DIB_wavelength.shape, np.mean(New_Noise), dtype=float)

        # Save results to dataframe
        df[target_date+"_data"] = DIB_flux
        df[target_date+"_error"] = uncertainty

        # Add observation to plot.
        if save_figure:
            plt.plot(DIB_wavelength, DIB_flux+offset, label=target_date)
            offset = offset+0.05

    # List to save stacking results.
    weig_avg = []
    weig_avg_err = []

    # Print all the data obtained.
    if verbose:
        print(df)

    # Compute weighted average.
    for i in range(len(df.index)):
        values = df[df.columns[::2]].iloc[i].values
        error = df[df.columns[1::2]].iloc[i].values
        avg_weig = weighted_average(values, error)
        weig_avg.append(avg_weig[0])
        weig_avg_err.append(avg_weig[1])

    # Normalization
    weig_avg = weig_avg/np.max(weig_avg)

    # Save to dataframe
    df["Weighted_Average"] = weig_avg
    df["Weighted_Average_Error"] = weig_avg_err

    # Add resulting profile to plot.
    if save_figure:
        plt.plot(DIB_wavelength, weig_avg, label='weighted_average')
        plt.legend()
        plt.savefig("AverageSpectra_"+str(DIB)+".pdf")
        plt.close()

    # Save restul to file.
    if save_to_file:
        final = pd.DataFrame({"Wavelength": DIB_wavelength,
                             "Flux": df["Weighted_Average"].to_numpy()})
        final.to_csv("Data/AverageSpectraData/"+str(DIB)+"/"+Target+"_avg_spectra.csv")

    # Return resulting profile.
    return(DIB_wavelength, df["Weighted_Average"].to_numpy())


if __name__ == "__main__":

    # Example
    CreateAverageSpectrum(6614, 'HD 170740')
