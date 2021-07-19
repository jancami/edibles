import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from edibles.utils.simulations.SRC.Functions import Signal_Noise_Calculator, weighted_average
from edibles.utils.edibles_oracle import EdiblesOracle
from edibles.utils.edibles_spectrum import EdiblesSpectrum


def CreateAverageSpectrum(DIB, Target, Sightline, save_figure=False,
                          save_to_file=False, show_figure=False):

    oracle = EdiblesOracle()
    List = oracle.getFilteredObsList([Target], Wave=DIB, MergedOnly=True)

    df = pd.DataFrame()
    i = 0.05
    for file in List:
        sp = EdiblesSpectrum(file)
        target_date = str(sp.datetime.date()).replace('-', '_')
        print(target_date)
        sp.getSpectrum(xmin=DIB-4, xmax=DIB+4)

        DIB_wavelength = np.asarray(sp.grid.byteswap().newbyteorder(), dtype='float64')

        DIB_flux = np.asarray(sp.interp_bary_flux.byteswap().newbyteorder(), dtype='float64')
        DIB_flux = DIB_flux/np.min(DIB_flux)
        cont_x1 = np.array(DIB_wavelength[-20:])
        cont_y1 = np.array(DIB_flux[-20:])
        SN1, Fit1 = Signal_Noise_Calculator(cont_x1, cont_y1)

        cont_x2 = np.array(DIB_wavelength[:15])
        cont_y2 = np.array(DIB_flux[:15])
        SN2, Fit2 = Signal_Noise_Calculator(cont_x2, cont_y2)

        Signal_Noise = 0.5*(SN1+SN2)
        New_Noise = DIB_flux/Signal_Noise

        uncertainty = np.full(DIB_wavelength.shape, np.mean(New_Noise), dtype=float)
        print(i)

        df[target_date+"_data"] = DIB_flux
        df[target_date+"_error"] = uncertainty
        if show_figure:
            plt.plot(DIB_wavelength, DIB_flux+i, label=target_date)
            i = i+0.05

    weig_avg = []
    weig_avg_err = []
    print(df)
    for i in range(len(df.index)):
        values = df[df.columns[::2]].iloc[i].values
        error = df[df.columns[1::2]].iloc[i].values
        avg_weig = weighted_average(values, error)
        weig_avg.append(avg_weig[0])
        weig_avg_err.append(avg_weig[1])

    df["Weighted_Average"] = weig_avg
    df["Weighted_Average_Error"] = weig_avg_err
    if save_figure:
        plt.plot(DIB_wavelength, weig_avg, label='weighted_average')
        plt.legend()
        plt.savefig("AverageSpectra_"+str(DIB)+".pdf")
        plt.close()

    if save_to_file == True:
        final = pd.DataFrame({"Wavelength": DIB_wavelength.to_numpy(),
                             "Flux": df["Weighted_Average"].to_numpy()})
        final.to_csv("Data/AverageSpectraData/"+str(DIB)+"/"+Sightline+"_avg_spectra.csv")
    return(DIB_wavelength.to_numpy(), df["Weighted_Average"].to_numpy())


if __name__ == "__main__":
    CreateAverageSpectrum(6614, 'HD 170740')
