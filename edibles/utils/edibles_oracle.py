import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from edibles import DATADIR
from edibles import PYTHONDIR
from edibles.utils.edibles_spectrum import EdiblesSpectrum


class EdiblesOracle:
    """
    This class will pocess the EDIBLES obs log and target info files.
    Users can then query the oracle for observations matching specific criteria.
    """

    def __init__(self):
        print(DATADIR)
        filename = PYTHONDIR + "/edibles/data/DR4_ObsLog.csv"
        self.obslog = pd.read_csv(filename)
        # print(self.obslog.dtypes)
        # total_rows = len(self.obslog.index)
        # print(total_rows)

    def GetObsListByWavelength(self, wave=None, MergedOnly=False, OrdersOnly=False):
        """
        This function filters the list of Observations to return only those
        that include the requested wavelength.
        We will create a set of boolean arrays that we will then combined
        as the filter.

        :param wave: Wavelength that the returned files will include
        :type wave: float
        :param MergedOnly: Only include spectra from merged orders
        :type MergedOnly: bool
        :param OrdersOnly: Only include individual spectrum orders
        :type OrdersOnly: bool

        """

        # Boolean matches for wavelength.
        if wave is None:
            wave = 5000
        bool_wave_matches = (self.obslog.WaveMin < wave) & (self.obslog.WaveMax > wave)

        # Do we have to filter out merged or single-order spectra? Note that if both
        # MergedOnly and OrdersOnly are True, only the Merged spectra will be returned.

        if MergedOnly and OrdersOnly:
            print("ONLY RETURNING MERGED SPECTRA")

        bool_order = self.obslog != "Z"
        if OrdersOnly is True:
            bool_order = self.obslog.Order != "ALL"
        if MergedOnly is True:
            bool_order = self.obslog.Order == "ALL"

        ind = np.where(bool_wave_matches & bool_order)
        # print(ind)
        return self.obslog.iloc[ind].Filename


if __name__ == "__main__":
    # print("Main")
    pythia = EdiblesOracle()
    List = pythia.GetObsListByWavelength(5000, MergedOnly=True)
    # print(List)
    for filename in List:
        sp = EdiblesSpectrum(filename)
        plt.figure()
        plt.title(filename)
        plt.xlabel("Wavelength (" + r"$\AA$" + ")")
        plt.xlim(5000, 5100)
        plt.plot(sp.wave, sp.flux)
        plt.show()
