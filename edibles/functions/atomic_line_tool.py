import numpy as np
from edibles.edibles import EDIBLES_PYTHONDIR


class AtomicLines:

    def __init__(self):
        self.filename = EDIBLES_PYTHONDIR + '/edibles/data/atomic_lines.txt'
        with open(self.filename) as f:

            first_line = f.readline()
            names = first_line.split('|')
            indices = []

            for i in range(len(names)):
                indices.append(len(names[i]))

            for i in range(1, len(indices)):
                indices[i] += indices[i - 1] + 1

            self.wavelength = []
            self.species = []
            self.TT = []
            self.Term = []
            self.J_ik = []
            self.f_ik = []
            self.TPF = []
            self.LVL_EN_CM_1 = []
            self.REF = []
            self.new_line = []
            data = [self.wavelength, self.species, self.TT, self.Term, self.J_ik, self.f_ik, self.TPF, self.LVL_EN_CM_1,
                    self.REF, self.new_line]

            for line in f:
                start = 0
                for i in range(len(indices)):
                    stop = indices[i]
                    data[i].append(line[start:stop].strip())
                    start = indices[i]

            data.insert(0, names)
            self.data = np.asarray(data)

    def findIndex(self, ion, wave):

        indices = []
        species = self.species
        for index in range(len(species)):
            if species[index] == ion:
                indices.append(index)

        diff = []
        wavelength = self.wavelength
        for i in indices:
            wave_table = float(wavelength[i])
            diff.append(np.abs(wave_table - wave))

        min_val_index = np.argmin(diff)
        index = indices[min_val_index]

        return index

    def get_f_known(self, ion, wave):

        index = self.findIndex(ion, wave)
        f_known = float(self.f_ik[index])
        return f_known

    def get_lvl_en_cm_1(self, ion, wave):

        index = self.findIndex(ion, wave)
        lvl_en_cm_1 = self.LVL_EN_CM_1[index]
        return lvl_en_cm_1

    def getLabWavelength(self, ion, wave):
        index = self.findIndex(ion, wave)
        lab_wavelength = float(self.wavelength[index])
        return lab_wavelength

    def getAllLabWavelength(self, ion):
        # returns list of all wavelengths of the ion
        return [float(self.data[1][i]) for i in range(len(self.data[1])) if self.data[2][i] == ion]


if __name__ == "__main__":
    obj = AtomicLines()

    ion = 'Na I'
    wave = 5000

    print(obj.get_f_known(ion, wave))
    print(obj.get_lvl_en_cm_1(ion, wave))
    print(obj.getLabWavelength(ion, wave))
    print(obj.getAllLabWavelength(ion))
