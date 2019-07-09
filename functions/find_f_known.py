import numpy as np

class atomic_lines:
    
       
    def __init__(self, filename):
        self.filename = filename
        with open(self.filename) as f:

            first_line = f.readline()
            names = first_line.split('|')
            indeces = []

            for i in range(len(names)):
                indeces.append(len(names[i]))

            for i in range(1, len(indeces)):
                indeces[i] += indeces[i-1] +1


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
            data = [self.wavelength, self.species, self.TT, self.Term, self.J_ik, self.f_ik, self.TPF, self.LVL_EN_CM_1, self.REF, self.new_line]
            
            for line in f:

                start = 0
                for i in range(len(indeces)):
                    stop = indeces[i] 
                    data[i].append(line[start:stop].strip())
                    start = indeces[i]

            data.insert(0, names)

            self.data = np.asarray(data)
   
    def find_index(self, ion, wave):

        indeces = []
        species = self.species
        for index in range(len(species)):

            if species[index] == ion:
                indeces.append(index)

        diff = []
        wavelength = self.wavelength
        for i in indeces:
            wave_table = float(wavelength[i])
            diff.append(np.abs(wave_table-wave))

        index = np.argmin(diff)
        return index
    
    def get_f_known(self, ion, wave):
        
        index = self.find_index(ion, wave)
        f_known = float(self.f_ik[index])
        return f_known

    def get_lvl_en_cm_1(self, ion, wave):

        index = self.find_index(ion, wave)
        lvl_en_cm_1 = self.LVL_EN_CM_1[index]
        return lvl_en_cm_1


if __name__ == "__main__":

    obj = atomic_lines('/home/ranjan/python/edibles/atomic_lines.txt')

    ion = 'Na I'
    wave = 3302.7

    print(obj.get_f_known(ion, wave))
    print(obj.get_lvl_en_cm_1(ion, wave))
        
#OLD CODE        
# import numpy as np

# def find_F(ion, wave):


#     filename = '/home/ranjan/python/edibles/atomic_lines.txt'

#     with open(filename) as f:

#         first_line = f.readline()
#         names = first_line.split('|')
#         indeces = []

#         for i in range(len(names)):
#             indeces.append(len(names[i]))
       
#         for i in range(1, len(indeces)):
#             indeces[i] += indeces[i-1] +1


#         wavelength = []
#         species = []
#         TT = []
#         Term = []
#         J_ik = []
#         f_ik = []
#         TPF = []
#         LVL_EN_CM_1 = []
#         REF = []
#         new_line = []
#         data = [wavelength, species, TT, Term, J_ik, f_ik, TPF, LVL_EN_CM_1, REF, new_line]
#         print(data)
#         for line in f:
            
#             start = 0
#             for i in range(len(indeces)):
#                 stop = indeces[i] 
#                 data[i].append(line[start:stop].strip())
#                 start = indeces[i]
        
#         data.insert(0, names)
        
#         data = np.asarray(data)


#     indeces = []
#     for index in range(len(species)):

#         if species[index] == ion:
#             indeces.append(index)

#     diff = []
#     for i in indeces:
#         wave_table = float(wavelength[i])
#         diff.append(np.abs(wave_table-wave))

#     index = np.argmin(diff)
#     f_known = float(f_ik[index])
#     return f_known


# if __name__ == '__main__':

#     ion = 'Na I'
#     wave = 3302.7
#     f = find_F(ion, wave)
#     print(f)