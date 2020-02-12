import os
import numpy as np
import csv
from edibles.edibles import DATADIR, EDIBLES_PYTHONDIR, DATARELEASE
from edibles.edibles.functions.edibles_spectrum import EdiblesSpectrum

def createObsList():
    """Set up of the list of FITS files in the data directory with the necessary
    information that we need to supply to the oracle to work. Essentially, we
    will recursively list all of the FITS files first, then open each of them to
    read the header and extract the information. Finally, store everything in a
    text file."""




    print(DATADIR)
    allfitsfiles = []
    for path, dirs, files in os.walk(DATADIR):
        for file in files:
            if file.endswith(".fits"):
                fullfilename = os.path.join(path, file)
                relative_filename = fullfilename[len(DATADIR):]
                print(relative_filename)
                allfitsfiles.append(relative_filename)

    print(len(allfitsfiles))


    class obsfile:
        def __init__(self):
            self.filename = ''
            self.object = ''
            self.date_obs = ''
            self.setting = ''
            self.wave_min = ''
            self.wave_max = ''
            self.ra = ''
            self.dec = ''
            self.exptime = ''
            # self.order=0
            # self.merged=False


    n_files = len(allfitsfiles)
    full_list = [obsfile() for i in range(n_files)]

    for count in range(len(allfitsfiles)):
        full_list[count].filename = allfitsfiles[count]
        # print(datadir + full_list[count].filename)
        spec = EdiblesSpectrum(full_list[count].filename)
        # print(spec.header)
        full_list[count].object = spec.header["OBJECT"]
        full_list[count].date_obs = spec.header["DATE-OBS"]
        full_list[count].ra = spec.header["RA"]
        full_list[count].dec = spec.header["DEC"]
        full_list[count].exptime = spec.header["EXPTIME"]
        print(allfitsfiles[count])
        if "HIERARCH ESO INS GRAT1 WLEN" in spec.header:
            full_list[count].setting = int(spec.header["HIERARCH ESO INS GRAT1 WLEN"])
        if "HIERARCH ESO INS GRAT2 WLEN" in spec.header:
            full_list[count].setting = int(spec.header["HIERARCH ESO INS GRAT2 WLEN"])
        wave, flux = spec.getSpectrum()
        full_list[count].wave_min = "{:.1f}".format(np.min(wave))
        full_list[count].wave_max = "{:.1f}".format(np.max(wave))
        del spec

    # Create arrays of formatted strings to print to a csv file now.
    pstrings = [['Object', 'RA', 'DEC', 'DateObs', 'setting', 'WaveMin', 'WaveMax',
                'Filename']]
    for count in range(n_files):
        pstrings.append([full_list[count].object, full_list[count].ra,
                         full_list[count].dec, full_list[count].date_obs,
                         full_list[count].setting, full_list[count].wave_min,
                         full_list[count].wave_max, full_list[count].filename])

    # Time to print things out! Let's use csv format to do that.
    outfile = EDIBLES_PYTHONDIR + "/data/" + DATARELEASE + "_ObsLog.csv"
    length_checker = np.vectorize(len)
    all_lengths = length_checker(allfitsfiles)
    print(np.max(all_lengths))
    with open(outfile, 'w') as csvFile:
        writer = csv.writer(csvFile)
        writer.writerows(pstrings)

    return



if __name__ == "__main__":
    createObsList()