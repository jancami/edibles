import os
import sys
import numpy as np
import csv
from astropy.io import fits
from edibles import DATADIR, PYTHONDIR, DATARELEASE
from edibles.utils.edibles_spectrum import EdiblesSpectrum


def createObsList(dryrun=True):
    """Set up of the list of FITS files in the data directory with the necessary
    information that we need to supply to the oracle to work. Essentially, we
    will recursively list all of the FITS files first, then open each of them
    to read the header and extract the information. Finally, store everything
    in a text file.

    :param dryrun: If True, dont save list. Default=True
    :type drtyrun: bool


    """

    print("Data Release in " + DATADIR)
    allfitsfiles = []
    for path, dirs, files in os.walk(DATADIR):
        for file in files:
            if file.endswith(".fits"):
                fullfilename = os.path.join(path, file)
                relative_filename = fullfilename[len(DATADIR):]
                # print(relative_filename)
                allfitsfiles.append(relative_filename)

    print(len(allfitsfiles))

    class obsfile:
        def __init__(self):
            self.filename = ""
            self.object = ""
            self.date_obs = ""
            self.setting = ""
            self.order = ""
            self.wave_min = ""
            self.wave_max = ""
            self.ra = ""
            self.dec = ""
            self.exptime = ""
            # self.order=0
            # self.merged=False

    n_files = len(allfitsfiles)
    full_list = [obsfile() for i in range(n_files)]

    for count in range(len(allfitsfiles)):
    #for count in range(1000):
        print(count, n_files, allfitsfiles[count])
        full_list[count].filename = allfitsfiles[count]
        #print(DATADIR + full_list[count].filename)
        #spec = EdiblesSpectrum(full_list[count].filename)
        # EdiblesSpectrum has become too slow -- takes 5 sec per file! 
        # So faster to just work with header! 
        hdu = fits.open(DATADIR + full_list[count].filename)
        header = hdu[0].header
        crval1 = header["CRVAL1"]
        cdelt1 = header["CDELT1"]
        nwave = header["NAXIS1"]

        objectstring = header["OBJECT"]
        # We need to make sure we have a consistent format. 
        # So replace kappa Ori with HD 38771; lambda Sco by HD 158926
        # And lots of inconsistencies with spaces! 
        if objectstring == 'kappa Ori': obectstring = 'HD 38771'
        if objectstring == 'lambda Sco': obectstring = 'HD 158926'
        # Idea: trim all whitespace, then add one after the HD characters. 
        objectstring = "".join(objectstring.split())
        objectstring = objectstring[:2] + ' ' + objectstring[2:]
        #print(objectstring)
        full_list[count].object = objectstring
        full_list[count].date_obs = header["DATE-OBS"]
        full_list[count].ra = header["RA"]
        full_list[count].dec = header["DEC"]
        full_list[count].exptime = header["EXPTIME"]
        idx_O = allfitsfiles[count].find("_O")
        if idx_O != -1:
            # print(idx)
            idx_dot = allfitsfiles[count].find(".")
            full_list[count].order = (allfitsfiles[count])[idx_O + 2: idx_dot]
        else:
            full_list[count].order = "ALL"
        if "HIERARCH ESO INS GRAT1 WLEN" in header:
            full_list[count].setting = int(header["HIERARCH ESO INS GRAT1 WLEN"])
        if "HIERARCH ESO INS GRAT2 WLEN" in header:
            full_list[count].setting = int(header["HIERARCH ESO INS GRAT2 WLEN"])
        full_list[count].wave_min = "{:.1f}".format(crval1)
        full_list[count].wave_max = "{:.1f}".format(crval1 + cdelt1 * nwave)
    

    # Create arrays of formatted strings to print to a csv file now.
    pstrings = [
        [
            "Object",
            "RA",
            "DEC",
            "DateObs",
            "Setting",
            "Order",
            "WaveMin",
            "WaveMax",
            "Filename",
        ]
    ]
    for count in range(n_files):
        pstrings.append(
            [
                full_list[count].object,
                full_list[count].ra,
                full_list[count].dec,
                full_list[count].date_obs,
                full_list[count].setting,
                full_list[count].order,
                full_list[count].wave_min,
                full_list[count].wave_max,
                full_list[count].filename,
            ]
        )

    # Time to print things out! Let's use csv format to do that.
    outfile = PYTHONDIR + "/data/" + DATARELEASE + "_ObsLog.csv"
    # length_checker = np.vectorize(len)
    # all_lengths = length_checker(allfitsfiles)
    # print(np.max(all_lengths))

    if dryrun is False:
        with open(outfile, "w") as csvFile:
            writer = csv.writer(csvFile)
            writer.writerows(pstrings)
            print('Wrote to file!')
    else:
        spamwriter=csv.writer(sys.stdout)
        spamwriter.writerows(pstrings)
        print('Wrote to stdout!')

    return


if __name__ == "__main__":
    createObsList(dryrun=False)
