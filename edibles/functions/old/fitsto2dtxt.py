from __future__ import print_function
from astropy.io import fits
import numpy as np
import os
import csv
import sys


def fitstotxt(target, filepath, writepath, xmin, xmax):
    """
    Reads a FITS file, applies barycentric and cloud corrections
    and writes a subsection within a given range to a .txt file.

    input:
        sightline
        path/to/fits/database/
        path/to/write/folder
        range max
        range min

    output:
        .txt file of spectrum within given range, tab separated, in a given folder.

        example output:

        # HD144470
        # Wavelength(1/cm)      Relative Intensity
        3303.1797537322504  1.0084753
        3303.199755852153   1.0125096
        3303.219757972055   1.012984
        3303.2397600919576  1.0122045

    example:

    fitstotxt('HD144470', '/data/DR3_fits/', '/home/txtfiles/, 6610, 6618)

    """

    c = 2.99 * (10 ** 8)  # m/s

    arms = ["BLUE_346", "BLUE_437", "REDL_564", "REDU_564", "REDL_860", "REDU_860"]

    if xmin <= 3876:
        l1 = 0
    if xmin <= 4990 and xmin >= 3754:
        l1 = 1
    if xmin <= 5667.9 and xmin >= 4616:
        l1 = 2
    if xmin <= 6693.9 and xmin >= 5668:
        l1 = 3
    if xmin <= 8649.9 and xmin >= 6694:
        l1 = 4
    if xmin >= 8650:
        l1 = 5

    if xmax <= 3876:
        l2 = 0
    if xmax <= 4990 and xmax >= 3754:
        l2 = 1
    if xmax <= 5667.9 and xmax >= 4616:
        l2 = 2
    if xmax <= 6693.9 and xmax >= 5668:
        l2 = 3
    if xmax <= 8649.9 and xmax >= 6694:
        l2 = 4
    if xmax >= 8650:
        l2 = 5

    if l1 == l2:
        warm = [arms[l1]]
    if l1 != l2:
        warm = arms[l1:l2]

    for i in range(len(warm)):
        os.chdir(filepath)
        loc = filepath + target + "/" + warm[i] + "/"
        if warm[i] == "REDL_564" or warm[i] == "REDU_564":
            loc = filepath + target + "/RED_564" + "/"
        if warm[i] == "REDL_860" or warm[i] == "REDU_860":
            loc = filepath + target + "/RED_860" + "/"

        if os.path.isdir(loc) is False:
            print("This object have not been observed yet!")
            return ()
        os.chdir(loc)

    # Min and Max values used in order to shift based on interstellar NaI

    NaIxmin = 3301.0
    NaIxmax = 3305.0
    # NaIxmin = 5895
    # NaIxmax = 5897.5

    # ymin = 0.5
    # ymax = 1

    lambda_res = 3302.8
    # lambda_res = 5897.5

    # first_peak_index = []
    # second_peak_index = []
    # third_peak_index = []
    visible_x = []
    x_to_plot = []
    visible_y = []
    y_to_plot = []
    points = []
    DIB_x = []
    DIB_y = []

    # file_list = [os.path.basename(q) for q in glob.glob(path + '*.fits')]
    path = loc
    file_list = os.listdir(path)

    for file_number in range(len(file_list)):

        file_name = file_list[file_number]

        if file_name.split(".")[-2] == "fits":

            hdulist = fits.open(path + file_name)
            hdu = hdulist[0]
            data_towork = hdu.data
            first_val = hdu.header["CRVAL1"]
            stepsize = hdu.header["CDELT1"]
            final_val = first_val + (stepsize * len(data_towork))
            x = np.linspace(first_val, final_val, len(data_towork))
            bcf = hdu.header["HIERARCH ESO QC VRAD BARYCOR"]

            # Analyze the NaI lines to get the proper shift
            for i in range(0, len(data_towork)):
                if NaIxmin <= x[i] <= NaIxmax:
                    visible_x.append(x[i])
                    visible_y.append(data_towork[i])
                if xmin <= x[i] <= xmax:
                    DIB_x.append(x[i])
                    DIB_y.append(data_towork[i])
                else:
                    continue

            if len(visible_x) == 0:
                continue
            else:
                x_to_plot.extend(visible_x)
                y_to_plot.extend(visible_y)

    for j in range(0, len(x_to_plot)):
        point = (((x_to_plot[j])), (y_to_plot[j]))
        points.append(point)
    points.sort()

    bary_corr = 1 + (bcf / c)

    xpoint = [x[0] for x in points]
    ypoint = [y[1] for y in points]
    lowest_min = np.argmin(ypoint / max(ypoint))
    lambda_obs = xpoint[lowest_min]
    cloud_vel = c * (lambda_obs - lambda_res) / lambda_res
    cloud_vel_corr = 1 - (cloud_vel / c)

    # correction for wavelenth taking into account barycentric velocity and cloud velocity.
    total_wavelength_correction = bary_corr * cloud_vel_corr

    # Scale the data in order to view desired DIB

    DIB_plot_x = [total_wavelength_correction * x_val for x_val in DIB_x]
    DIB_plot_y_1 = [(y / max(DIB_y)) for y in DIB_y]
    a = np.mean(DIB_plot_y_1[:25])
    DIB_plot_y = DIB_plot_y_1 / a
    NaI_plot_x = [total_wavelength_correction * x_val for x_val in xpoint]
    NaI_plot_y = [y / max(ypoint) for y in ypoint]

    # write to .txt file
    os.chdir(writepath)
    f = open(target + "_subrange.txt", "w+")
    f.write("# " + target)
    f.write("\n")
    title = "# Wavelength(1/cm)      Relative Intensity \n"
    f.write(title)

    writer = csv.writer(f, delimiter="\t")
    writer.writerows(zip(DIB_plot_x, DIB_plot_y))

    f.close()


# fullCmdArguments = sys.argv
# args = fullCmdArguments[1:]
# arN = len(sys.argv)

# print(args)
# if len(args) != 5:
#     print('\nSyntax:    python fitsto2dtxt.py target, filepath, writepath, xmin, xmax\n')

# else:
#     fitstotxt(args[0], args[1], args[2], args[3], args[4])
