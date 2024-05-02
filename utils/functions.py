from __future__ import print_function

import os
import csv
import numpy as np
from scipy.signal import find_peaks
import astropy.constants as cst
from astropy.io import fits
import astropy.units as u
from heapq import nsmallest


def barycorrectSpectrum(wave_array, v_bary):
    """Barycentric wavelength correction tool.

    :param wave_array: The wavelength array of the spectrum
    :type wave_array: ndarray
    :param v_bary: The barycentric velocity
    :type v_bary: float64

    :return: The corrected wavelength array
    :rtype: ndarray

    """
    wave_array = wave_array + (v_bary / cst.c.to("km/s").value) * wave_array

    return wave_array


def fitstotxt(target, filepath, writepath, xmin, xmax):
    """
    Reads a FITS file, applies barycentric and cloud corrections
    and writes a subsection within a given range to a .txt file.

    :param target: Name of Star
    :type target: str
    :param filepath: Data storage location
    :type filepath: str
    :param writepath: Location of wanted txt file
    :type writepath: str
    :param xmin: Minimum wavelength
    :type xmin: float64
    :param xmax: Maximum Wavelength
    :type xmax: float64

    :return: spectrum within given range, tab separated, in a given folder.
    :rtype: txt file

    :Example:

    >>> fitstotxt('HD144470', '/data/DR3_fits/', '/home/txtfiles/, 6610, 6618)
    # HD144470
    # Wavelength(1/cm)      Relative Intensity
    3303.1797537322504  1.0084753
    3303.199755852153   1.0125096
    3303.219757972055   1.012984
    3303.2397600919576  1.0122045

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

    # correction for wavelenth taking into account
    # barycentric velocity and cloud velocity.
    total_wavelength_correction = bary_corr * cloud_vel_corr

    # Scale the data in order to view desired DIB

    DIB_plot_x = [total_wavelength_correction * x_val for x_val in DIB_x]
    DIB_plot_y_1 = [(y / max(DIB_y)) for y in DIB_y]
    a = np.mean(DIB_plot_y_1[:25])
    DIB_plot_y = DIB_plot_y_1 / a

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
    #     print('\nSyntax:    python fitsto2dtxt.py target, filepath,
    #                   writepath, xmin, xmax\n')

    # else:
    #     fitstotxt(args[0], args[1], args[2], args[3], args[4])


def make_grid(lambda_start, lambda_end, resolution=None, oversample=None):

    # check keywords
    if oversample is None:
        oversample = 40.0
    if resolution is None:
        resolution = 1500.0

    lambda_start = float(lambda_start)
    lambda_end = float(lambda_end)

    # produce grid
    R = resolution * oversample
    n_points = (
        round(
            (np.log(lambda_end / lambda_start)) / (np.log(-(1 + 2 * R) / (1 - 2 * R)))
        )
        + 1
    )
    f = -(1 + 2 * R) / (1 - 2 * R)
    factor = f ** np.arange(n_points)

    wave = np.full(int(n_points), lambda_start, dtype=float)
    grid = wave * factor

    return grid


def param_convert(params):
    """
    Function to convert voigt parameteres from astronomical
    to mathematical version.

    :param params: (cent, b_eff, Gamma, scaling)
    :type params: tuple

    :return: converted parameters
    :rtype: tuple

    """

    cent, b_eff, Gamma, scaling = params

    cent = cent * u.AA
    b_eff = b_eff * u.km / u.s
    Gamma = Gamma * 1 / u.s

    # cent
    nu0 = cent.to(u.Hz, equivalencies=u.spectral())
    # b_eff
    delta_nu_D = nu0 * b_eff / cst.c.to("km/s")  # freq * km/s / km/s
    sigma = delta_nu_D / np.sqrt(2)
    alpha_Hz = sigma * np.sqrt(2 * np.log(2))
    alpha = alpha_Hz * cst.c / (nu0 ** 2)  # Hz * m/s  / Hz^2
    alpha = alpha.decompose().to(u.AA)
    # Gamma
    gam = Gamma / (4 * np.pi) * 1e-10

    converted_params = (cent.value, alpha.value, gam.value, scaling)
    return converted_params


def peak_wavelength_largest(wave, flux, n=1):
    """
    Function that returns the wavelengths of the lowest n peaks of the equation

    :param wave: wavelength grid
    :type wave: ndarray
    :param flux: flux grid
    :type flux: ndarray
    :param n: number of peaks to return, default=1
    :type n: int

    :return: wavelengths of lowest n peaks
    :rtype: list

    """
    peaks, _ = find_peaks(-flux)  # indices of peaks
    peak_flux = nsmallest(n, flux[peaks])  # smallest n flux values at peaks
    peak_wavelength = [
        wave[np.where(flux == peak)][0] for peak in peak_flux
    ]  # corresponding wavelengths
    return peak_wavelength


def peak_wavelength_all_prominent(wave, flux, prominence):
    """
    Function that returns the wavelengths of the lowest n peaks of the equation

    :param wave: wavelength grid
    :type wave: ndarray
    :param flux: flux grid
    :type flux: ndarray
    :param prominence: prominence of points (see scipy.signal.find_peaks)
    :type prominence: float, between 0 and 1

    :return: wavelengths of all prominent peaks
    :rtype: list

    """

    yrange = max(flux) - min(flux)
    peaks, _ = find_peaks(-flux, prominence=prominence * yrange)
    return [wave[peak] for peak in peaks]


def printHeader(input_fits):
    """
    A function to print out the header of a FITS file

    :param input_fits: Path to FITS file starting from DATADIR
    :type input_fits: str

    """
    from edibles import DATADIR

    path = DATADIR + input_fits

    hdu = fits.open(path)
    print(hdu[0].header)


def read_array(filename, dtype, separator=","):
    """ Read a file with an arbitrary number of columns.
         The type of data in each column is arbitrary
         It will be cast to the given dtype at runtime
     """
    cast = np.cast
    data = [[] for dummy in range(len(dtype))]
    for line in open(filename, "r"):
        fields = line.strip().split(separator)
        for i, number in enumerate(fields):
            data[i].append(number)
    for i in range(len(dtype)):
        data[i] = cast[dtype[i]](data[i])
    return np.rec.array(data, dtype=dtype)


def read_line_catalog(input_catalog):
    """
    (xpos_atoms, labels_atoms) =
        read_line_catalog('auxilarary_data/line_catalogs/edibles_linelist_atoms.csv')
    """

    x = np.genfromtxt(
        input_catalog, dtype=None, skip_header=1, delimiter=",", unpack=True
    )
    xpos = x["f1"]
    labels = x["f2"]

    """ alternative:
    xpos,labels = np.genfromtxt(input_catalog, usecols=(0,1), skip_header=1,
                                dtype=None, delimiter=',')
    """

    return xpos, labels


def read_sky_transmission(transmission_dat, scale_factor=1.0):

    import vac2air_ciddor

    """ transmission input spectrum in nm
    """
    wt, ft = np.loadtxt(transmission_dat, unpack=True)
    wtc = vac2air_ciddor(wt * 10.0)
    ft = ft ** scale_factor
    return wtc, ft


def smooth(y, box_pts):
    import numpy as np

    box = np.ones(box_pts) / box_pts
    y_smooth = np.convolve(y, box, mode="same")
    return y_smooth


def parseTextFile(file_name, delimiter=",", header=0):
    """ Parse a text file to a list. The file contents are delimited
    and have a header.


    :param file_name: The path to the file
    :type file_name: str
    :param delimiter: The delimiter to use to parse the file
    :type delimiter: str
    :param header: The number of lines at the top of the file to ignore
    :type header: int

    :return: Text file parsed into a list
    :rtype: list

    """

    with open(file_name) as f:

        # Skip the header
        for i in range(header):
            next(f)

        data = []

        # Parse file contents
        for line in f:

            # Remove the newline char
            line = line.replace("\n", "").replace("\r", "")

            # Split the line by the delimiter
            line = line.split(delimiter)

            # Strip whitespaces from individual entries in the line
            for i, entry in enumerate(line):
                line[i] = entry.strip()

            # Add the contents of the line to the data list
            data.append(line)

        return data


def vac2air_ciddor(vacw):
    """ Convert vacuum wavelengths in Angstroms to air wavelengths.

    This uses the relation from Ciddor 1996, Applied Optics LP,
    vol. 35, Issue 9, p.1566. Only valid for wavelengths > 2000 Ang.

    example:

    wt,ft=np.loadtxt("transmission.dat", unpack=True)
    wtc = vac2air_ciddor(wt*10.0) #make sure wave is in Angstroms.

    """
    k0 = 238.0185
    k1 = 1e-8 * 5792105.0
    k2 = 57.362
    k3 = 1e-8 * 167917.0
    s2 = (1e4 / vacw) ** 2
    n = 1 + k1 / (k0 - s2) + k3 / (k2 - s2)
    airw = vacw / n

    return airw


def vac2air_morton(vacw):
    """ Convert vacuum wavelengths in Angstroms to air wavelengths.

    This uses the relation from Morton 1991, ApJS, 77, 119. Only valid
    for wavelengths > 2000 Ang.  Use this for compatibility with older
    spectra that may have been corrected using the (older) Morton
    relation.  The Ciddor relation used in vac2air_ciddor() is claimed
    to be more accurate at IR wavelengths.

    example:
    wtc = vac2air_morton(wt*10.0) #make sure wave is in Angstroms.

    """
    temp = (1e4 / vacw) ** 2
    airw = (
        1.0
        / (1.0 + 6.4328e-5 + 2.94981e-2 / (146 - temp) + 2.5540e-4 / (41 - temp))
        * vacw
    )
    return airw


def write_spectrum_ascii(output_name, x, y, yerr, header):
    f = open(output_name, "w")
    if header is not None:
        f.write(header)
    for i in range(len(x) - 1):
        if yerr is not None:
            # print x
            f.write(
                str("%6.4f" % x[i])
                + "\t"
                + str("%7.4f" % y[i])
                + "\t"
                + str("%7.4f" % yerr[i])
                + "\n"
            )
        else:
            # print y
            f.write(str("%6.4f" % x[i]) + "\t" + str("%7.4f" % y[i]) + "\n")
    f.close()
    return
