#!/usr/bin/env python

# +
#
# DESCRIPTION:
#        To simply plot any target in the EDIBLES DR3 and upper versions
#
# CALLING SEQUENCE:
#        ipython edibles_plt.py target_name obs_date wtype vel:wave xrng
#
#    Arguments:
#        target_name  ->  target Name
#        obs_date     ->  date of observation YYYYMMDD
#        wtype        ->  type of velocity correction
#                         geo: Geocentric correction
#                         hel: Heliocentric correction
#                         bar: Barycentric correction
#        vel:wave     ->  convert wavelength into velocity around 'wave'
#        xrng         ->  the wavelength range. should be in format of
#                         x:wi-wf e.g. x:5880-5910
#                         x:wc e.g. x:5780 will plot the spectrum around wc+-10A
#
#
#        dates        ->  plot multiple dates of the target
#                         all_dates
#        save         ->  save the plotted data to a text file
#                         One file will be created for each plot
#
#
# EXMPLE:
#     ALL DATES of HD170740 in velocity and apply the Barycentric correction
#     python edit_edibles_plt.py HD170740 bar x:7662-7669 vel:7667 all_dates
#
#     SINGLE date of HD170740 in velocity and apply the Barycentric correction
#     python edit_edibles_plt.py HD170740 20140915 bar x:7662-7669 vel:7667
#
#
# Author:
#        Amin Farhang
#
# Date:
#        Dec 2017
#
# COPY RIGHT:
#        This codes are produced to be used in EDIBLES project
#        PI: J. Cami, N. Cox please for use it first contact the
#        author or one of the PI's
#
# +

from __future__ import print_function
import os
import glob
import sys
import getopt
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt


# determine the DR directory
DR_address = '/data/FITS/'


# ****************
# read arguments *
# ****************
fullCmdArguments = sys.argv
args = fullCmdArguments[1:]
arN = len(sys.argv)


# check out the args
obs_date = ''
wtype = ''
wtypev = ''
xrng = ''
all_dates = False
save = False
for i in range(arN - 1):
    if args[i][:2] == 'HD':
        targetName = args[i]
    if args[i][:3] == '201':
        obs_date = args[i]
    if args[i][:2] == 'x:':
        xrng = args[i]
    if args[i] == 'geo':
        wtype = 'Geocentric'
    if args[i] == 'hel':
        wtype = 'Heliocentric'
    if args[i] == 'bar':
        wtype = 'Barycentric'
    if args[i][:3] == 'vel':
        wtypev = 'velocity'
        cen_vel = float(args[i].split(':')[1])
    if args[i] == 'all_dates':
        all_dates = True
    if args[i] == 'save':
        save = True



if len(xrng) != 0:
    if len(xrng.split('-')) == 2:
        xrng = xrng[2:]
        xmin_init = float(xrng.split('-')[0])
        xmax_init = float(xrng.split('-')[1])
    if len(xrng.split('-')) == 1:
        xrng = xrng[2:]
        xmin_init = float(xrng.split('-')[0]) - 2.0
        xmax_init = float(xrng.split('-')[0]) + 2.0
else:
    print('x range is not defined !')
    print('We are plotting the REDL 564-nm arm instead')
    xmin_init = 4617
    xmax_init = 5666

arms = ['BLUE_346', 'BLUE_437', 'REDL_564', 'REDU_564', 'REDL_860', 'REDU_860']

xmin = xmin_init
xmax = xmax_init

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
    os.chdir(DR_address)
    loc = DR_address + args[0] + '/' + warm[i]
    if warm[i] == 'REDL_564' or warm[i] == 'REDU_564':
        loc = DR_address + args[0] + '/RED_564'
    if warm[i] == 'REDL_860' or warm[i] == 'REDU_860':
        loc = DR_address + args[0] + '/RED_860'

    if os.path.isdir(loc) is False:
        print('This object have not been observed yet!')
        sys.exit()
    os.chdir(loc)

    # read wave and flux
    if warm[i] == 'BLUE_346' or warm[i] == 'BLUE_437':
        ar = 'B'
    if warm[i] == 'REDL_564' or warm[i] == 'REDL_860':
        ar = 'L'
    if warm[i] == 'REDU_564' or warm[i] == 'REDU_860':
        ar = 'U'
    obsNames = '%s*_%s.fits.gz' % (targetName, ar)

    n_science = 0
    all_names = []

    # choose the highest S/N
    for file in glob.glob(obsNames):
        n_new = float(file.split('_')[2][1:])
        if n_new > n_science:
            n_science = n_new
            trName = file
            all_names.append(trName)

    # look for a specific date
    if len(obs_date) != 0:
        if warm[i] == 'BLUE_346':
            w_ar = 'w346'
        if warm[i] == 'BLUE_437':
            w_ar = 'w437'
        if warm[i] == 'REDL_564' or warm[i] == 'REDU_564':
            w_ar = 'w564'
        if warm[i] == 'REDL_860' or warm[i] == 'REDU_860':
            w_ar = 'w860'
        all_names = glob.glob(targetName + '_' + w_ar +
                              '*' + obs_date + '_' + ar + '.fits.gz')
        trName = all_names[0]

    # plot all dates
    if all_dates is True:
        if warm[i] == 'BLUE_346':
            w_ar = 'w346'
        if warm[i] == 'BLUE_437':
            w_ar = 'w437'
        if warm[i] == 'REDL_564' or warm[i] == 'REDU_564':
            w_ar = 'w564'
        if warm[i] == 'REDL_860' or warm[i] == 'REDU_860':
            w_ar = 'w860'
        all_names = glob.glob(targetName + '_' + w_ar +
                              '*' + obs_date + '_' + ar + '.fits.gz')



    fig = plt.figure()
    for index, trName in enumerate(all_names):

        xmin = xmin_init
        xmax = xmax_init

        # read the wavelength solution
        os.chdir(loc)
        data = fits.open(trName)
        flux = data[0].data
        hdr = fits.getheader(trName).copy()
        CRVAL1 = hdr['CRVAL1']
        CDELT1 = hdr['CDELT1']
        wave_air = CRVAL1 + CDELT1 * np.arange(len(flux))
        w_air = []
        flx = []

        # Geocentric (diurnal) velocity
        if wtype == 'Geocentric':
            ra = hdr['RA']
            dec = hdr['DEC']
            lmst = hdr['LST']
            deg2rad = np.pi / 180.0
            obs_lat = hdr['HIERARCH ESO TEL GEOLAT']
            obs_lon = hdr['HIERARCH ESO TEL GEOLON']
            obs_alt = hdr['HIERARCH ESO TEL GEOELEV']
            dlat = -(11. * 60. + 32.743) * np.sin(2 * obs_lat * deg2rad) + \
                1.1633 * np.sin(4 * obs_lat * deg2rad) - 0.0026 * \
                np.sin(6 * obs_lat * deg2rad)
            lat = obs_lat + dlat / 3600.0
            r = 6378160.0 * (0.998327073 + 0.001676438 * np.cos(2 * lat * deg2rad) -
                             0.00000351 * np.cos(4 * lat * deg2rad) + 0.000000008 * np.cos(6 * lat * deg2rad)) + obs_alt
            v = 2. * np.pi * (r / 1000.) / (23.934469591229 * 3600.)
            v_diurnal = v * np.cos(lat * deg2rad) * np.cos(dec * deg2rad) * np.sin((ra - lmst * 15.0) * deg2rad)
            wave_air = wave_air * (1.0 + v_diurnal / 299792.458)


        # Heliocentric correction (in air)
        if wtype == 'Heliocentric':
            v_heli = hdr['HIERARCH ESO QC VRAD BARYCOR']
            wave_air = wave_air * (1.0 + v_heli / 299792.458)
            print(trName.split("_")[3] + ': ' + str(v_heli))

        # Barycentric correction (in air)
        if wtype == 'Barycentric':
            v_bari = hdr['HIERARCH ESO QC VRAD BARYCOR']
            wave_air = wave_air * (1.0 + v_bari / 299792.458)
            print(trName.split("_")[3] + ': ' + str(v_bari))

        idx = [i for i, e in enumerate(wave_air) if (e > xmin and e < xmax)]
        # print('len(idx) = ' + str(len(idx)))
        flx_new = []
        for j in idx:
            flx_new.append(flux[j])
        ymin = min(flx_new)
        ymax = max(flx_new)

        # velocity
        if wtypev == 'velocity':
            wave_air = ((wave_air - cen_vel) / cen_vel) * 299792.458

        # put them into one file
        for k in range(len(wave_air)):
            w_air.append(wave_air[k])
            flx.append(flux[k])


        # velocity range
        if wtypev == 'velocity':
            if len(xrng) != 0:
                xmin = ((xmin - cen_vel) / cen_vel) * 299792.458
                xmax = ((xmax - cen_vel) / cen_vel) * 299792.458
            else:
                xmin = (((cen_vel - 20) - cen_vel) / cen_vel) * 299792.458
                xmax = (((cen_vel + 20) - cen_vel) / cen_vel) * 299792.458


        os.chdir(os.getenv("HOME"))

        rows = int(np.ceil(len(all_names)/2.))
        if len(all_names) == 1:
            cols = 1
        else:
            cols = 2

        plt.subplot(rows, cols, index+1)
        plt.plot(w_air, flx)
        plt.xlim([xmin, xmax])
        plt.ylim([ymin - 0.1 * (ymax - ymin), ymax + 0.1 * (ymax - ymin)])
        xtitle = 'Wavelength'
        if wtypev == 'velocity':
            xtitle = 'Velocity (km/s)'
        # plt.xlabel(xtitle)
        plt.ylabel('Flux')
        plt.title(trName.split("_")[3])
        plt.ticklabel_format(useOffset=False)
        file_save = trName.split("_")[0] + '_' + trName.split("_")[3] + '.txt'
        data = np.column_stack((w_air, flux))

        # save subset of data to txt file
        
        if save is True:
            data_tosave = []
            for i, val in enumerate(data):
                if xmin < val[0] < xmax:
                    data_tosave.append(data[i])
            data_tosave = np.array(data_tosave)

            os.chdir('/export/home/klay/github/edibles/HD170740/txt/')
            np.savetxt(file_save, data_tosave)

fig.canvas.set_window_title(trName.split("_")[0])

plt.tight_layout(h_pad=-0.1)
plt.show()
