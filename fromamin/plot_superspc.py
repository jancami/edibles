import numpy as np
from astropy.io import fits
from scipy import interpolate
from continuum_sigma import continuum_sigma
import matplotlib.pylab as plt

#loc_dr4 = '/Users/amin/DR4/'
loc_dr4 = '/data/DIBs/data/EDIBLES/DR3/FITS/'
loc_ord = '/Users/amin/Documents/UV_DIB/finding/order_ascii/'
colors = ['#00014C', '#00017F', '#000099', '#0001FF', '#5E05E5', '#8805E4', '#B104E3', '#DB04E2', '#E103BE', '#E00393', '#DF0268', '#DE023E']


# ==========================================================
#     average spectrum - this code average multi spectrums
# ==========================================================
def average_spec(wave, flux):

    # find the min and max for grid intervals
    if len(wave) > 1:
        minv = min(wave[0])
        maxv = max(wave[0])
        for loop in range(len(wave)):
            tmp_wvl = wave[loop]
            tmp_flx = flux[loop]
            if min(tmp_wvl) > minv: minv = min(tmp_wvl)
            if max(tmp_wvl) < maxv: maxv = max(tmp_wvl)
    else:
        minv = min(wave[0])
        maxv = max(wave[0])


    # make a grid
    grid = np.arange(minv+0.02, maxv-0.02, 0.02)

    # interpolate all spectra into a portion of grid
    wave_gridded = []
    flux_gridded = []
    sigma_grid = []
    for lop in range(len(wave)):
        tmp_wvl = wave[lop]
        tmp_flx = flux[lop]
        f = interpolate.interp1d(wave[lop], np.array(flux[lop]), kind='cubic') #
        flux_temp = f(grid)

        # print lop
        sigma = continuum_sigma(grid, flux_temp, plot=None)
        flux_gridded.append(flux_temp)
        sigma_grid.append(1.0/sigma)  # (np.std(flux_temp))  sigma**2


    wave_av = grid.copy()
    flux_av = np.average(flux_gridded, weights=sigma_grid, axis=0)

    return wave_av, flux_av


# ==========================================================
#     average spectrum - this code average multi spectrums
# ==========================================================
def averaging(target_list):
    lst_wvl_geo = []
    lst_wvl_str = []
    lst_wvl_ism = []
    lst_flx = []
    for loop in range(len(target_list)):
        # find index
        star_name = target_list[loop]
        idx = [i for i, xx in enumerate(targets) if xx==star_name][0]

        # read observing date
        date = obs_date[idx]
        fits_name = loc_dr4 + star_name + '/BLUE_346/' + star_name + '_w346_blue_' + date + '_' + order + '.fits'
        ascii_name = loc_dr4 + star_name + '/BLUE_346/' + star_name + '_w346_blue_' + date + '_' + order + '.ascii'

        # read baricentric correction from header
        data = fits.open(fits_name)
        hdr = data[0].header
        bari_corr = hdr['HIERARCH ESO QC VRAD BARYCOR']

        # read spectrum
        wave, flux = np.loadtxt(ascii_name, unpack=True)

        # normalize
        z = np.polyfit(wave, flux, 2)
        p = np.poly1d(z)
        continuum_fit = p(wave)
        flux = flux / continuum_fit

        # move to other references
        wave_bari = wave * (1.0 + bari_corr/299792.458)
        wave_str  = wave_bari * (1.0 - vel_str[idx]/299792.458)
        wave_ism  = wave_bari * (1.0 - vel_ism[idx]/299792.458)

        # append to list
        lst_wvl_geo.append(wave)
        lst_wvl_str.append(wave_str)
        lst_wvl_ism.append(wave_ism)
        lst_flx.append(flux)

    return lst_wvl_geo, lst_wvl_str, lst_wvl_ism, lst_flx




# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#     target based on E(B-V)
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
targets = ['HD147889', 'HD167971', 'HD112272', 'HD186841', 'HD73882' , 'HD149404', 'HD152424', 'HD148937',
           'HD185859', 'HD80558' , 'HD167838', 'HD153919', 'HD124314', 'HD147888', 'HD170740', 'HD185418',
           'HD150136', 'HD147933', 'HD41117' , 'HD37367' , 'HD27778' , 'HD147165', 'HD172694', 'HD37903' ,
           'HD93205' , 'HD149757', 'HD151804', 'HD203532', 'HD24398' , 'HD37022' , 'HD114886', 'HD23180' ,
           'HD57061' , 'HD66811' , 'HD180554', 'HD40111' , 'HD93030' , 'HD81188' , 'HD164353', 'HD36695' ,
           'HD157246', 'HD23016' , 'HD183143', 'HD169454']

obs_date = ['20150816', '20140921', '20170501', '20160910', '20170425', '20170418', '20160411', '20150817',
            '20150920', '20150603', '20160808', '20160707', '20170501', '20140924', '20150424', '20140923',
            '20150628', '20170615', '20160222', '20141011', '20141029', '20150531', '20150816', '20140929',
            '20150227', '20150720', '20150720', '20140923', '20140923', '20140924', '20170503', '20141029',
            '20141010', '20141010', '20160612', '20141029', '20141210', '20170418', '20170615', '20150822',
            '20150720', '20140920', '20170418', '20170418']

vel_ism = [-8.896, -13.435,    1.997, -11.619,  20.153,  -0.726,  -11.619,  -16.159,
           -9.488,  18.337,  -14.343,  -4.357,   1.089,  -9.804,  -10.712,  -11.619,
          -13.435,  -8.896,    12.89,  15.614,  14.706,  -6.173,   -7.988,    26.60,
            5.628, -14.851,    0.181,   12.89,   12.89,  19.245,   38.309,   11.983,
           31.954, -13.435, -13.4355,  13.798,   6.536,   7.444,  -15.251,   27.415,
            0.181,  17.429,  -20.153,  -22.695]

vel_str = [ -19.98, 13.72,  -23.21,  20.15,  36.20,  -109.09,  -44.07,   10.51,
              4.09, 20.15,   -6.34, -67.35,   5.10,   -17.58,  -20.79,  -10.82,
            -36.53,  3.40,   10.51,  44.23,  12.12,    25.76,    8.10,   15.33,
            -62.53, 72.33,  -48.89,   2.48,  17.74,    25.76,   55.47,   61.00,
             88.38, 15.33,  -14.36,  -7.14,   4.89,    28.17,   -2.32,  -74.57,
             -2.32,  9.71, -26.591,  1.703]



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# read superspectrums
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
redName = 'superspec_red.ascii'
bluName = 'superspec_blue.ascii'
geoName = 'superspec_geo.ascii'
strName = 'superspec_str.ascii'
xred, yred = np.loadtxt(redName, unpack=True)
xblue, yblue = np.loadtxt(bluName, unpack=True)
xgeo, ygeo = np.loadtxt(geoName, unpack=True)
xstr, ystr = np.loadtxt(strName, unpack=True)



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# coadding orders
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
order = 'O1'
vetical_line = 3078

vlflag = False
if vetical_line is not None:
    vlflag = True

red_list = ['HD167971','HD112272','HD186841','HD73882','HD149404','HD152424','HD148937','HD185859','HD80558','HD167838','HD153919']
blue_list = ['HD180554', 'HD40111', 'HD93030', 'HD81188', 'HD164353', 'HD36695', 'HD157246', 'HD23016']

if order == 'O1':
    blue_list = ['HD180554', 'HD81188', 'HD36695', 'HD157246']
    red_list = ['HD73882','HD149404','HD152424','HD148937','HD185859','HD80558','HD167838','HD153919']

if order == 'O2' or order == 'O3':
    blue_list = ['HD180554', 'HD81188', 'HD164353', 'HD36695', 'HD157246']
    red_list = ['HD112272','HD186841','HD73882','HD149404','HD152424','HD148937','HD185859','HD80558','HD167838','HD153919']

if order == 'O4' or order == 'O5' or order == 'O6' or order == 'O7':
    red_list = ['HD112272','HD186841','HD73882','HD149404','HD152424','HD148937','HD185859','HD80558','HD167838','HD153919']
    blue_list = ['HD180554', 'HD81188', 'HD164353', 'HD36695', 'HD157246']

if order == 'O8' or order == 'O9':
    blue_list = ['HD180554', 'HD81188', 'HD164353', 'HD36695', 'HD157246']

if order == 'O32':
    red_list = ['HD167971','HD167971','HD112272','HD186841','HD73882','HD149404','HD152424','HD148937','HD185859','HD80558','HD167838','HD153919']



# superspectrum for blue sightlines
lst_wvl_geo_blue, lst_wvl_str_blue, lst_wvl_ism_blue, lst_flx_blue = averaging(blue_list)
super_wvl_ism_blue, super_flx_ism_blue = average_spec(lst_wvl_ism_blue, lst_flx_blue)

# superspectrum for red sightlines
lst_wvl_geo_red, lst_wvl_str_red, lst_wvl_ism_red, lst_flx_red = averaging(red_list)
super_wvl_geo_red, super_flx_geo_red = average_spec(lst_wvl_geo_red, lst_flx_red)
super_wvl_str_red, super_flx_str_red = average_spec(lst_wvl_str_red, lst_flx_red)
super_wvl_ism_red, super_flx_ism_red = average_spec(lst_wvl_ism_red, lst_flx_red)






# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Blue superspectrum
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig = plt.figure()
gs = plt.GridSpec(3, 2)
ax1 = fig.add_subplot(gs[:4])
ax2 = fig.add_subplot(gs[4:])
ax1.title.set_text('Blue individual sightlines')

for loop in range(len(lst_wvl_ism_blue)):
    x = lst_wvl_ism_blue[loop]
    y = lst_flx_blue[loop]
    if loop == 0:
        miny = min(y) - 0.05
        xmin = min(x)
        xmax = max(x)
    if loop == len(lst_wvl_ism_blue)-1:
        maxy = max(y) + loop*0.1 + 0.03
    if min(x) >= xmin: xmin = min(x)
    if max(x) <= xmax: xmax = max(x)
    ax1.plot(x, y + loop*0.1, color=colors[loop])

ymin = min(lst_flx_blue[0][( (lst_wvl_ism_blue[0]>=xmin) & (lst_wvl_ism_blue[0]<=xmax) )])

ax2.plot(super_wvl_ism_blue, super_flx_ism_blue, color='gray')
ax2.plot(xblue, yblue + 0.03, color='blue')
ax1.set_xlim(xmin, xmax)
ax1.set_ylim(ymin, maxy)
ax2.set_xlim(xmin, xmax)
ax2.set_ylim(min(super_flx_ism_blue)-0.02, max(yblue) + 0.04)



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Geo superspectrum
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig = plt.figure()
gs = plt.GridSpec(3, 2)
ax1 = fig.add_subplot(gs[:4])
ax2 = fig.add_subplot(gs[4:])
ax1.title.set_text('Geo individual sightlines')

for loop in range(len(lst_wvl_geo_red)):
    x = lst_wvl_geo_red[loop]
    y = lst_flx_red[loop]
    if loop == 0:
        miny = min(y) - 0.05
        xmin = min(x)
        xmax = max(x)
    if loop == len(lst_wvl_geo_red)-1:
        maxy = max(y) + loop*0.1 + 0.03
    if min(x) >= xmin: xmin = min(x)
    if max(x) <= xmax: xmax = max(x)
    ax1.plot(x, y + loop*0.1, color=colors[loop])

ymin = min(lst_flx_red[0][( (lst_wvl_geo_red[0]>=xmin) & (lst_wvl_geo_red[0]<=xmax) )])

ax2.plot(super_wvl_geo_red, super_flx_geo_red, color='gray')
ax2.plot(xgeo, ygeo + 0.03, color='green')
ax1.set_xlim(xmin, xmax)
ax1.set_ylim(ymin, maxy)
ax2.set_xlim(xmin, xmax)
ax2.set_ylim(min(super_flx_geo_red)-0.02, max(ygeo) + 0.04)



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# stellar superspectrum
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig = plt.figure()
gs = plt.GridSpec(3, 2)
ax1 = fig.add_subplot(gs[:4])
ax2 = fig.add_subplot(gs[4:])
ax1.title.set_text('Stellar individual sightlines')

for loop in range(len(lst_wvl_str_red)):
    x = lst_wvl_str_red[loop]
    y = lst_flx_red[loop]
    if loop == 0:
        miny = min(y) - 0.05
        xmin = min(x)
        xmax = max(x)
    if loop == len(lst_wvl_str_red)-1:
        maxy = max(y) + loop*0.1 + 0.03
    if min(x) >= xmin: xmin = min(x)
    if max(x) <= xmax: xmax = max(x)
    ax1.plot(x, y + loop*0.1, color=colors[loop])

ymin = min(lst_flx_red[0][( (lst_wvl_str_red[0]>=xmin) & (lst_wvl_str_red[0]<=xmax) )])

ax2.plot(super_wvl_str_red, super_flx_str_red, color='gray')
ax2.plot(xstr, ystr + 0.03, color='brown')
ax1.set_xlim(xmin, xmax)
ax1.set_ylim(ymin, maxy)
ax2.set_xlim(xmin, xmax)
ax2.set_ylim(min(super_flx_str_red)-0.02, max(ystr) + 0.04)




# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Red superspectrum
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig = plt.figure()
gs = plt.GridSpec(3, 2)
ax1 = fig.add_subplot(gs[:4])
ax2 = fig.add_subplot(gs[4:])
ax1.title.set_text('Red individual sightlines')

for loop in range(len(lst_wvl_ism_red)):
    x = lst_wvl_ism_red[loop]
    y = lst_flx_red[loop]
    if loop == 0:
        miny = min(y) - 0.05
        xmin = min(x)
        xmax = max(x)
    if loop == len(lst_wvl_ism_red)-1:
        maxy = max(y) + loop*0.1 + 0.03
    if min(x) >= xmin: xmin = min(x)
    if max(x) <= xmax: xmax = max(x)
    ax1.plot(x, y + loop*0.1, color=colors[loop])

ymin = min(lst_flx_red[0][( (lst_wvl_ism_red[0]>=xmin) & (lst_wvl_ism_red[0]<=xmax) )])

ax2.plot(super_wvl_ism_red, super_flx_ism_red, color='gray')
ax2.plot(xred, yred + 0.03, color='red')
if vlflag is True:
    ax1.plot([vetical_line, vetical_line], [0,2], linestyle='--', color='green')
    ax2.plot([vetical_line, vetical_line], [0,2], linestyle='--', color='green')
ax1.set_xlim(xmin, xmax)
ax1.set_ylim(ymin, maxy)
ax2.set_xlim(xmin, xmax)
ax2.set_ylim(min(super_flx_ism_red)-0.02, max(yred) + 0.04)




plt.show()
