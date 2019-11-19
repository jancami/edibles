import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib import ticker

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from scipy import interpolate

import numpy as np
from astropy.io import fits

#font = {'family' : 'arial',
#        'weight' : 'normal',
#        'size'   : 16}

SMALL_SIZE = 6
MEDIUM_SIZE = 8
BIGGER_SIZE = 10


plt.style.use('seaborn-white')
import seaborn as sns
sns.set()

plt.ion()

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

#plt.rcParams['font.family'] = 'serif'
#plt.rcParams['font.serif'] = 'Ubuntu'
plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})


def vac2air_ciddor(vacw):
    """ Convert vacuum wavelengths in Angstroms to air wavelengths.

    This uses the relation from Ciddor 1996, Applied Optics LP,
    vol. 35, Issue 9, p.1566. Only valid for wavelengths > 2000 Ang.
    """
    k0 = 238.0185
    k1 = 1e-8 * 5792105.
    k2 = 57.362
    k3 = 1e-8 * 167917.
    s2 = (1e4 / vacw)**2
    n = 1 + k1/(k0 - s2) + k3/(k2 - s2)
    airw = vacw / n

    return airw


files = ['HD170740_346_rct.asc','HD170740_437x_rct.asc',
         'HD170740_564L_rct.asc','HD170740_564U_rct.asc',
         'HD170740_860L_rct.asc','HD170740_860U_rct.asc']

vbary = [8.057485, -28.454783, 8.041374, 8.041374, 2.17133, 2.17133]

#location_dib_spectrum = "/Users/nick/Dropbox/EDIBLES/DATA_RELEASES/AUXILLARY_DATA/DIB_Catalogs/"
location_dib_spectrum = "/Users/nick/DropboxOffline/EDIBLES/python4github/DIB_Catalogs/"
wave_dibspec,dibspec=np.loadtxt(location_dib_spectrum+"DIBspectrum_JD94_AverageLOS_step0.01.dat",unpack=True) #4000-9700

location_dib_spectrum = "/Users/nick/DropboxOffline/EDIBLES/python4github/DIB_Catalogs/"
wave_dibspec2,dibspec2=np.loadtxt(location_dib_spectrum+"DIBspectrum_Hobbs2009_HD183143.dat",unpack=True) #4400-8800

location_sky="/Users/nick/DropboxOffline/EDIBLES/python4github/SkyTransmission/" ## adjust if necessary.
wt,ft=np.loadtxt(location_sky+"transmission_300to1000nm.dat", unpack=True)
wtc = vac2air_ciddor(wt*10.0)

#### START FIGURE ####

fig = plt.figure(num=None, figsize=(5,5), dpi=150, facecolor='w', edgecolor='k')
gs1 = GridSpec(3, 3, left=0.05, right=0.95, top=0.975, bottom=0.05, wspace=0.1, hspace=0.1)

ax_main = fig.add_subplot(gs1[0,:])
ax_middle = fig.add_subplot(gs1[1,:])

### main panel (top)
ax_main.cla()
ax_main.set_xlim(3100,10550)
ax_main.set_ylim(0.5,1.1)
#ax_main.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=True,
#                     bottom=True, top=True, left=False, right=False)
ax_main.tick_params(axis="y",direction="in", pad=2, labelsize=6, labelleft=True, labelright=True)
ax_main.tick_params(axis="x",direction="in", pad=2, labelsize=6)

xmajloc = ticker.MultipleLocator(base=1000.0)
xminloc = ticker.MultipleLocator(base=200.0)
ax_main.yaxis.set_major_locator(xmajloc)
ax_main.yaxis.set_minor_locator(xminloc)

ymajloc = ticker.MultipleLocator(base=0.20)
yminloc = ticker.MultipleLocator(base=0.05)
ax_main.yaxis.set_major_locator(ymajloc)
ax_main.yaxis.set_minor_locator(yminloc)

for s,v in zip(files,vbary):
    print s,v
    w1,n=np.loadtxt(s,unpack=True)
    w = w1+(v+13.0)/299792.458*w1
    ax_main.plot(w, n, color='Blue', linewidth=0.5)
idx=np.where( (wtc > 4000) & ( wtc < 11000) )

"""
wtc2 = wtc[idx]+(v+13.0)/299792.458*wtc[idx]
ax_main.plot(wtc2,(ft[idx]/np.nanmax(ft[idx]))**2.0, color='Orange', linestyle='solid', label="Telluric", linewidth=1, alpha=0.7)
"""

ax_main.text(3150, 1.05, "HD170740", size=6, rotation=0.,
         ha="left", va="center",
         bbox=dict(boxstyle="square", ec="black", fc="white", alpha=0.8) )

ax_main.axvspan(5772,5802, alpha=0.75, color='gray')
ax_main.axvspan(6610,6615, alpha=0.75, color='gray')
ax_main.axvspan(6372,6382, alpha=0.75, color='gray')

ax_main.axvspan(9565, 9645, alpha=0.75, color='gray')

#idx=np.where( (wave_dibspec2 > 6100) & ( wave_dibspec2 < 6300) )
ax_main.plot(wave_dibspec2,(dibspec2**0.5), linewidth=2.0, alpha=0.9, linestyle='solid', color='Green', label="")


### middle panel ###
ax_middle.cla()
ax_middle.set_xlim(9565,9645)
ax_middle.set_ylim(0.8,1.02)
ax_middle.yaxis.tick_right()
ymajloc = ticker.MultipleLocator(base=0.10)
yminloc = ticker.MultipleLocator(base=0.02)
ax_middle.yaxis.set_major_locator(ymajloc)
ax_middle.yaxis.set_minor_locator(yminloc)
#ax_middle.yaxis.set_label_position("right")
ax_middle.tick_params(axis='y', direction='out', pad=2, labelright=True, labelleft=True)
ax_middle.tick_params(axis="x", direction="in", pad=2)

# plot spectrum
s=files[5]
v=vbary[5]
w1,n=np.loadtxt(s,unpack=True)
w = w1+(v+13.0)/299792.458*w1
ax_middle.plot(w, n, color='Blue', linewidth=1.0)

ax_middle.text(9577, 0.91, "C60+", size=6, rotation=0.,
         ha="center", va="center",
         bbox=dict(boxstyle="square", ec="black", fc="white", alpha=0.8) )

ax_middle.text(9632, 0.91, "C60+", size=6, rotation=0.,
         ha="center", va="center",
         bbox=dict(boxstyle="square", ec="black", fc="white", alpha=0.8) )

# plot telluric
idx=np.where( (wtc > 9400) & ( wtc < 9700) )
wtc2 = wtc[idx]+(v+13.0)/299792.458*wtc[idx]
ax_middle.plot(wtc2,(ft[idx]/np.nanmax(ft[idx]))**2.0, color='Orange', linestyle='solid', label="Telluric", linewidth=0.8, alpha=0.9)

# plot DIB synthetic
idx=np.where( (wave_dibspec > 9500) & ( wave_dibspec < 9900) )
ax_middle.plot(wave_dibspec[idx],(dibspec[idx]**0.3), linewidth=1.0, alpha=0.9, linestyle='solid', color='Green', label="")

#### bottom panels
ax_bottom1 = fig.add_subplot(gs1[2,0])
ax_bottom2 = fig.add_subplot(gs1[2,2])
ax_bottom3 = fig.add_subplot(gs1[2,1])

"""
# plot telluric
idx=np.where( (wtc > 5750) & ( wtc < 6650) )
wtc2 = wtc[idx]+(v+13.0)/299792.458*wtc[idx]
ax_bottom1.plot(wtc2,(ft[idx]/0.885), color='Orange', linestyle='solid', label="Telluric", linewidth=0.8, alpha=0.9)
ax_bottom2.plot(wtc2,(ft[idx]/0.93), color='Orange', linestyle='solid', label="Telluric", linewidth=0.8, alpha=0.9)
ax_bottom3.plot(wtc2,(ft[idx]/0.92), color='Orange', linestyle='solid', label="Telluric", linewidth=0.8, alpha=0.9)
"""

### bottom1 ###
ax_bottom1.cla()
ax_bottom1.set_xlim(5772,5803)
ax_bottom1.set_ylim(0.83,1.03)
ax_bottom1.yaxis.set_label_position("left")
ax_bottom1.tick_params(axis='y', direction='out', pad=2)
ax_bottom1.tick_params(axis="x",direction="in", pad=2)
s=files[3]
v=vbary[3]

ymajloc = ticker.MultipleLocator(base=0.05)
yminloc = ticker.MultipleLocator(base=0.01)
ax_bottom1.yaxis.set_major_locator(ymajloc)
ax_bottom1.yaxis.set_minor_locator(yminloc)
ymajloc = ticker.MultipleLocator(base=10.0)
yminloc = ticker.MultipleLocator(base=2.0)
ax_bottom1.xaxis.set_major_locator(ymajloc)
ax_bottom1.xaxis.set_minor_locator(yminloc)

#ax_bottom1.plot(wave_dibspec2,(dibspec2**0.3), linewidth=2.0, alpha=0.9, linestyle='solid', color='Green', label="")

w1,n=np.loadtxt(s,unpack=True)
w = w1+(v+13.0)/299792.458*w1
ax_bottom1.plot(w, n, color='Blue', linewidth=1.0)

ax_bottom1.text(5780, 0.85, "5780 DIB", size=6, rotation=0.,
         ha="center", va="center",
         bbox=dict(boxstyle="square", ec="black", fc="white", alpha=0.8) )

ax_bottom1.text(5797, 0.85, "5797 DIB", size=6, rotation=0.,
         ha="center", va="center",
         bbox=dict(boxstyle="square", ec="black", fc="white", alpha=0.8) )

#ax_bottom1.plot(wsigma,fsigma/consigma-0.04, color='Gray', linewidth=1.0)
#ax_bottom1.annotate(r"$\sigma$", [5775,0.93])

#ax_bottom1.plot(wzeta,fzeta/conzeta, color='Blue', linewidth=1.0)
#ax_bottom1.annotate(r"$\zeta$", [5775,1.01])

### bottom3 ###

ax_bottom3.cla()
ax_bottom3.set_xlim(6374.01,6381.99)
ax_bottom3.set_ylim(0.83,1.03)
ax_bottom3.yaxis.set_label_position("left")
ax_bottom3.tick_params(axis='y', direction='out', pad=2, labelright=False, labelleft=False)
ax_bottom3.tick_params(axis="x",direction="in", pad=2)

ymajloc = ticker.MultipleLocator(base=0.05)
yminloc = ticker.MultipleLocator(base=0.01)
ax_bottom3.yaxis.set_major_locator(ymajloc)
ax_bottom3.yaxis.set_minor_locator(yminloc)
xmajloc = ticker.MultipleLocator(base=2.0)
xminloc = ticker.MultipleLocator(base=0.5)
ax_bottom3.xaxis.set_major_locator(xmajloc)
ax_bottom3.xaxis.set_minor_locator(xminloc)

#ax_bottom3.plot(wave_dibspec2,(dibspec2**0.3), linewidth=2.0, alpha=0.9, linestyle='solid', color='Green', label="")

s=files[3]
v=vbary[3]
w1,n=np.loadtxt(s,unpack=True)
w = w1+(v+13.0)/299792.458*w1
ax_bottom3.plot(w, n, color='Blue', linewidth=1.0)

ax_bottom3.text(6376, 0.85, "6376 DIB", size=6, rotation=0.,
         ha="center", va="center",
         bbox=dict(boxstyle="square", ec="black", fc="white", alpha=0.8) )

ax_bottom3.text(6379, 0.85, "6379 DIB", size=6, rotation=0.,
         ha="center", va="center",
         bbox=dict(boxstyle="square", ec="black", fc="white", alpha=0.8) )


### bottom2 ###
s=files[3]
v=vbary[3]

ax_bottom2.cla()
ax_bottom2.set_xlim(6611.5,6615.99)
ax_bottom2.set_ylim(0.83,1.03)
ax_bottom2.yaxis.set_label_position("left")
#ax_bottom2.tick_params(axis='y', direction='out', pad=2, labelright=True, labelleft=False)
ax_bottom2.tick_params(axis="x",direction="in", pad=2)

ymajloc = ticker.MultipleLocator(base=0.05)
yminloc = ticker.MultipleLocator(base=0.01)
ax_bottom2.yaxis.set_major_locator(ymajloc)
ax_bottom2.yaxis.set_minor_locator(yminloc)
xmajloc = ticker.MultipleLocator(base=1.0)
xminloc = ticker.MultipleLocator(base=0.5)
ax_bottom2.xaxis.set_major_locator(xmajloc)
ax_bottom2.xaxis.set_minor_locator(xminloc)

#ax_bottom2.plot(wave_dibspec2,(dibspec2**0.3), linewidth=2.0, alpha=0.9, linestyle='solid', color='Green', label="")

w1,n=np.loadtxt(s,unpack=True)
w = w1+(v+13.0)/299792.458*w1
ax_bottom2.plot(w, n, color='Blue', linewidth=1.0)

ax_bottom2.text(6613.5, 0.85, "6613 DIB", size=6, rotation=0.,
         ha="center", va="center",
         bbox=dict(boxstyle="square", ec="black", fc="white", alpha=0.8) )

### SAVING TO PDF!!!
plt.savefig("panelplot.pdf", dpi=600)

##### END #####
