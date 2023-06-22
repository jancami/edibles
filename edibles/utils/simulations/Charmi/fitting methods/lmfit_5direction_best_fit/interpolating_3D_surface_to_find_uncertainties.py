import numpy as np
import pandas as pd
import astropy.constants as const
import matplotlib.pyplot as plt
import timeit
import scipy.stats as ss
from scipy.signal import argrelextrema
from matplotlib import cm
from scipy import interpolate
import scipy.ndimage.filters as filters
import scipy.ndimage.morphology as morphology
from scipy.interpolate import griddata

def detect_local_minima(arr):
    # https://stackoverflow.com/questions/3684484/peak-detection-in-a-2d-array/3689710#3689710
    """
    Takes an array and detects the troughs using the local maximum filter.
    Returns a boolean mask of the troughs (i.e. 1 when
    the pixel's value is the neighborhood maximum, 0 otherwise)
    """
    # define an connected neighborhood
    # http://www.scipy.org/doc/api_docs/SciPy.ndimage.morphology.html#generate_binary_structure
    neighborhood = morphology.generate_binary_structure(len(arr.shape),2)
    # apply the local minimum filter; all locations of minimum value 
    # in their neighborhood are set to 1
    # http://www.scipy.org/doc/api_docs/SciPy.ndimage.filters.html#minimum_filter
    local_min = (filters.minimum_filter(arr, footprint=neighborhood)==arr)
    # local_min is a mask that contains the peaks we are 
    # looking for, but also the background.
    # In order to isolate the peaks we must remove the background from the mask.
    # 
    # we create the mask of the background
    background = (arr==0)
    # 
    # a little technicality: we must erode the background in order to 
    # successfully subtract it from local_min, otherwise a line will 
    # appear along the background border (artifact of the local minimum filter)
    # http://www.scipy.org/doc/api_docs/SciPy.ndimage.morphology.html#binary_erosion
    eroded_background = morphology.binary_erosion(
        background, structure=neighborhood, border_value=1)
    # 
    # we obtain the final mask, containing only peaks, 
    # by removing the background from the local_min mask
    detected_minima = local_min ^ eroded_background
    return np.where(detected_minima)  


Bmin = 0.0028
Bmax = 0.005
stepsize_B = 0.0001
Bs  = np.arange(Bmin, Bmax, stepsize_B)
print(Bs)
print(len(Bs))

Tmin = 40
Tmax = 70
stepsize_T= 1
Ts =  np.arange(Tmin, Tmax, stepsize_T)
print(Ts)
print(len(Ts))

BB, TT = np.meshgrid(Bs, Ts)
points = np.empty((660, 2))
points[:, 0] = BB.flatten()
points[:, 1] = TT.flatten()


BB = pd.read_csv(r'//Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/fitting methods/BB_Bmin_0.0028_Bmax_0.0034900000000000018_stepsize_1e-05_.txt', delim_whitespace=(True), header = None)
TT = pd.read_csv(r'/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/fitting methods/Bmax_0.05_TT_Tmin_60.0_Tmax_67.5_stepsize_0.5_.txt', delim_whitespace=(True), header = None)
red_chi = pd.read_csv(r'/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/fitting methods/red_chi_finer_B_0.0028_to_0.0035.txt', delim_whitespace=(True), header = None)

#correct origin 166937
BB = pd.read_csv(r'/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/fitting methods/166937_BB_Bmin_0.0028_Bmax_0.004899999999999996_stepsize_0.0001_.txt', delim_whitespace=(True), header = None)
TT = pd.read_csv(r'/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/fitting methods/166937_Bmax_0.005_TT_Tmin_40_Tmax_69_stepsize_1_.txt', delim_whitespace=(True), header = None)
red_chi = pd.read_csv(r'/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/fitting methods/166937_red_chi_Bmin_0.0028_Bmax_0.004899999999999996_Tmin_40_Tmax_69_.txt', delim_whitespace=(True), header = None)


print(BB.shape)
print(TT.shape)
print(red_chi.shape)

ax = plt.axes(projection = '3d')
#ax.plot_surface(BB, TT, red_chi)

red_chi = np.array(red_chi)

f = interpolate.interp2d(BB, TT, red_chi, kind='quintic')


#import matplotlib.pyplot as plt
# stepsize_B = 0.000001
# BBnew = np.arange(0.0028, 0.0035, stepsize_B)
# stepsize_T= 0.1
# TTnew =  np.arange(60, 68, stepsize_T)
# BBM, TTM = np.meshgrid(BBnew, TTnew)

#red_chinew = f(BBnew, TTnew)

Bmin = 0.0028
Bmax = 0.005
stepsize_B = 0.000001
Bs  = np.arange(Bmin, Bmax, stepsize_B)
print(Bs)
print(len(Bs))

Tmin = 40
Tmax = 70
stepsize_T= 0.1
Ts =  np.arange(Tmin, Tmax, stepsize_T)
print(Ts)
print(len(Ts))
BBM, TTM = np.meshgrid(Bs, Ts)
print(BBM.shape)
print(TTM.shape)
#print(red_chinew.shape)

red_chinew = griddata(points, red_chi.flatten(), (BBM, TTM), method='linear')

#ax.plot_surface(BBM, TTM, red_chinew)





value = 197.7
chi = 179 * red_chinew 
ax.plot_surface(BBM, TTM, chi, rstride=1, cstride=1, alpha=0.2, linewidth=0, cmap='summer', zorder = 0.5) #antialiased=True)

local_minima_locations = detect_local_minima(red_chinew)
print(local_minima_locations)
lml_i = local_minima_locations[0]
lml_j= local_minima_locations[1]

red_chi_mins = []

#ax.scatter3D(BBM[25][249], TTM[25][249], chi[25][249], marker = 'o', c = 'black')

# print('------------')
for i,j in zip(lml_i, lml_j):
    ax.scatter3D(BBM[i][j], TTM[i][j], chi[i][j], marker = 'o', c = 'black')
    print(BBM[i][j])
    #print(BB[i][j])
    print(TTM[i][j])
    print(red_chinew[i][j])
    red_chi_mins.append(red_chinew[i][j])
    print(chi[i][j])
    print('----///////--------')


chi = 179 * red_chinew 
indices = []   
lower_value = 197.6
upper_value = 197.7

for i in range(len(chi)):
        for j in range(len(chi[i])):
            if lower_value < chi[i][j] < upper_value:
                indices.append((i, j))
                #print(chi[i][j])

#print(indices)

# mask = chi < 214.26
# indices = np.where(mask)
# print(chi.shape)
# BBM_f = BBM[indices]
# TTM_f = TTM[indices]
# chi_f = chi[indices]
# print(chi_f.shape)

chi[chi >= 197.7] = None
chi_local = chi
#print(chi)
#print(max(chi.shape[0]))



ax.plot_surface(BBM, TTM, chi_local, cmap = 'autumn', alpha = 1, zorder = 1)
indi = np.argwhere(~np.isnan(chi_local))

BBM_local = []
TTM_local = []
for i in indi:
    ii = i[0]
    jj = i[1]
    BBM_local.append(BBM[ii][jj])
    TTM_local.append(TTM[ii][jj])
    
print(min(BBM_local))
print(max(BBM_local))
print(min(TTM_local))
print(max(TTM_local))
    
    
# B_minus = (min(BBM_local) - BBM[25][249])
# B_plus = (max(BBM_local) - BBM[25][249])

# print(B_minus)
# print(B_plus)
# print(min(BBM_local))
# print(BBM[25][249])
# print(max(BBM_local))


# T_minus = (min(TTM_local) - TTM[25][249])
# T_plus = (max(TTM_local) - TTM[25][249])

# print(T_minus)
# print(T_plus)
# print(min(TTM_local))
# print(TTM[25][249])
# print(max(TTM_local))



# ax.view_init(20, 120)

# ax.set_xlim(0.0029,0.0032)
# ax.set_ylim(60, 66.4)


# ax.set_xlabel('B', size = 20, labelpad = 20)
# ax.set_ylabel('T', size = 20, labelpad = 20)
# ax.set_zlabel(r'$\chi^{2}$',  size = 20, labelpad = 20, rotation   = 90)


