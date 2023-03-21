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


stepsize_B = 0.00001
Bs  = np.arange(0.0028, 0.0032, stepsize_B)
#Bs  = np.arange(0.0005, 0.01, stepsize_B)
#Bs  = np.arange(1, 10, 4)
print(Bs)
print(len(Bs))

stepsize_T= 0.5
Ts =  np.arange(57, 68, stepsize_T)
#Ts =  np.arange(5, 100, 40)
print(Ts)
print(len(Ts))
BB, TT = np.meshgrid(Bs, Ts)
points = np.empty((902, 2)) #902
points[:, 0] = BB.flatten()
points[:, 1] = TT.flatten()

#1666937 (broadest)
BB = pd.read_csv(r'/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/fitting methods/BB_Bmin_0.0028_Bmax_0.003200000000000001_stepsize_1e-05_.txt', delim_whitespace=(True), header = None)
TT = pd.read_csv(r'/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/fitting methods/Bmax_0.0032_TT_Tmin_57.0_Tmax_67.5_stepsize_0.5_.txt', delim_whitespace=(True), header = None)
red_chi = pd.read_csv(r'/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/fitting methods/red_chi_finer_B_0.0028_to_0.0032_57_to_68.txt', delim_whitespace=(True), header = None)

#185418 (narrowest)
BB = pd.read_csv(r'/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/fitting methods/185418_BB_Bmin_0.0028_Bmax_0.003200000000000001_stepsize_1e-05_.txt', delim_whitespace=(True), header = None)
TT = pd.read_csv(r'/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/fitting methods/185418_Bmax_0.0032_TT_Tmin_57.0_Tmax_67.5_stepsize_0.5_.txt', delim_whitespace=(True), header = None)
red_chi = pd.read_csv(r'/Users/charmibhatt/Desktop/Local_GitHub/edibles/edibles/utils/simulations/Charmi/fitting methods/185418_red_chi_finer_B_0.0028_to_0.0032_57_to_68.txt', delim_whitespace=(True), header = None)


    
print(BB.shape)
print(TT.shape)
print(red_chi.shape)

ax = plt.axes(projection = '3d', computed_zorder=False)
#ax.plot_surface(BB, TT, red_chi)

red_chi = np.array(red_chi)

#f = interpolate.interp2d(BB, TT, red_chi, kind='linear')


#import matplotlib.pyplot as plt
stepsize_B = 0.000001
BBnew = np.arange(0.0028, 0.0032, stepsize_B)
stepsize_T= 0.05
TTnew =  np.arange(57, 68, stepsize_T)
BBM, TTM = np.meshgrid(BBnew, TTnew)

#red_chinew = f(BBnew, TTnew)

print(BBM.shape)
print(TTM.shape)
#print(red_chinew.shape)

red_chinew = griddata(points, red_chi.flatten(), (BBM, TTM), method='linear')

#ax.plot_surface(BBM, TTM, red_chinew)



chi = 179 * red_chinew 
ax.plot_surface(BBM, TTM, chi, rstride=1, cstride=1, alpha=0.5, linewidth=0, cmap='summer', zorder = 1) #antialiased=True)


local_minima_locations = detect_local_minima(red_chinew)
print(local_minima_locations)
lml_i = local_minima_locations[0]
lml_j= local_minima_locations[1]

#red_chi_mins = []


# print('------------')
for i,j in zip(lml_i, lml_j):
    #ax.scatter3D(BBM[i][j], TTM[i][j], chi[i][j], marker = 'o', c = 'black')
    print(BBM[i][j])
    #print(BB[i][j])
    print(TTM[i][j])
    print(red_chinew[i][j])
    #red_chi_mins.append(red_chinew[i][j])
    print(chi[i][j])
    print(i)
    print(j)
    print('----///////--------')


chi = 179 * red_chinew 
indices = []   
lower_value = 247.25 #214.2674
upper_value = 247.255 #214.2676

for i in range(len(chi)):
        for j in range(len(chi[i])):
            if lower_value < chi[i][j] < upper_value:
                indices.append((i, j))
                print(chi[i][j])

print(indices)

# mask = chi < 214.26
# indices = np.where(mask)
# print(chi.shape)
# BBM_f = BBM[indices]
# TTM_f = TTM[indices]
# chi_f = chi[indices]
# print(chi_f.shape)

#chi[chi >= 214.26] = None
chi[chi >= 247.25] = None
chi_local = chi
print(chi)
#print(max(chi.shape[0]))



ax.plot_surface(BBM, TTM, chi_local, cmap = 'pink_r', alpha = 1, zorder = 2)
indi = np.argwhere(~np.isnan(chi_local))

BBM_local = []
TTM_local = []
for i in indi:
    ii = i[0]
    jj = i[1]
    BBM_local.append(BBM[ii][jj])
    TTM_local.append(TTM[ii][jj])


i = 50 #110
j = 399 #250 
print('----///////--------')
B_minus = (min(BBM_local) - BBM[i][j])
B_plus = (max(BBM_local) - BBM[i][j])

# print(B_minus)
# print(B_plus)
print(min(BBM_local))
print(BBM[i][j])
print(max(BBM_local))


T_minus = (min(TTM_local) - TTM[i][j])
T_plus = (max(TTM_local) - TTM[i][j])

print(T_minus)
print(T_plus)
print(min(TTM_local))
print(TTM[i][j])
print(max(TTM_local))



ax.view_init(20, 120)

#ax.set_xlim(0.00285,0.0032)
# ax.set_ylim(60, 66.4)
ax.scatter3D(BBM[i][j], TTM[i][j], chi[i][j], marker = 'o', c = 'black', zorder = 3, s= 5)


ax.set_xlabel('B', size = 10, labelpad = 30)
ax.set_ylabel('T', size = 10, labelpad = 30)
ax.set_zlabel(r'$\chi^{2}$',  size = 10, labelpad = 10, rotation   = 90)
ax.tick_params(axis='x', which='major', labelsize=10, rotation = 80)
ax.tick_params(axis='y', which='major', labelsize=8, rotation = 0)

ax.tick_params(axis='both', which='minor', labelsize=8)

