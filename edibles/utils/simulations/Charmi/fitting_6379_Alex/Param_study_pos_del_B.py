# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 15:49:04 2023

@author: alexr
"""

from functions import allowed_parallel_transitions,get_rotational_spectrum
    
import matplotlib.pyplot as plt

import beepy as bp

#%%
Jmax = 300

combinations = allowed_parallel_transitions(Jmax)
# print(combinations)

delta_B = 0.8
delta_C = delta_B
zeta = -0.55
sigma = 0.1953
origin = 0

#Ranges of temperatures and Bs
Ts = (2.7, 10, 30, 70, 100)  
Bs = (0.000501, 0.001584, 0.005011, 0.0158489, 0.0501187)
B_label = ('5.0 x 10$^{-4}$', '1.6 x 10$^{-3}$', '5.0 x 10$^{-3}$', '1.6 x 10$^{-2}$', '5.0 x 10$^{-2}$' )

#Create figure
fig, axes = plt.subplots(5, 5, figsize=(19,10), sharex=(True), sharey=(True))
fig.suptitle(' 'r'$\Delta B =$ ' + str(delta_B) + '% ,  'r'$\Delta C =$ ' + str(delta_C) + '% , 'r'$\zeta^{\prime}  = $' + str(zeta) + ', 'r'$\sigma = $'+ str(sigma) + 'cm$^{-1}$ \n' , size ='xx-large')

#Labelling T and B values for rows and columns
rows = ['T = {} K'.format(row) for row in Ts ]
cols = ['B = {} cm$^{}$ '.format(col, {-1}) for col in B_label]

#Labelling figure
for ax, col in zip(axes[0], cols):
    ax.set_title(col, fontsize = 15)
    
for ax, col in zip(axes[4], cols):
    ax.set_xlabel('Wavenumber (cm$^{-1}$)', labelpad =10, fontsize = 13)

for ax, row in zip(axes[:,0], rows):
    ax.set_ylabel('Intensity', rotation=90, labelpad=7, fontsize = 13)

pad = 25 # in points    
for ax, row in zip(axes[:,0], rows):
    ax.annotate(row, xy=(0, 0.5), xytext=(-ax.yaxis.labelpad - pad, 0),
                xycoords=ax.yaxis.label, textcoords='offset points', ha='right', va='center', fontsize=15)
   
fig.tight_layout()


n = 0
for B in Bs:
    m = 0
    for T in Ts:
        xs, ys = get_rotational_spectrum(B, delta_B, zeta, T, sigma, origin, combinations, Jmax = Jmax, bell = False, transition = 'parallel')
        
        axes[m,n].plot(xs, ys, linewidth = 1)
        axes[m,n].xaxis.set_major_locator(plt.MultipleLocator(5))
        axes[m,n].xaxis.set_minor_locator(plt.MultipleLocator(1))
        axes[m,n].set_xlim(-12,12)
        axes[m,n].axhline(y=1, linestyle = '--', color = 'gray')    
        axes[m,n].yaxis.set_major_locator(plt.MultipleLocator(0.05))
        axes[m,n].yaxis.set_minor_locator(plt.MultipleLocator(0.01))
        
        m = m + 1
        
    n = n + 1
    
# winsound.Beep(440,500)
plt.savefig('Parallel_transition_param_study/del_B={}_zeta={}.png'.format(delta_B,zeta), bbox_inches = 'tight')

bp.beep(sound='ready')

plt.show()
