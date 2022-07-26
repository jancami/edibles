# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.8
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

import numpy as np
import matplotlib.pyplot as plt
import copy


def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))


data = np.loadtxt(r'C:\Users\hkhan\Downloads\coronene.txt', skiprows = 1)
#print(datas)

# +
waves = (1/copy.deepcopy(data[:, 0]))*1e8
#print(np.min(waves))
#print(np.max(waves))

x = np.linspace(6613.5, 6614.1, 10000)
y = np.zeros(x.shape)

amp = copy.deepcopy(data[:, 1])

fig1, ax1 = plt.subplots(figsize = (12,5))

for it1 in range(len(waves)):
    ax1.plot([waves[it1], waves[it1]], [0, amp[it1]], color = 'red')

ax1.set(xlabel = r'Wavelength ($\AA$)', ylabel = 'Normalised intensities')
fig1.show()
# -

"""
waves = np.linspace(0.5, 1, 400)

x = np.linspace(0.5, 1, 10000)

y = np.zeros(x.shape)

amp = np.ones(waves.shape)
#print(waves)
for it1 in range(len(waves)):
    amp[it1] = np.abs(np.random.normal(1, 10, 1))
    #plt.plot(waves[it1], amp[0])
    #y = y + amp*gaussian(x, waves[it1], 0.095/2.355)

plt.plot(waves, amp)
"""

for it2 in range(len(waves)):
    y = y + amp[it2]*gaussian(x, waves[it2], waves[it2]/(2.355*110000))

# +
# #%matplotlib inline

fig2, ax2 = plt.subplots(figsize = (12,5))

ax2.plot(x, y/np.max(y))
ax2.set(xlabel = r'Wavelength ($\AA$)', ylabel = 'Normalised intensities')
fig2.show()
#plt.set
# -


