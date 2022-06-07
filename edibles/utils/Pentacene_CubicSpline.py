#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.constants as cst
import astropy.units as u
from heapq import nsmallest
from scipy.signal import find_peaks
from functions import vac2air_ciddor
from scipy.interpolate import interp1d
from scipy.interpolate import CubicSpline


# In[2]:


snr5000 = np.loadtxt("Pentacene_air_snr5000.txt").transpose()
snr2000 = np.loadtxt("Pentacene_air_snr2000.txt").transpose()
snr1000 = np.loadtxt("Pentacene_air_snr1000.txt").transpose()
snr500 = np.loadtxt("Pentacene_air_snr500.txt").transpose()
snr100 = np.loadtxt("Pentacene_air_snr100.txt").transpose()


# In[3]:


abs_peak = np.array([5267.577214, 5288.886049, 5305.721707, 5337.921532, 5361.114095])
peak_start = np.array([5261.009267, 5284.022748, 5303.017072, 5334.004866, 5358.019634])
peak_end = np.array([5274.000218, 5291.019457, 5310.017477, 5344.020308, 5366.001094])


# In[4]:


def spectral_stacker(snr_val):
    
    snr_val = str(snr_val)
    varName = "snr"+snr_val
    fdata=globals()[varName]
    
    fpeak1 = fdata[:, np.logical_and(fdata[0]>=peak_start[0], fdata[0]<=peak_end[0])]
    fpeak2 = fdata[:, np.logical_and(fdata[0]>=peak_start[1], fdata[0]<=peak_end[1])]
    fpeak3 = fdata[:, np.logical_and(fdata[0]>=peak_start[2], fdata[0]<=peak_end[2])]
    fpeak4 = fdata[:, np.logical_and(fdata[0]>=peak_start[3], fdata[0]<=peak_end[3])]
    fpeak5 = fdata[:, np.logical_and(fdata[0]>=peak_start[4], fdata[0]<=peak_end[4])]

    fD1 = np.zeros(fpeak1.shape)
    fD1[0] = fpeak1[0]-abs_peak[0]
    fD1[1] = fpeak1[1]
    
    fD2 = np.zeros(fpeak2.shape)
    fD2[0] = fpeak2[0]-abs_peak[1]
    fD2[1] = fpeak2[1]

    fD3 = np.zeros(fpeak3.shape)
    fD3[0] = fpeak3[0]-abs_peak[2]
    fD3[1] = fpeak3[1]
    
    fD4 = np.zeros(fpeak4.shape)
    fD4[0] = fpeak4[0]-abs_peak[3]
    fD4[1] = fpeak4[1]

    fD5 = np.zeros(fpeak5.shape)
    fD5[0] = fpeak5[0]-abs_peak[4]
    fD5[1] = fpeak5[1]
    
    start = (peak_start-abs_peak).max()
    end = (peak_end-abs_peak).min()
    fpoints = max(fD1.shape[1], fD2.shape[1], fD3.shape[1], fD4.shape[1], fD5.shape[1])
    fwavelength = np.linspace(start, end, num=fpoints)

    
    fp1 = CubicSpline(fD1[0], fD1[1])
    fp2 = CubicSpline(fD2[0], fD2[1])
    fp3 = CubicSpline(fD3[0], fD3[1])
    fp4 = CubicSpline(fD4[0], fD4[1])
    fp5 = CubicSpline(fD5[0], fD5[1])

    ffig, faxs = plt.subplots(1, 2, figsize=(18,6))
    plt.tight_layout(rect=[0,0.03,1,0.95])
    
    fp1in = np.zeros((2,fpoints))
    fp1in[0] = fwavelength
    fp1in[1] = fp1(fwavelength)
    faxs[0].plot(fp1in[0],fp1in[1], label='Peak1')
    
    fp2in = np.zeros((2,fpoints))
    fp2in[0] = fwavelength
    fp2in[1] = fp2(fwavelength)
    faxs[0].plot(fp2in[0],fp2in[1], label='Peak2')
    

    fp3in = np.zeros((2,fpoints))
    fp3in[0] = fwavelength
    fp3in[1] = fp3(fwavelength)
    faxs[0].plot(fp3in[0],fp3in[1], label='Peak3')
    
    
    fp4in = np.zeros((2,fpoints))
    fp4in[0] = fwavelength
    fp4in[1] = fp4(fwavelength)
    faxs[0].plot(fp4in[0],fp4in[1], label='Peak4')
    
    
    fp5in = np.zeros((2,fpoints))
    fp5in[0] = fwavelength
    fp5in[1] = fp5(fwavelength)
    faxs[0].plot(fp5in[0],fp5in[1], label='Peak5')
    faxs[0].set_xlabel("Relative Wavelength (Å)", size=16)
    faxs[0].set_ylabel("Relative Intensity", size=16)
    faxs[0].grid()
    faxs[1].grid()
    faxs[0].legend()
    

    fresult = np.zeros((2,fpoints))
    fresult[0] = fwavelength
    fresult[1] = (fp1in[1]+fp2in[1]+fp3in[1]+fp4in[1]+fp5in[1])/5.0
    
    
    faxs[1].plot(fresult[0],fresult[1], label="Final Stacked Image")
    
    plt.xlabel("Wavelength (Å)",size=16)
    plt.ylabel("Relative Intensity", size=16)
    plt.subplots_adjust(wspace=0.2)
    title="Stacked Image of SNR"+snr_val
    plt.suptitle(title, size=22)
    plt.legend()
    
    


# In[5]:


spectral_stacker(5000)


# In[6]:


spectral_stacker(2000)


# In[7]:


spectral_stacker(1000)


# In[8]:


spectral_stacker(500)


# In[9]:


spectral_stacker(100)


# In[ ]:




