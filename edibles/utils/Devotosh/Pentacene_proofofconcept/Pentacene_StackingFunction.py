import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.constants as cst
import astropy.units as u
from heapq import nsmallest
from scipy.signal import find_peaks
from functions import vac2air_ciddor
from scipy.interpolate import interp1d



# Loading the files which consists of the lab simulated spectra of Pentacene. (Here, we are importing five files, 
#since each has different SNR values)
snr5000 = np.loadtxt("Pentacene_air_snr5000.txt")
snr2000 = np.loadtxt("Pentacene_air_snr2000.txt")
snr1000 = np.loadtxt("Pentacene_air_snr1000.txt")
snr500 = np.loadtxt("Pentacene_air_snr500.txt")
snr100 = np.loadtxt("Pentacene_air_snr100.txt")


# Creating arrays manually for the absorptions.
peaks = np.zeros((5,3))
peaks[0,0] = 5261.009267           #starting value of abs
peaks[0,1] = 5267.577214           #peak values of abs
peaks[0,2] = 5274.000218           #ending values of abs

peaks[1,0] = 5284.022748
peaks[1,1] = 5288.886049
peaks[1,2] = 5291.019457

peaks[2,0] = 5303.017072
peaks[2,1] = 5305.721707
peaks[2,2] = 5310.017477

peaks[3,0] = 5334.004866
peaks[3,1] = 5337.921532
peaks[3,2] = 5344.020308

peaks[4,0] = 5358.019634
peaks[4,1] = 5361.114095
peaks[4,2] = 5366.001094


#Function for Spectral stacking, which takes in the input as an interger which is the SNR value of the file with the 
# spectral data

def spectral_stacker(fdata,peaks):
    
    No_Peaks = peaks.shape[0]
    fpeaks = np.empty(shape=No_Peaks, dtype=object)
    
    ffig, faxs = plt.subplots(3, 2, figsize=(14,10))
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    
    faxs[0, 0].plot(fdata[:, 0], fdata[:, 1])
    faxs[0, 0].set_title('Lab Spectrum')
    faxs[0, 0].set(xlabel = 'Wavelength (Å)', ylabel = 'Relative intensity')
    
    
    # Creating separate arrays for each absorption in the spectra
    for i1 in range(No_Peaks):
        fpeaks[i1] = fdata[np.logical_and(fdata[:,0]>=peaks[i1,0], fdata[:,0]<=peaks[i1,2]), :]
        p = 'Peak'+str(i1+1)
        faxs[0,1].plot(fpeaks[i1][:,0], fpeaks[i1][:,1], label=p)
    faxs[0, 1].set_title('Peaks')
    faxs[0, 1].set(xlabel = 'Wavelength (Å)', ylabel = 'Relative intensity')
    faxs[0, 1].legend()
    
    fDshft = fpeaks
    for j1 in range (No_Peaks):
        fDshft[j1][:,0] = fpeaks[j1][:,0]-peaks[j1,1]
        fDshft[j1][:,1] = fpeaks[j1][:,1]
        shp = 'Peak ' + str(j1+1)
        faxs[1, 0].plot(fDshft[j1][:, 0], fDshft[j1][:, 1], label = shp)
        
    faxs[1, 0].set_title('Spectrum with peaks shifted')
    faxs[1, 0].set(xlabel = 'Relative Wavelength (Å) (shifted by peak)', ylabel = 'Relative intensity')
    faxs[1, 0].legend()
    
    fpoints = fDshft[0].shape[0]
    fstart = peaks[0][0]-peaks[0][1]
    fend = peaks[0][2]-peaks[0][1]
    for k1 in range(No_Peaks):
        if fDshft[k1].shape[0]>=fpoints:
            fpoints = fDshft[k1].shape[0]
        if peaks[k1][0]-peaks[k1][1] > fstart:
            fstart = peaks[k1][0]-peaks[k1][1]
        if peaks[k1][2]-peaks[k1][1] < fend:
            fend = peaks[k1][2]-peaks[k1][1]
    fwavelength = np.linspace(fstart, fend, num = fpoints)  
    
            
    fitpn = np.empty(shape=No_Peaks, dtype=object)  
    for l1 in range(No_Peaks):
        fitpn[l1] = interp1d(fDshft[l1][:,0], fDshft[l1][:,1])
        
        
    ffin = np.zeros((fpoints,2))
    ffin[:,0] = fwavelength
    for m1 in range(No_Peaks):
        ffin[:,1] = ffin[:,1]+(fitpn[m1](fwavelength)/No_Peaks)
        pk1 = 'Peak' +str(m1+1)
        faxs[1,1].plot(fwavelength, fitpn[m1](fwavelength), label=pk1)         
        faxs[1,1].set_title('Peaks with same widths')
        faxs[1,1].set(xlabel = 'Relative Wavelength (Å) (shifted by peak)', ylabel = 'Relative intensity')
        faxs[1,1].legend()
        faxs[2,0].plot(fwavelength, fitpn[m1](fwavelength),label=pk1)
        faxs[2,0].set_title('Individual Peaks with the stacked peak')
        faxs[2,0].set(xlabel = 'Relative Wavelength (Å) (shifted by peak)', ylabel = 'Relative intensity')
        faxs[2,0].legend()
        
    faxs[2, 0].plot(ffin[:, 0], ffin[:, 1], label='Stacked spectra')
    faxs[2,0].legend()
    faxs[2, 1].plot(ffin[:, 0], ffin[:, 1], label = 'Final Stacked peak')
    faxs[2, 1].set_title('Only Stacked peak')
    faxs[2, 1].set(xlabel = 'Relative Wavelength (Å)(shifted by peak)', ylabel = 'Relative intensity')
    faxs[2, 1].legend()
    plt.subplots_adjust(hspace=0.5, wspace=0.4)
    
    
    
    
  spectral_stacker(snr5000,peaks)
  spectral_stacker(snr2000,peaks)
  spectral_stacker(snr1000,peaks)
  spectral_stacker(snr500,peaks)
  spectral_stacker(snr100,peaks)
