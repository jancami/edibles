import numpy as np
import math
import matplotlib.pyplot as plt

from edibles.utils.edibles_oracle import EdiblesOracle
from edibles.utils.edibles_spectrum import EdiblesSpectrum
from edibles.utils.voigt_profile import voigt_absorption_line

objectlist=["HD 63804", "HD 147084","HD 147933", "HD 147889","HD 148184","HD 154368","HD 161056","HD 169454","HD 185859","HD 186745","HD 186841"]
filename = "HD 147889"
  	 
pythia = EdiblesOracle()
List_series = pythia.getFilteredObsList(Wave=6700, OrdersOnly=True, object=[filename])
List = List_series.tolist()

print("Number of objects = ", len(List))

sp = EdiblesSpectrum(List[0])
sp.getSpectrum(xmin=6706.5, xmax=6709.5)
sp.flux = sp.flux / np.median(sp.flux)
plt.plot(sp.wave, sp.flux)

plt.title(filename)
plt.xlabel(r'Wavelength ($\AA$)')
plt.ylabel('Flux')


############################################################################################################################################
#manual fitting

lambda0 = [6707.761, 6707.912]
f = [8.6e-03,5.2e-03]   
gamma = [5e8,5e8]
b = [0.4]
N = [1e12]
v_rad = 19            #check samuel's report for this parameter
v_resolution = 4


wave = sp.wave
flux = sp.flux
        
AbsorptionLine = voigt_absorption_line(
            wave,
            lambda0=lambda0,
            b=b,
            N=N,
            f=f,
            gamma=gamma,
            v_rad=v_rad,
            v_resolution=v_resolution)

plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
plt.plot(wave, AbsorptionLine, color="orange", marker="*")
plt.xlim(6707.5,6709) 


############################################################################################################################
#calculating chi square

x1=6707.86
x2=6708.52

plt.axvline(x=x1, color = 'g', linestyle = '-')
plt.axvline(x=x2, color = 'g', linestyle = '-')

idx = np.where((wave > x1) & (wave < x2))
xvalues = wave[idx]
yvalues = flux[idx]
line_yvalue=np.average(yvalues) #ie average of the yvalues

plt.axhline(y=line_yvalue, color = 'r', linestyle = '-')


sigma=abs(yvalues-line_yvalue)
print("sigma = ", sigma)

avsigma=np.average(sigma)

print("avsigma = ", avsigma)
m=line_yvalue+avsigma
n=line_yvalue-avsigma
plt.axhline(y=m, color = 'g', linestyle = '-')
plt.axhline(y=n, color = 'g', linestyle = '-') 

Od=flux[idx]             #observed data
Md=AbsorptionLine[idx]   #model data 
print(Od) 
print(Md)    

subtracted_array = np.subtract(Od,Md)
subtracted = np.square(subtracted_array)
subtracted = list(subtracted)

chi_sqrd=sum(subtracted/sigma)

print("chi_sqrd = ",chi_sqrd)

n=len(Od)    #number of observations
m=8          #number of fitted parameters
v=n-m
reduced_chi_sqrd=chi_sqrd/v
print("reduced_chi_sqrd = ", reduced_chi_sqrd)


plt.show()


