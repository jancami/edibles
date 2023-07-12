import numpy as np
import math
import matplotlib.pyplot as plt

from edibles.utils.edibles_oracle import EdiblesOracle
from edibles.utils.edibles_spectrum import EdiblesSpectrum
from edibles.utils.voigt_profile import voigt_absorption_line

pythia = EdiblesOracle()
objectlist=["HD 63804", "HD 147084","HD 147933", "HD 147889","HD 148184","HD 154368","HD 161056","HD 169454","HD 185859","HD 186745","HD 186841"]
List = pythia.getFilteredObsList(Wave=6700, OrdersOnly=True, object=["HD 147889"])

print("Number of objects = ", len(List))
    
    
#changing from pandas series to python list as index of observations in pandas is incorrect
files = []
for filename in List:
      files.append(filename)

"""
targetlist=[]
filenamelist=[]

for filename in files:
    sp = EdiblesSpectrum('/' + filename)
    sp.getSpectrum(xmin=6707, xmax=6710)
    targetlist.append(sp.target)
    filenamelist.append(filename)

print(targetlist)
"""

sp = EdiblesSpectrum(filename)
sp.getSpectrum(xmin=6707, xmax=6710)
sp.flux = sp.flux / np.median(sp.flux)
plt.plot(sp.wave, sp.flux)

plt.title(filename)
plt.xlabel(r'Wavelength ($\AA$)')
plt.ylabel('Flux')


############################################################################################################################################


lambda0 = [6707.761, 6707.912]
f = [8.5e-03,5.2e-03]   
gamma = [5e8,5e8]
b = [0.4]
N = [1e12]
v_rad = 19                  #check samuel's report for this parameter
v_resolution = 4


wrange = [6705, 6710]
wave = sp.wave
flux = sp.flux
idx = np.where((wave > wrange[0]) & (wave < wrange[1]))
wave = wave[idx]
flux = flux[idx]
flux = flux / np.median(flux)
       
        
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


############################################################################################################################


x1=6708.95
x2=6709.39

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
Md=AbsorptionLine[idx]   #model data - but this is the straight line part so how does that work?
#print(Od) 
#print(Md)    

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


