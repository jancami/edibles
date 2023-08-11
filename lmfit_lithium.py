import numpy as np
import math
import matplotlib.pyplot as plt
from lmfit import Model, Parameters, report_fit
from edibles.utils import voigt_profile as vp

from edibles.utils.edibles_oracle import EdiblesOracle
from edibles.utils.edibles_spectrum import EdiblesSpectrum
from edibles.utils.voigt_profile import multi_voigt_absorption_line, fit_multi_voigt_absorptionlines

filename = "HD 157978"
  	 
pythia = EdiblesOracle()
List_series = pythia.getFilteredObsList(Wave=6700, OrdersOnly=True, object=[filename])
List = List_series.tolist()

sp = EdiblesSpectrum(List[0])
sp.getSpectrum(xmin=6706, xmax=6710)
sp.flux = sp.flux / np.median(sp.flux)

wave = sp.wave
flux = sp.flux
wavegrid=wave
ydata=flux

restwave = [6707.761, 6707.912, 6707.921, 6708.072]
f = [0.4982, 0.2491, 0.4982, 0.2491]
gamma = [3.69e7, 3.69e7, 3.69e7, 3.69e7]
b=4
N=1e10
v_rad=6
v_resolution=4
n_step=25
chi_sqrd=[]
perc=[]
print("f = ",f)

percentage_list=np.arange(0,60.1,1)

for percentage in percentage_list:
	percentage=round(percentage,1)
	f_use=[0.4982, 0.2491, 0.4982, 0.2491]
	f_use[2]=f[2]*percentage/100
	f_use[3]=f[3]*percentage/100
	print("f corrected=",f_use)
	perc.append(percentage)

	fitresult = fit_multi_voigt_absorptionlines(wavegrid=wave, ydata=flux, restwave=restwave, f=f_use, gamma=gamma, 
                       b=b, N=N, v_rad=v_rad, v_resolution=v_resolution, n_step=n_step)
	fitresult.params.pretty_print()
	
	chi_sqrd.append(fitresult.chisqr)
	
	plt.plot(wave,flux)
	plt.xlim(6706,6710) 
	plt.title(filename + " - All Iterations")
	plt.plot(wave, fitresult.best_fit, color='b')

plt.show()

#finding & plotting chi square
print("chi_sqrd = ",chi_sqrd)
print("perc = ",perc)
print("min chi_sqrd = ",min(chi_sqrd))

idx=chi_sqrd.index(min(chi_sqrd))
percentage=perc[idx]
print("percentage = ", percentage)

x=percentage_list
plt.plot(x,chi_sqrd,color='green', marker='o', linestyle='dashed')
plt.xlabel('Percentage')
plt.ylabel('Chi Squared')
plt.title("Chi Squared: " + filename)
plt.show()	


#plotting best fit
f_use=[0.4982, 0.2491, 0.4982, 0.2491]
f_use[2]=f[2]*percentage/100
f_use[3]=f[3]*percentage/100
print("f corrected=",f_use)

fitresult = fit_multi_voigt_absorptionlines(wavegrid=wave, ydata=flux, restwave=restwave, f=f_use, gamma=gamma, b=b, N=N, v_rad=v_rad, v_resolution=v_resolution, n_step=n_step)
fitresult.params.pretty_print()
	
plt.plot(wave,flux)
plt.xlim(6706,6710)
plt.title("Best Fit: " + filename)
plt.plot(wave, fitresult.best_fit, color='g')
plt.show()


