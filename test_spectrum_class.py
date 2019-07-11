from edibles_settings import *
from edibles_spectrum import *

sp = edibles_spectrum(datadir+"/HD170740/RED_860/HD170740_w860_n20_20140916_L.fits")
print("Barycentric Velocity is", sp.v_bary)
wave,flux = sp.GetSpectrum()
plt.plot(wave, flux)
axes = plt.gca()
axes.set_xlim([7660,7705])
axes.set_ylim([0,160])
plt.vlines((7667.021,7701.093), 0, 160, linestyles='dashed', colors='r')
plt.show()