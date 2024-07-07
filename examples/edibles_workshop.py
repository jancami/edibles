import matplotlib.pyplot as plt

from edibles.utils.edibles_spectrum import EdiblesSpectrum
from edibles.sightline import Sightline


# # #########################################################################
# # Find files

# %%%%%%%
# GUI!!!!
# %%%%%%%


# #########################################################################
# Create a spectrum from a fits file

file = "/HD170740/RED_860/HD170740_w860_redl_20140915_O12.fits"

print()
sp = EdiblesSpectrum(file)

# print(sp.target)
# print(sp.date)
# print(sp.v_bary)
# print(sp.wave_units)
# print(sp.flux_units)
# print(sp.wave)
# print(sp.bary_wave)
# print(sp.flux)

# ##############################
# Plotting

# plt.plot(sp.wave, sp.flux, label='Geocentric')
# plt.plot(sp.bary_wave, sp.flux, label='Barycentric')
# plt.title(sp.target + ", " + sp.date[:10])
# plt.legend()
# plt.show()



# print()
# print("#" * 80)
# print()

# ##############################

sp.getSpectrum(7657, 7670)


plt.plot(sp.wave, sp.flux, label='Geocentric')
# plt.plot(sp.bary_wave, sp.flux, label='Barycentric')
plt.title(sp.target + ", " + sp.date[:10])
plt.legend()
plt.show()

# #########################################################################


# Create sightline model (with continuum)
sightline = Sightline(sp)


# Fit the sightline model
# sightline.fit(plot=True)


# Add source of absorption
sightline.add_source('Telluric')


# Add line with all guessed parameters
sightline.add_line(name='line1', source='Telluric')
# sightline.fit(plot=True)


# # # # Add lines with given starting parameters
pars = {'d': 0.01, 'tau_0': 0.3, 'lam_0': 7665.9}
sightline.add_line(name='line2', pars=pars, source='Telluric')
pars = {'d': 0.01, 'tau_0': 0.3, 'lam_0': 7664.8}
sightline.add_line(name='line3', pars=pars, source='Telluric')
pars = {'d': 0.01, 'tau_0': 0.3, 'lam_0': 7659.2}
sightline.add_line(name='line4', pars=pars, source='Telluric')
sightline.fit(plot=True)


# evaluate the model with a set of input parameters (result of fit here)
out = sightline.complete_model.eval(data=sp.flux, params=sightline.result.params, x=sp.wave)
resid = sp.flux - out


# Plotting evaluation and residual
plt.plot(sp.wave, sp.flux)
plt.plot(sp.wave, out)
plt.plot(sp.wave, resid)
plt.show()


# Add line with guessed parameters from residual data
sightline.add_source('Interstellar')
sightline.add_line(name='K_line', source='Interstellar', guess_data=resid)

sightline.fit(report=True, plot=True)
