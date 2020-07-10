import matplotlib.pyplot as plt

from edibles.sightline import Sightline
from edibles.utils.edibles_spectrum import EdiblesSpectrum



FILE1 = "/HD170740/RED_860/HD170740_w860_redl_20140915_O12.fits"
xmin = 7661.5
xmax = 7669

sp1 = EdiblesSpectrum(FILE1)
sp1.getSpectrum(xmin=7661, xmax=7670)

sightline = Sightline(sp1)

# Add source
sightline.add_source('telluric')

# Add line with auto-guessed params
sightline.add_line(name='line1', source='telluric')

# Add line with user defined params
d = {'d': 0.01, 'tau_0': 0.6, 'lam_0': 7664.8}
sightline.add_line(name='line2', pars=d, source='telluric')

# Add line with undefined source
d = {'d': 0.01, 'tau_0': 0.1, 'lam_0': 7665.2}

sightline.add_line(name='line3', source='interstellar')

# Add line with no source & user defined pars
d = {'d': 0.01, 'tau_0': 0.1, 'lam_0': 7662}
sightline.add_line(name='line4', pars=d)

# ###############################################################
# Fit and plot
sightline.fit(report=True, plot=True)

out = sightline.model.eval(data=sp1.flux, params=sightline.result.params, x=sp1.wave)
resid = sp1.flux - out

plt.plot(sp1.wave, sp1.flux)
plt.plot(sp1.wave, out)
plt.plot(sp1.wave, resid)
plt.show()
# ###############################################################

# Add line using guess_pars, and link parameters together
sightline.add_line(name='line5', source='interstellar', guess_data=resid)
sightline.model_pars['interstellar_line5_lam_0'].set(expr='interstellar_line3_lam_0 + 0.091')

# ###############################################################
# Fit and plot
sightline.fit(report=True, plot=True)

out = sightline.model.eval(data=sp1.flux, params=sightline.result.params, x=sp1.wave)
resid = sp1.flux - out

plt.plot(sp1.wave, sp1.flux)
plt.plot(sp1.wave, out)
plt.plot(sp1.wave, resid)
plt.show()
# ###############################################################
