import matplotlib.pyplot as plt

from edibles.utils.edibles_spectrum import EdiblesSpectrum
from edibles.models import ContinuumModel, VoigtModel


filename = "/HD170740/RED_860/HD170740_w860_redl_20140915_O12.fits"

method = 'least_squares'

sp = EdiblesSpectrum(filename)
print(sp.target)
subset = sp.getSpectrum(xmin=7661, xmax=7670)

# #################################################################################

cont_model = ContinuumModel(n_anchors=4)
cont_pars = cont_model.guess(subset.flux, x=subset.wave)
model = cont_model
pars = cont_pars

result = model.fit(data=subset.flux, params=pars, x=subset.wave, method=method)
out = cont_model.eval(data=subset.flux, params=result.params, x=subset.wave)
resid = subset.flux - out

# result.plot_fit()
# plt.show()

# #################################################################################

voigt1 = VoigtModel(prefix='voigt1_')
voigt1_pars = voigt1.guess(resid, x=subset.wave)


model = model * voigt1
pars = result.params + voigt1_pars

result = model.fit(data=subset.flux, params=pars, x=subset.wave, method=method)
out = model.eval(data=subset.flux, params=result.params, x=subset.wave)
resid = subset.flux - out

# result.plot_fit()
# plt.show()

# #################################################################################

voigt2 = VoigtModel(prefix='voigt2_')
voigt2_pars = voigt2.guess(resid, x=subset.wave)


model = model * voigt2
pars = result.params + voigt2_pars

result = model.fit(data=subset.flux, params=pars, x=subset.wave, method=method)
out = model.eval(data=subset.flux, params=result.params, x=subset.wave)
resid = subset.flux - out

# result.plot_fit()
# plt.show()

# #################################################################################

voigt3 = VoigtModel(prefix='voigt3_')
voigt3_pars = voigt3.guess(resid, x=subset.wave)


model = model * voigt3
pars = result.params + voigt3_pars

result = model.fit(data=subset.flux, params=pars, x=subset.wave, method=method)
out = model.eval(data=subset.flux, params=result.params, x=subset.wave)
resid = subset.flux - out

# result.plot_fit()
# plt.show()

# #################################################################################

voigt4 = VoigtModel(prefix='voigt4_')
voigt4_pars = voigt4.guess(resid, x=subset.wave)


model = model * voigt4
pars = result.params + voigt4_pars

result = model.fit(data=subset.flux, params=pars, x=subset.wave, method=method)
out = model.eval(data=subset.flux, params=result.params, x=subset.wave)
resid = subset.flux - out


# #################################################################################

print(result.fit_report())

plt.plot(subset.wave, subset.flux)
plt.plot(subset.wave, out)
plt.plot(subset.wave, resid)
plt.show()
