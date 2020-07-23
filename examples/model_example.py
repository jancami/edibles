import matplotlib.pyplot as plt

from edibles.utils.edibles_spectrum import EdiblesSpectrum
from edibles.models import ContinuumModel, VoigtModel


filename = "/HD170740/RED_860/HD170740_w860_redl_20140915_O12.fits"

method = 'least_squares'

sp = EdiblesSpectrum(filename)
print(sp.target)
sp.getSpectrum(xmin=7661, xmax=7670)

# #################################################################################

cont_model = ContinuumModel(n_anchors=4)
cont_pars = cont_model.guess(sp.flux, x=sp.wave)
model = cont_model
pars = cont_pars

result = model.fit(data=sp.flux, params=pars, x=sp.wave, method=method)
out = cont_model.eval(data=sp.flux, params=result.params, x=sp.wave)
resid = sp.flux - out

# result.plot_fit()
# plt.show()

# #################################################################################

voigt1 = VoigtModel(prefix='voigt1_')
voigt1_pars = voigt1.guess(resid, x=sp.wave)


model = model * voigt1
pars = result.params + voigt1_pars

result = model.fit(data=sp.flux, params=pars, x=sp.wave, method=method)
out = model.eval(data=sp.flux, params=result.params, x=sp.wave)
resid = sp.flux - out

# result.plot_fit()
# plt.show()

# #################################################################################

voigt2 = VoigtModel(prefix='voigt2_')
voigt2_pars = voigt2.guess(resid, x=sp.wave)


model = model * voigt2
pars = result.params + voigt2_pars

result = model.fit(data=sp.flux, params=pars, x=sp.wave, method=method)
out = model.eval(data=sp.flux, params=result.params, x=sp.wave)
resid = sp.flux - out

# result.plot_fit()
# plt.show()

# #################################################################################

voigt3 = VoigtModel(prefix='voigt3_')
voigt3_pars = voigt3.guess(resid, x=sp.wave)


model = model * voigt3
pars = result.params + voigt3_pars

result = model.fit(data=sp.flux, params=pars, x=sp.wave, method=method)
out = model.eval(data=sp.flux, params=result.params, x=sp.wave)
resid = sp.flux - out

# result.plot_fit()
# plt.show()

# #################################################################################

voigt4 = VoigtModel(prefix='voigt4_')
voigt4_pars = voigt4.guess(resid, x=sp.wave)


model = model * voigt4
pars = result.params + voigt4_pars

result = model.fit(data=sp.flux, params=pars, x=sp.wave, method=method)
out = model.eval(data=sp.flux, params=result.params, x=sp.wave)
resid = sp.flux - out


# #################################################################################

print(result.fit_report())

plt.plot(sp.wave, sp.flux)
plt.plot(sp.wave, out)
plt.plot(sp.wave, resid)
plt.show()
