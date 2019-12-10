from __future__ import print_function

from edibles.functions.edibles_spectrum import EdiblesSpectrum
from edibles.fit.models.create_model import createCont
from edibles.fit.models.models import Sightline
from edibles.fit.fit import fit


file1 = '/HD170740/RED_860/HD170740_w860_redl_20140915_O12.fits'
sp1 = EdiblesSpectrum(file1)
xmin = 7661.5
xmax = 7669.0


data = sp1.getSpectrum(xmin=xmin,xmax=xmax, bary=False)

cont = createCont(data, n_points=4)
sightline = Sightline(star_name=sp1.target, cont=cont)
sightline.addSource(source_name='Telluric', b=0.001, d=0.05)
sightline.addLine(name='tell1', lam_0=7664.5, tau_0=0.6)
sightline.addLine(name='tell2', lam_0=7666, tau_0=0.6)


sightline.addSource(source_name='Interstellar', b=0.01, d=0.02)
sightline.addLine(name='int1', lam_0=7665.2, tau_0=0.2)
sightline.addLine(name='int2', lam_0=7665.3, tau_0=0.01)

fit_m=fit(sp1.target,data,sightline.model,breakdown=False,silent=False)

print([line for line in sightline.lines['all']])

