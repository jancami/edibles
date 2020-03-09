import edibles.src.edibles_spectrum as eS
import edibles.src.model as eM
import numpy as np
import matplotlib.pyplot as plt
from sherpa import models
import edibles.edibles_settings as Setting

## Step 1: load spectrum, cut to desire wavelength, and add telluric reference
sightline = "HD36483"
date="20140924"
order="O12"
filenames = "/HD36486/RED_860/HD36486_w860_redl_20140924_O12.fits"
sp = eS.EdiblesSpectrum(filenames, panel_name=order)
sp.cutSpectrum(xmin=7667, xmax=7687)
sp.addLinelist(species="Telluric", wavelengths=[7670.602, 7671.664, 7676.563, 7677.617, 7682.755, 7683.800])

## Step 2: add masks, separate into 3 panels (for 3 pairs of telluric), and continuum fitting
sp.addMask(n=3)
sp.duplicatePanels(ntimes=2)
sp.cutSpectrum_visual()
sp.renamePanels(["Telluric_1", "Telluric_2", "Telluric_3"])
sp.fitContinuum(mode="p", n=3, apply_mask=True)

# Step 3: check spectrum before fitting
sp.resetMask()
sp.converToVelocity(center=sp.linelist["Telluric"][np.array([0,2,4])])
sp.cutSpectrum(xmin=-50, xmax=100)
sp.showSpectrum()

# Step 4: make the model and fit
# the model contains 6 telluric lines with known wav and a varying v_offset
CstCont = models.Const1D()
model_working = eM.Sightline(sightline, CstCont)
model_working.addLineSeries(sp.linelist["Telluric"], "known_wave", line_name="Telluric")
model2fit = model_working.model
sp.importModel(model2fit)
sp.fitModel()

# Step 5: output results
# build-in outputs
sp.outputResult(filename="HD170740_"+date+"_"+order+".txt")
sp.showSpectrum(x_label=sightline,save=True,filename=sightline+"_"+order)
# output the v_offset to a separate file
(name, v_offset) = sp.raiseParameter(par="v_offset")
result_data = Setting.resultdir + "telluric.txt"
f = open(result_data,"a")
#f.write("Date, Order, 7670.602, 7671.664, 7676.563, 7677.617, 7682.755, 7683.800\n")
v_str = np.array2string(v_offset, separator=", ", precision=3)[1:-2]
f.write(", ".join([date,order,v_str])+"\n")
f.close()
