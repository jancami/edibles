import numpy as np
import matplotlib.pyplot as plt
import os
from edibles.utils.edibles_spectrum import EdiblesSpectrum
from edibles.utils.edibles_oracle import EdiblesOracle
from edibles.continuum import Continuum
from edibles import DATADIR

pythia = EdiblesOracle()
dirs = [f for f in os.listdir(DATADIR) if not f.startswith(".")]
figure_output = DATADIR.replace("/data", "/Na3302Plots")
if not os.path.isdir(figure_output):
    os.mkdir(figure_output)

for i, sightline in enumerate(dirs):
    print("="*20)
    print("Now working on Sightline " + sightline)
    print("Sightline {i}/{j}".format(i=i+1, j=len(dirs)))
    print("="*20)

    matched_files = pythia.getFilteredObsList(object=[sightline.replace("HD", "HD ")], OrdersOnly=True, Wave=3302)
    figure_output_sightline = os.path.join(figure_output, sightline)
    if not os.path.isdir(figure_output_sightline):
        os.mkdir(figure_output_sightline)

    for file in matched_files:
        sp = EdiblesSpectrum(file)
        sp.getSpectrum(xmin=3301, xmax=3305)
        date = sp.datetime.date()
        order = file.split("_")[-1].replace(".fits", "")

        # build a 4 anchor points spline
        cont = Continuum(sp, method="spline", n_anchors=4, plot=False, verbose=0)
        params = cont.model.guess(sp.flux, x=sp.wave)
        out = cont.model.eval(params=params, x=sp.wave)
        cont.add_to_csv(
            user="Na Worker", comments="1st order continuum for Na 3302"
        )

        # save plots
        filename = "_".join([sightline, str(date), order]) + ".png"
        filename = "/".join([figure_output_sightline, filename])
        plt.plot(sp.wave, sp.flux)
        plt.plot(sp.wave, out)
        plt.ylabel(sightline)
        plt.xlabel("  ".join([str(date), order]))
        plt.savefig(filename)
        plt.close()



