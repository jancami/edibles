import matplotlib.pyplot as plt
import os
from edibles.projects.DIBCode_Empirical.importdata import importdata
from edibles.projects.DIBCode_Empirical.continuum_window import continuum_window
from edibles.projects.DIBCode_Empirical.perform_fit import perform_fit
from edibles.projects.DIBCode_Empirical.run_model_fitting import run_model_fitting
from edibles.projects.DIBCode_Empirical.plot_continuum_removed import plot_continuum_removed
from edibles.projects.DIBCode_Empirical.plotall import plotall
from edibles.projects.DIBCode_Empirical.try_orders import try_orders
from edibles import EDIBLES_OUTPUTDIR

from edibles.projects.DIBCode_Empirical.format_params_grouped import format_params_grouped

minrange=4758
maxrange=4764
c0=1.0272
c1=0.0004
data_piece="564nm_redl_O5"
ContinuumMin=4768
ContinuumMax=4772
wavelength_target=4762
csv_name=f"HD170740components{wavelength_target}.csv"
removed_components=[]
#try_orders(target="HD 170740",minrange=minrange, maxrange=maxrange, ContinuumMin = ContinuumMin, ContinuumMax = ContinuumMax, filestring="564nm_redl_O")
#try_orders is used to check each order, useful for finding the correct order for new DIBs
figpath = f"{EDIBLES_OUTPUTDIR}/DIB {wavelength_target}"
if not os.path.exists(figpath):
    os.makedirs(figpath)
# If the figure path and directories do not exist they will have to be made
recalibrated_wavelength, coadded_flux, coadd_SNR, target, common_wave_full, coadded_flux_full = importdata(target="HD 170740",minrange=minrange, maxrange=maxrange, data_piece = data_piece, ContinuumMin = ContinuumMin, ContinuumMax = ContinuumMax)
plt.plot(recalibrated_wavelength, coadded_flux)
plt.xlabel("Wavelength [Å]")
plt.ylabel("Normalized Flux")
plt.title(f"DIB {wavelength_target} Å")
plt.savefig(f"{figpath}/Initial_DIB_Data.png")
plt.show()

continuum_window(common_wave_full,coadded_flux_full,ContinuumMin = ContinuumMin, ContinuumMax = ContinuumMax, wavelength_target=wavelength_target, figpath=figpath)

run_model_fitting(minrange=minrange, maxrange=maxrange, folder = f"{EDIBLES_OUTPUTDIR}/component_data", csv_name = csv_name, removed_components=removed_components, recalibrated_wavelength=recalibrated_wavelength,coadded_flux=coadded_flux, coadd_SNR=coadd_SNR, figpath=figpath, wavelength_target=wavelength_target)
# Change the folder component to your path that you want component data to go into

listt = ['HD 170740', 'HD 23180', 'HD 24398', 'HD 147683', 'HD 149757', 'HD 166937', 'HD 185418', 'HD 185859', 'HD 203532']
#listt = ['HD 170740', 'HD 23180', 'HD 24398', 'HD 144470', 'HD 147165', 'HD 147683', 'HD 149757', 'HD 166937', 'HD 184915', 'HD 185418', 'HD 185859', 'HD 203532']
data_results = plotall(targets = listt, minrange=minrange, maxrange=maxrange, csv_name = f'{EDIBLES_OUTPUTDIR}/component_data/{csv_name}', c0=c0, c1=c1, model_v_shift=10.67, data_piece = data_piece, ContinuumMin = ContinuumMin, ContinuumMax = ContinuumMax, figpath=figpath)
plot_continuum_removed(data_results,figpath)
# Add output for Git repository
# Add documentation for newly added files


#DIB 6270 - minrange=6267.5, maxrange=6272, c0=1.0272, c1=0.0004, data_piece = "564nm_redu_O10", ContinuumMin = 6253, ContinuumMax = 6256, removed_components=["g3","g4","g5"]
#DIB 6203 - minrange=6199, maxrange=6208.25, c0=1.016831287649225, c1=0, data_piece= "564nm_redu_O9", ContinuumMin = 6190, ContinuumMax = 6194, removed_components=["l1","l2","g3"]
#DIB 6613 - minrange=6611, maxrange=6616, c0=1.0315, c1=0.01, data_piece = "564nm_redu_O16", ContinuumMin = 6600, ContinuumMax = 6604, removed_components=["g4","g5"]
#DIB 4762 - minrange=4758, maxrange=4766, c0=1.0035, c1=-0.0045, data_piece = "564nm_redl_O5", ContinuumMin = 4768, ContinuumMax = 4772, removed_components=["g1"]