import matplotlib.pyplot as plt
from edibles.projects.DIBCode_Empirical.importdata import importdata
from edibles.projects.DIBCode_Empirical.continuum_window import continuum_window
from edibles.projects.DIBCode_Empirical.perform_fit import perform_fit
from edibles.projects.DIBCode_Empirical.run_model_fitting import run_model_fitting
from edibles.projects.DIBCode_Empirical.plot_continuum_removed import plot_continuum_removed
from edibles.projects.DIBCode_Empirical.plotall import plotall
from edibles.projects.DIBCode_Empirical.try_orders import try_orders
from edibles import DATADIR

from edibles.projects.DIBCode_Empirical.format_params_grouped import format_params_grouped

minrange=6200
maxrange=6208.5
c0=1.016831287649225
c1=0
data_piece="564nm_redu_O9"
ContinuumMin=6190
ContinuumMax=6194
wavelength_target=6203
csv_name=f"OliverHD170740components{wavelength_target}.csv"
removed_components=["l0","g4"]
#try_orders(target="HD 170740",minrange=minrange, maxrange=maxrange, ContinuumMin = ContinuumMin, ContinuumMax = ContinuumMax)
#try_orders is used to check each order
figpath = f"/Users/oliverridge/PycharmProjects/EDIBLES/Final_Figures/DIB {wavelength_target}"
recalibrated_wavelength, coadded_flux, coadd_SNR, target, common_wave_full, coadded_flux_full = importdata(target="HD 170740",minrange=minrange, maxrange=maxrange, data_piece = data_piece, ContinuumMin = ContinuumMin, ContinuumMax = ContinuumMax)
plt.plot(recalibrated_wavelength, coadded_flux)
plt.xlabel("Wavelength [Å]")
plt.ylabel("Normalized Flux")
plt.title(f"DIB {wavelength_target} Å")
plt.savefig(f"{figpath}/Normalized_Flux.png")
plt.show()

continuum_window(common_wave_full,coadded_flux_full,ContinuumMin = ContinuumMin, ContinuumMax = ContinuumMax, wavelength_target=wavelength_target)

run_model_fitting(minrange=minrange, maxrange=maxrange, folder = "/Users/oliverridge/PycharmProjects/EDIBLES/component_data", csv_name = csv_name, removed_components=removed_components, recalibrated_wavelength=recalibrated_wavelength,coadded_flux=coadded_flux, coadd_SNR=coadd_SNR, figpath=figpath) # Change the folder component to your path that you want component data to go into

listt = ['HD 170740', 'HD 23180', 'HD 24398', 'HD 144470', 'HD 147165', 'HD 147683', 'HD 149757', 'HD 166937', 'HD 184915', 'HD 185418', 'HD 185859', 'HD 203532']
data_results = plotall(targets = listt, minrange=minrange, maxrange=maxrange, csv_name = f'/Users/oliverridge/PycharmProjects/EDIBLES/component_data/{csv_name}', c0=c0, c1=c1, model_v_shift=10.67, data_piece = data_piece, ContinuumMin = ContinuumMin, ContinuumMax = ContinuumMax, figpath=figpath)
plot_continuum_removed(data_results,figpath)
# Make output for Git repository
#DIB 6270 - minrange=6268, maxrange=6272, c0=1.0226, c1=0.0037, data_piece = "564nm_redu_O10", ContinuumMin = 6253, ContinuumMax = 6256
#DIB 6203 - minrange=6200, maxrange=6208.25, c0=1.016831287649225, c1=0, data_piece= "564nm_redu_O9", ContinuumMin = 6190, ContinuumMax = 6194
#DIB 6613 - minrange=6611, maxrange=6616, c0=1.0315, c1=0.01, data_piece = "564nm_redu_O16", ContinuumMin = 6600, ContinuumMax = 6604