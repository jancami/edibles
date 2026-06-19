import matplotlib.pyplot as plt
from edibles.projects.DIBCode_Empirical.importdata import importdata
from edibles.projects.DIBCode_Empirical.continuumwindow import continuum_window
from edibles.projects.DIBCode_Empirical.run_model_fitting import run_model_fitting
from edibles.projects.DIBCode_Empirical.plot_continuum_removed import plot_continuum_removed
from edibles.projects.DIBCode_Empirical.plotall import plotall
from edibles.projects.DIBCode_Empirical.try_orders import try_orders

minrange=6268
maxrange=6272
c0=1.0226
c1=0.0037
data_piece="564nm_redu_O10"
ContinuumMin=6253
ContinuumMax=6256
wavelength_target=6270
csv_name=f"OliverHD170740components6270.csv"
#try_orders(target="HD 170740",minrange=minrange, maxrange=maxrange, ContinuumMin = ContinuumMin, ContinuumMax = ContinuumMax)

recalibrated_wavelength, coadded_flux, coadd_SNR, target, common_wave_full, coadded_flux_full = importdata(target="HD 170740",minrange=minrange, maxrange=maxrange, data_piece = data_piece, ContinuumMin = ContinuumMin, ContinuumMax = ContinuumMax)
plt.plot(recalibrated_wavelength, coadded_flux)
plt.xlabel("Wavelength [Å]")
plt.ylabel("Normalized Flux")
plt.title("DIB 6270 Å")
plt.show()

continuum_window(common_wave_full,coadded_flux_full,ContinuumMin = ContinuumMin, ContinuumMax = ContinuumMax, wavelength_target=wavelength_target)

run_model_fitting(minrange=minrange, maxrange=maxrange, folder = "/Users/oliverridge/PycharmProjects/EDIBLES/component_data", csv_name = csv_name)

listt = ['HD 170740', 'HD 23180', 'HD 24398', 'HD 144470', 'HD 147165', 'HD 147683', 'HD 149757', 'HD 166937', 'HD 184915', 'HD 185418', 'HD 185859', 'HD 203532']
data_results = plotall(targets = listt, minrange=minrange, maxrange=maxrange, csv_name = f'/Users/oliverridge/PycharmProjects/EDIBLES/component_data/{csv_name}', c0=c0, c1=c1, model_v_shift=10.67, data_piece = data_piece, ContinuumMin = ContinuumMin, ContinuumMax = ContinuumMax)
plot_continuum_removed(data_results)
