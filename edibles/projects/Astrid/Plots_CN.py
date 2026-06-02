from edibles.utils.edibles_oracle import EdiblesOracle
from edibles.utils.edibles_spectrum import EdiblesSpectrum
import matplotlib.pyplot as plt
import numpy as np
import os


print("Inicializando Oracle...")
pythia = EdiblesOracle()


file_list = pythia.getFilteredObsList(
   object=None,
   Wave=3874,
   closest_order=True
)


print("Total espectros:", len(file_list))


CN_lines = [3873.99, 3874.61, 3875.76]
useful_files = []


# -------- ANALISIS --------
for file in file_list:


   try:
       spec = EdiblesSpectrum(file)
   except:
       continue


   mask = (spec.wave > 3870) & (spec.wave < 3878)
   wave = spec.wave[mask]
   flux = spec.flux[mask]


   if len(wave) < 10:
       continue


   dips_detected = 0


   for line in CN_lines:
       region = (wave > line - 0.2) & (wave < line + 0.2)


       if np.sum(region) > 5:
           local_flux = flux[region]
           mean_flux = np.mean(flux)


           if np.min(local_flux) < 0.98 * mean_flux:
               dips_detected += 1


   if dips_detected >= 2:
       useful_files.append(file)


print("Total útiles:", len(useful_files))




# -------- GUARDAR GRAFICAS --------
output_folder = "CN_plots"
os.makedirs(output_folder, exist_ok=True)


print("Guardando espectros...")


for i, file in enumerate(useful_files):


   try:
       spec = EdiblesSpectrum(file)
   except:
       continue


   mask = (spec.wave > 3870) & (spec.wave < 3878)
   wave = spec.wave[mask]
   flux = spec.flux[mask]


   plt.figure(figsize=(8,5))
   plt.plot(wave, flux, color='black')


   # marcar CN
   for line in CN_lines:
       plt.axvline(x=line, color='red', linestyle='--')


   plt.xlabel("Wavelength (Å)")
   plt.ylabel("Flux")
   plt.title(file)


   filename = file.replace("/", "_") + ".png"
   filepath = os.path.join(output_folder, filename)


   plt.savefig(filepath)
   plt.close()


   print(f"{i+1}/{len(useful_files)} guardado")


print("✅ TODOS LOS ESPECTROS GUARDADOS")

