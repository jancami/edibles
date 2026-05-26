"""Loads one spectrum from the Edibles archive.

wave_flux_data is what the GUI actually calls (once per selected FITS file in
the co-add step). It returns a wavelength-cropped, median-normalized chunk of
spectrum, with the file looked up either by name (the usual case from the GUI)
or by index into the EdiblesOracle results for a star + transition wavelength.

wave_flux_info is just the index-lookup helper — wave_flux_data uses it
internally when you pass a file_number instead of a file_name.

The output (wave, flux) pairs get fed into co_adding_flux.coadd_spectra next.
"""

import numpy as np
from edibles.utils.edibles_spectrum import EdiblesSpectrum
from edibles.utils.edibles_oracle import EdiblesOracle


def wave_flux_info(star_name=None, wave=None):
    pythia = EdiblesOracle()
    List = pythia.getFilteredObsList(object=[star_name], MergedOnly=False, Wave=wave)
    test = List.values.tolist()
    
    return test

def wave_flux_data(star_name=None, wave=None, wrange=None, file_number=None, file_name=None):
    if file_name is None:
        if file_number is None:
            print("Error: Must provide either file_name or file_number")
            return None, None
        test = wave_flux_info(star_name=star_name, wave=wave)
        if file_number >= len(test):
            print(f"Error: file_number {file_number} out of range for {len(test)} files found.")
            return None, None
        file_name = test[file_number]
    
    # SANITIZE FILENAME FOR WINDOWS - replace colons with hyphens
    if isinstance(file_name, str):
        file_name = file_name.replace(':', '_')
    
    sp = EdiblesSpectrum(file_name)
    try:
        sp.getSpectrum(xmin=wrange[0], xmax=wrange[1])
    except Exception as e:
        print(f"Warning: getSpectrum failed with range {wrange}: {e}")
        # Fallback: Check bounds and clamp
        if hasattr(sp, 'wave') and sp.wave is not None:
             file_min = np.min(sp.wave)
             file_max = np.max(sp.wave)
        elif hasattr(sp, 'raw_wave') and sp.raw_wave is not None:
             file_min = np.min(sp.raw_wave)
             file_max = np.max(sp.raw_wave)
        else:
             print("Error: Could not determine file bounds.")
             return None, None
             
        # Epsilon buffer to ensure we are strictly inside bounds
        epsilon = 0.05 
        req_min = max(wrange[0], file_min + epsilon)
        req_max = min(wrange[1], file_max - epsilon)
        
        if req_min < req_max:
             print(f"Retrying with clamped range: {req_min:.3f} - {req_max:.3f}")
             try:
                 sp.getSpectrum(xmin=req_min, xmax=req_max)
             except Exception as e2:
                 print(f"Error in retry: {e2}")
                 return None, None
        else:
             print(f"Error: No overlap between requested {wrange} and file {file_min}-{file_max}")
             return None, None
    wave_new = sp.bary_wave
    flux_new = sp.bary_flux
    idx = np.where((wave_new > wrange[0]) & (wave_new < wrange[1]))
    wave = wave_new[idx]
    flux = flux_new[idx]

    # Normalize the flux (use the subset 'flux', not 'flux_new')
    if len(flux) > 0:
        flux = flux / np.median(flux)
    else:
        return None, None
    
    return wave, flux
    
    
    
    
    



