from edibles.utils.edibles_oracle import EdiblesOracle
from edibles.utils.edibles_spectrum import EdiblesSpectrum
from astrovoigtfit import *

# calling standard python libraries
import numpy as np
import matplotlib.pyplot as plt
import re


def observations(star,molecule,file_no, wave_range,species_file='species.txt'):
    
    with open(species_file) as f:
        lines = f.readlines()

    headers = lines[0].split()
    line_index = headers.index('line')

    for line in lines[1:]:
        parts = line.split()
        if parts[0] == molecule:
            line_val = parts[line_index].strip('[]')
            line_val = float(line_val)
            # print(line_val)
            # print(wave_range)
            
    pythia = EdiblesOracle()
    List = pythia.getFilteredObsList(object=[star], MergedOnly=False, Wave=line_val)

    # print(List)
    test = List.values.tolist()
    # print(test)
    
    filename = test[file_no]
    print(filename)
    wrange = wave_range
    sp = EdiblesSpectrum(filename)
    sp.getSpectrum(wrange[0],wrange[1])
    wave = sp.bary_wave
    flux = sp.bary_flux
    idx = np.where((wave > wrange[0]) & (wave < wrange[1]))
    wave = wave[idx]
    flux = flux[idx]
    flux = flux / np.median(flux)

    # Filter the spectrum within the specified wavelength range
    idx = np.where((wave > wrange[0]) & (wave < wrange[1]))
    wave = wave[idx]
    flux = flux[idx]

    # Normalize the flux
    flux = flux / np.median(flux)
    
    return wave, flux

def get_available_files(star, molecule, species_file='species.txt'):
    """Retrieves the list of available observation files for a given star and molecule."""
    with open(species_file) as f:
        lines = f.readlines()

    headers = lines[0].split()
    line_index = headers.index('line')
    
    line_val = None
    for line in lines[1:]:
        parts = line.split()
        if parts[0] == molecule:
            line_val = parts[line_index].strip('[]')
            line_val = float(line_val)
            break
            
    if line_val is None:
        return []

    pythia = EdiblesOracle()
    List = pythia.getFilteredObsList(object=[star], MergedOnly=False, Wave=line_val)
    return List.values.tolist()

def get_species_data(species_file='species.txt'):
    """Parses species.txt to return available species and their wavelength ranges."""
    species_data = {}
    with open(species_file) as f:
        lines = f.readlines()
        
    # Assuming header is line 0
    # Format: species line lambda f_value gamma wrange
    
    for line in lines[1:]:
        parts = line.split()
        if len(parts) < 6: continue
        
        name = parts[0]
        # Parse wrange: use last element as it's the last column
        wrange_str = parts[-1].strip('[]')
        try:
            wrange = [float(x) for x in wrange_str.split(',')]
            line_val = parts[1]
            species_data[name] = {'wrange': wrange, 'line': line_val}
        except:
            pass
            
    return species_data
    
    
def continuum_fit(star, molecule,file_no, wave_range, absorption_range,degree=3):
    
    # Get the observed spectrum
    wave, flux = observations(star, molecule[0], file_no, wave_range)

    absorption_range = absorption_range  # Wavelength range of the absorption feature
    

    #  Fit continuum and normalize
    continuum_normalized_flux, continuum, poly, std_dev = fit_continuum(
        wave, flux, absorption_range, degree, return_std=True
    )  
    
    return wave,continuum_normalized_flux ,flux, continuum  

    # # plotting the normalized spectrum and the continuum curve
    # plt.figure(figsize=(8, 8))
    # plt.subplot(2, 1, 1)
    # plt.plot(wave, flux, 'b-', label='Original Flux')
    # plt.plot(wave, continuum, 'r-', label='Continuum Fit')
    # plt.axvspan(absorption_range[0], absorption_range[1], color='gray', alpha=0.3, label='Absorption Region')
    # plt.xlabel('Wave')
    # plt.ylabel('Flux')
    # plt.legend()
    # plt.title('Continuum Fitting')

    # # fitting continuum Normalized flux
    # plt.subplot(2, 1, 2)
    # plt.plot(wave, continuum_normalized_flux, 'g-', label='Normalized Flux')
    # plt.axhline(1.0, color='k', linestyle='--', label='Continuum = 1.0')
    # plt.axvspan(absorption_range[0], absorption_range[1], color='gray', alpha=0.3)
    # plt.xlabel('Wavelength')
    # plt.ylabel('Normalized Flux')
    # plt.legend()
    # plt.title('Normalized Spectrum')


def _parse_bracketed_list(s):
    inner = s.strip().strip('[]')
    if not inner:
        return []
    return [float(x.strip()) for x in inner.split(',') if x.strip()]

def get_species_params(species_file, species_params, molecule):
    with open(species_file) as f:
        lines = f.readlines()

    for idx, species in enumerate(molecule):
        for line in lines[1:]:
            if not line.strip():
                continue
            # species name is first whitespace-separated token
            parts = line.split()
            if parts[0] != species:
                continue

            # extract the three bracketed columns (lambda, f_value, gamma)
            bracketed = re.findall(r'\[.*?\]', line)
            if len(bracketed) < 3:
                raise ValueError(f"Line for species {species} does not have expected bracketed fields: {line!r}")

            lambda_vals = _parse_bracketed_list(bracketed[0])
            f_vals = _parse_bracketed_list(bracketed[1])
            gamma_vals = _parse_bracketed_list(bracketed[2])

            # Handle case where lambda is [0] (e.g. K_4044)
            if len(lambda_vals) == 1 and lambda_vals[0] == 0:
                try:
                    line_val = float(parts[1])
                    lambda_vals = [line_val]
                    print(f"Warning: Lambda is 0 for {species}, using line value {line_val}")
                except:
                    pass

            print(f"Species: {species}, Lambda: {lambda_vals}, f: {f_vals}, Gamma: {gamma_vals}")

            reordered = {
                'lambda': lambda_vals,
                'f': f_vals,
                'gamma': gamma_vals,
            }

            for key, value in species_params[idx].items():
                reordered[key] = value

            species_params[idx] = reordered
            break  # found species line, move to next

    return species_params




def astrovoigtfit_run(star,molecule, wave_range,species_params,absorption_range,file_no,species_file='species.txt', degree=3):
    
    # Get the observed spectrum and fit the continuum
    
    wave, flux = observations(star,molecule[0],file_no, wave_range,species_file)
    species_param_updated = get_species_params(species_file, species_params, molecule)
    # print("Species parameters:", species_param_updated)
    
    absorption_range = absorption_range  # Wavelength range of the absorption feature
    # degree = 3  # Degree of Chebyshev polynomial - Removed hardcoding

    #  Fit continuum and normalize
    continuum_normalized_flux, continuum, poly, std_dev = fit_continuum(
        wave, flux, absorption_range, degree, return_std=True
    )     
    

    # fitting the data using astro_voigt_fit function
    fitresult= astro_voigt_fit(
        wavegrid=wave, 
        ydata=continuum_normalized_flux,
        species_params=species_param_updated,
        v_resolution=3, 
        n_step=25, 
        std_dev=0.002
    )
    fitresult.params.pretty_print() #printing the fitting parameters



    print("chi-square value ",fitresult.chisqr)
    print("reduced chi-square value ",fitresult.redchi)
    print("FITTING RESULT :", fitresult.success)


    plt.plot(wave,fitresult.best_fit,color ='purple',label ="fit")
    plt.plot(wave, continuum_normalized_flux,color ='gray',label ='data',alpha = 0.7)
    plt.xlabel("Wavelength ($\AA$)")
    plt.ylabel("Normalised flux")
    plt.title("Multi cloud single line model fit for CH+",color = 'darkgreen')
    plt.grid()
    plt.legend()
    plt.show()

def astrovoigtfit_general_run(wave, flux, molecules, species_params, absorption_range, species_file='species.txt', degree=3):
    """
    General run function for fitting user-supplied wave and flux data.
    """
    
    # Populate species parameters with atomic data (lambda, f, gamma)
    species_param_updated = get_species_params(species_file, species_params, molecules)
    
    # degree = 3  # Degree of Chebyshev polynomial - Removed hardcoding

    # Fit continuum and normalize
    continuum_normalized_flux, continuum, poly, std_dev = fit_continuum(
        wave, flux, absorption_range, degree, return_std=True
    )     
    
    # Fitting the data using astro_voigt_fit function
    fitresult = astro_voigt_fit(
        wavegrid=wave, 
        ydata=continuum_normalized_flux,
        species_params=species_param_updated,
        v_resolution=3, 
        n_step=25, 
        std_dev=0.002 # You might want to allow this to be passed in or estimated
    )
    
    fitresult.params.pretty_print() # Printing the fitting parameters

    print("chi-square value ", fitresult.chisqr)
    print("reduced chi-square value ", fitresult.redchi)
    print("FITTING RESULT :", fitresult.success)

    # Plotting
    plt.figure(figsize=(10, 6))
    plt.plot(wave, fitresult.best_fit, color='purple', label="Fit")
    plt.plot(wave, continuum_normalized_flux, color='gray', label='Data', alpha=0.7)
    plt.xlabel("Wavelength ($\AA$)")
    plt.ylabel("Normalized Flux")
    plt.title(f"General Model Fit for {', '.join(molecules)}", color='darkgreen')
    plt.grid(True)
    plt.legend()
    # plt.show()  # Commented out - GUI displays plots, no need for popup
    plt.close()  # Close the figure to free memory
    
    return fitresult, continuum_normalized_flux, continuum, std_dev

    import os
import numpy as np

def save_data(star, molecule, file_no, wave_range,
              species_file='species.txt',
              filetype='txt',
              out_dir=None,
              out_name=None):
    """
    Get wave & flux from observations() and save as 2-column text file.

    Parameters
    ----------
    star : str
    molecule : str
    file_no : int
    wave_range : list or tuple, [min_wave, max_wave]
    species_file : str, optional
        Path to species file (passed to observations()).
    filetype : {'txt', 'dat'}, optional
        Output file type/extension. Default is 'txt'.
    out_dir : str or None, optional
        Directory to save the file in.
        - If None (default), saves in the current working directory.
        - If given, saves in that directory (creates it if needed).
    out_name : str or None, optional
        Output filename (without directory). If None, a name like
        '{star}_{molecule}_{file_no}.{filetype}' is used.
    """

    # Get wave and flux from your existing function
    wave, flux = observations(star, molecule, file_no, wave_range, species_file=species_file)

    # Make sure filetype is valid and clean
    filetype = filetype.lower().lstrip('.')
    if filetype not in ('txt', 'dat'):
        raise ValueError("filetype must be 'txt' or 'dat'")

    # Default filename if user doesn't pass one
    if out_name is None:
        out_name = f"{star}_{molecule}_{file_no}.{filetype}"
    else:
        # If user gave a name without extension, add it
        root, ext = os.path.splitext(out_name)
        if ext == "":
            out_name = f"{out_name}.{filetype}"

    # Decide where to save
    if out_dir is None:
        # Same location as where you run the code (current working directory)
        out_path = out_name
    else:
        os.makedirs(out_dir, exist_ok=True)
        out_path = os.path.join(out_dir, out_name)

    # Stack into 2 columns: wave, flux
    data = np.column_stack((wave, flux))

    # Save as text file with a header line
    np.savetxt(out_path, data, fmt="%.8f", header="wave flux", comments='')

    print(f"Saved data to: {out_path}")
    return out_path

def astrovoigtfit_general_run_GM2(wave, flux, molecules, species_params, absorption_range, n_knots=10, knots_x_array=None, species_file='species.txt', spline_order=3):
    """
    General Mode 2 run function for simultaneous continuum and Voigt fitting.
    """
    
    # Populate species parameters with atomic data (lambda, f, gamma)
    species_param_updated = get_species_params(species_file, species_params, molecules)
    
    # Perform Simultaneous Fit
    fitresult = astro_simultaneous_fit(
        wavegrid=wave, 
        ydata=flux,
        species_params=species_param_updated,
        n_knots=n_knots,
        knots_x_array=knots_x_array,
        v_resolution=3, 
        n_step=25, 
        std_dev=0.002,
        absorption_ranges=[absorption_range] if absorption_range else None,
        spline_order=spline_order
    )
    
    fitresult.params.pretty_print() # Printing the fitting parameters

    print("chi-square value ", fitresult.chisqr)
    print("reduced chi-square value ", fitresult.redchi)
    print("FITTING RESULT :", fitresult.success)

    # Extract Continuum from the model
    # The model is Continuum * Absorption. 
    # We can reconstruct the continuum using the spline parameters.
    
    used_knots_x = fitresult.userkws.get('knot_x_array')
    
    # 2. Get knot_y values from params
    knot_y_values = []
    i = 0
    while True:
        name = f'knot_y_{i}'
        if name in fitresult.params:
            knot_y_values.append(fitresult.params[name].value)
            i += 1
        else:
            break
            
    # 3. Evaluate Spline
    # 3. Evaluate Spline
    from scipy.interpolate import CubicSpline, make_interp_spline
    if used_knots_x is not None and len(knot_y_values) == len(used_knots_x):
        if spline_order == 3:
            spline = CubicSpline(used_knots_x, knot_y_values)
        else:
            spline = make_interp_spline(used_knots_x, knot_y_values, k=spline_order)
            
        continuum = spline(wave)
    else:
        # Fallback if something is weird
        print("Warning: Could not reconstruct continuum perfectly.")
        continuum = np.ones_like(wave)

    normalized_flux = flux / continuum
    
    # Calculate std dev of residuals (approx)
    std_dev = np.std(flux - fitresult.best_fit)

    # Plotting
    plt.figure(figsize=(10, 6))
    plt.plot(wave, flux, color='gray', label='Data', alpha=0.7)
    plt.plot(wave, fitresult.best_fit, color='purple', label="Total Fit")
    plt.plot(wave, continuum, color='orange', linestyle='--', label="Continuum")
    plt.xlabel("Wavelength ($\AA$)")
    plt.ylabel("Flux")
    plt.title(f"GM2 Simultaneous Fit for {', '.join(molecules)}", color='darkgreen')
    plt.grid(True)
    plt.legend()
    # plt.show()
    plt.close()
    
    return fitresult, normalized_flux, continuum, std_dev
