import numpy as np
from numpy.polynomial.chebyshev import Chebyshev
from scipy.interpolate import CubicSpline
from lmfit import Parameters, Model
from main import master_function




# first we will fit the continuum
def fit_continuum(wavelength, flux, absorption_range, degree, return_std=False):
    """
    Fit a Chebyshev polynomial to the continuum, excluding an absorption region.
    
    Parameters:
    - wavelength: array of wavelengths
    - flux: array of flux values
    - absorption_range: tuple (start, end) defining the absorption region to exclude
    - degree: degree of the Chebyshev polynomial
    - return_std: if True, returns standard deviation of continuum points outside absorption region
    
    Returns:
    - normalized_flux: array of normalized flux values (flux/continuum)
    - continuum: array of fitted continuum values for all wavelengths
    - poly: the fitted Chebyshev polynomial object
    - std (optional): standard deviation of points outside absorption region
    """
    # Mask the absorption region
    mask = (wavelength < absorption_range[0]) | (wavelength > absorption_range[1])
    wl_continuum = wavelength[mask]
    flux_continuum = flux[mask]
    
    # Fit Chebyshev polynomial
    poly = Chebyshev.fit(wl_continuum, flux_continuum, deg=degree)
    
    # Evaluate continuum over the entire wavelength range
    continuum = poly(wavelength)
    normalized_flux = flux / continuum
    
    if return_std:
        # Calculate standard deviation only for points outside absorption region
        std = np.std(normalized_flux[mask] - 1)  # Subtract 1 because continuum is normalized to 1
        return normalized_flux, continuum, poly, std
    else:
        return normalized_flux, continuum, poly



# --- Wrapper for model evaluation ---
def Voigt_fit_wrapper(**params_list):
    """
    Generalized Voigt fitting wrapper that works with any number of species.
    
    Expected parameters:
    - n_species: number of species
    - wavegrid: wavelength grid
    - v_resolution: velocity resolution
    - n_step: number of steps
    - For each species i: n_trans_i, n_component_i
    - For each transition j in species i: lambda_i_j, f_i_j, gamma_i_j
    - For each component k in species i: v_rad_i_k, b_i_k, N_i_k
    """
    n_species = params_list['n_species']
    wavegrid = params_list['wavegrid']
    v_resolution = params_list['v_resolution']
    n_step = params_list['n_step']
    
    # Prepare data structures for all species
    species_data = {}
    
    for species_idx in range(int(n_species)):
        n_trans = params_list[f'n_trans_{species_idx}']
        n_component = params_list[f'n_component_{species_idx}']
        n_component = int(n_component)
        
        # Extract transition parameters for this species
        n_trans = int(float(n_trans)) 
        all_lambda = np.empty(n_trans)
        all_f = np.empty(n_trans)
        all_gamma = np.empty(n_trans)
        
        for i in range(n_trans):
            all_lambda[i] = params_list[f'lambda_{species_idx}_{i}']
            all_f[i] = params_list[f'f_{species_idx}_{i}']
            all_gamma[i] = params_list[f'gamma_{species_idx}_{i}']
        
        # Extract component parameters for this species
        all_v_rad = np.empty(n_component)
        all_b = np.empty(n_component)
        all_N = np.empty(n_component)
        
        for i in range(n_component):
            all_v_rad[i] = params_list[f'v_rad_{species_idx}_{i}']
            all_b[i] = params_list[f'b_{species_idx}_{i}']
            all_N[i] = params_list[f'N_{species_idx}_{i}']
        
        species_data[species_idx] = {
            'lambda': all_lambda,
            'f': all_f,
            'gamma': all_gamma,
            'v_rad': all_v_rad,
            'b': all_b,
            'N': all_N
        }
    
    # Call the existing master_function with the old parameter format
    # Convert species_data back to the old format that master_function expects
    master_kwargs = {}
    
    for species_idx in species_data:
        # Convert species index to ordinal suffix
        if species_idx == 0:
            suffix = '1st'
        elif species_idx == 1:
            suffix = '2nd'
        elif species_idx == 2:
            suffix = '3rd'
        else:
            suffix = f'{species_idx + 1}th'
        
        # Add parameters in the old format
        master_kwargs[f'lambda_{suffix}'] = species_data[species_idx]['lambda']
        master_kwargs[f'f_{suffix}'] = species_data[species_idx]['f']
        master_kwargs[f'gamma_{suffix}'] = species_data[species_idx]['gamma']
        master_kwargs[f'b_{suffix}'] = species_data[species_idx]['b']
        master_kwargs[f'N_{suffix}'] = species_data[species_idx]['N']
        master_kwargs[f'v_rad_{suffix}'] = species_data[species_idx]['v_rad']
    
    model = master_function(
        wavegrid, 
        v_resolution=v_resolution, 
        n_step=n_step,
        **master_kwargs
    )
    
    return model

def astro_voigt_fit(
    wavegrid, 
    ydata, 
    species_params,
    v_resolution=0.0, 
    n_step=25, 
    std_dev=0.002
):
    """
    Generalized fitting function for multiple species with v_rad constraints.
    
    Parameters:
    -----------
    wavegrid : array
        Wavelength grid
    ydata : array
        Observed data to fit
    species_params : dict
        Dictionary containing parameters for each species.
        Structure: {
            0: {  # species index (Reference species)
                'lambda': [list of wavelengths],
                'f': [list of f values],
                'gamma': [list of gamma values],
                'b': [list of b values for each component],
                'N': [list of N values for each component],
                'v_rad': [list of radial velocities for each component]
            },
            1: {  # next species
                ... parameters ...
                'tie_b': True/False,      # Optional: if True, tie b values to species 0
                'tie_v_rad': True/False   # Optional: if True, tie v_rad values to species 0
            },
            ...
        }

    Example of Tying:
    -----------------
    To tie the b and v_rad parameters of species 1 to species 0:
    
    species_params = {
        0: { ... params for species 0 ... },
        1: { 
            ... params for species 1 ... ,
            'tie_b': True,       # b_1_i will be constrained to b_0_i
            'tie_v_rad': True    # v_rad_1_i will be constrained to v_rad_0_i
        }
    }
    
    Note: Tying assumes 1:1 mapping of components between the tied species and species 0.

    v_resolution : float
        Velocity resolution
    n_step : int
        Number of steps
    std_dev : float
        Standard deviation for weighting
        
    Returns:
    --------
    result : lmfit.model.ModelResult
        Fitting result
    """
    
    n_species = len(species_params)
    
    # converingggg to numpyyy array.......
    for species_idx in species_params:
        for key in ['lambda', 'f', 'gamma', 'b', 'N', 'v_rad']:
            species_params[species_idx][key] = np.asarray(species_params[species_idx][key])
    
    
    params = Parameters()
    params.add('n_species', value=n_species, vary=False)
    params.add('v_resolution', value=v_resolution, vary=False)
    params.add('n_step', value=n_step, vary=False)
    
    # Add parameters for each species
    for species_idx in range(n_species):
        species_data = species_params[species_idx]
        
        # Transition parameters
        n_trans = species_data['lambda'].size
        params.add(f'n_trans_{species_idx}', value=n_trans, vary=False)
        
        for i in range(n_trans):
            params.add(f'lambda_{species_idx}_{i}', value=species_data['lambda'][i], vary=False)
            params.add(f'f_{species_idx}_{i}', value=species_data['f'][i], vary=False)
            params.add(f'gamma_{species_idx}_{i}', value=species_data['gamma'][i], vary=False)
        
        # Component parameters
        n_component = species_data['v_rad'].size
        params.add(f'n_component_{species_idx}', value=n_component, vary=False)
        
        # Check tying flags (default to False if not present)
        # Note: Species 0 cannot be tied to itself in this logic, so we force False for idx 0
        tie_b = species_data.get('tie_b', False) if species_idx > 0 else False
        tie_v_rad = species_data.get('tie_v_rad', False) if species_idx > 0 else False
        
        # Free parameters for fitting
        for i in range(n_component):
            # --- b parameter ---
            if tie_b:
                # Tie to species 0, component i
                # We assume species 0 has at least as many components as this species
                params.add(f'b_{species_idx}_{i}', expr=f'b_0_{i}')
            else:
                # Independent parameter
                params.add(f'b_{species_idx}_{i}', value=species_data['b'][i], min=0.2, max=5.5, vary=True)

            # --- N parameter (always independent) ---
            params.add(f'N_{species_idx}_{i}', value=species_data['N'][i], min=0, vary=True)

            # --- v_rad parameter ---
            if tie_v_rad:
                # Tie to species 0, component i
                params.add(f'v_rad_{species_idx}_{i}', expr=f'v_rad_0_{i}')
            else:
                # Independent parameter
                params.add(f'v_rad_{species_idx}_{i}', value=species_data['v_rad'][i], vary=True)
    
    # Creating model and fitting by passing the inputs to the wrapper........
    voigtmod = Model(Voigt_fit_wrapper, independent_vars=['wavegrid'])
    result = voigtmod.fit(ydata, params, wavegrid=wavegrid, weights=1/std_dev)
    
    return result  # great that you are reading this :)



 
# # for global spectra......
# def dual_spectrum_voigt_fit(
#     wavegrid1, ydata1, species_params1,
#     wavegrid2, ydata2, species_params2,
#     v_resolution=3, 
#     n_step=25, 
#     std_dev1=0.002,
#     std_dev2=0.002
# ):
#     """
#     Wrapper function to fit two different Voigt spectra simultaneously 
#     with constraints that N, b, and v_rad parameters vary equally across both spectra.
    
#     Parameters:
#     -----------
#     wavegrid1, wavegrid2 : array
#         Wavelength grids for spectrum 1 and 2
#     ydata1, ydata2 : array
#         Observed data for spectrum 1 and 2
#     species_params1, species_params2 : dict
#         Species parameters for each spectrum (same structure as astro_voigt_fit)
#     v_resolution : float
#         Velocity resolution
#     n_step : int
#         Number of steps
#     std_dev1, std_dev2 : float
#         Standard deviations for weighting each spectrum
        
#     Returns:
#     --------
#     result : lmfit.model.ModelResult
#         Combined fitting result for both spectra
#     """
    
#     from lmfit import Parameters, Model
#     import numpy as np
    
#     # Validate that both spectra have the same structure......
#     if len(species_params1) != len(species_params2):
#         raise ValueError("Both spectra must have the same number of species")
    
#     for species_idx in species_params1:
#         if species_idx not in species_params2:
#             raise ValueError(f"Species {species_idx} not found in second spectrum")
        
        
#         n_comp1 = species_params1[species_idx]['v_rad'].size if hasattr(species_params1[species_idx]['v_rad'], 'size') else len(species_params1[species_idx]['v_rad'])
#         n_comp2 = species_params2[species_idx]['v_rad'].size if hasattr(species_params2[species_idx]['v_rad'], 'size') else len(species_params2[species_idx]['v_rad'])
        
#         if n_comp1 != n_comp2:
#             raise ValueError(f"Species {species_idx} has different number of components in the two spectra")
    
#     # let's  Convert lists to numpy arrays.......
#     for species_idx in species_params1:
#         for key in ['lambda', 'f', 'gamma', 'b', 'N', 'v_rad']:
#             species_params1[species_idx][key] = np.asarray(species_params1[species_idx][key])
#             species_params2[species_idx][key] = np.asarray(species_params2[species_idx][key])
    
#     n_species = len(species_params1)
    
    
#     params = Parameters()
    
#     # Add global parameters for both spectra
#     params.add('n_species', value=n_species, vary=False)
#     params.add('v_resolution', value=v_resolution, vary=False)
#     params.add('n_step', value=n_step, vary=False)
    
#     # Create master parameters for N, b, and v_rad that will be shared
#     # We'll use the first spectrum's values as initial values
#     master_params = {}
    
#     for species_idx in species_params1:
#         species_data1 = species_params1[species_idx]
#         n_component = species_data1['v_rad'].size
        
#         # Create master parameters for this species
#         for i in range(n_component):
#             # Master b parameter
#             b_param_name = f'master_b_{species_idx}_{i}'
#             params.add(b_param_name, value=species_data1['b'][i], min=0.5, max=5.5)
            
#             # Master N parameter  
#             N_param_name = f'master_N_{species_idx}_{i}'
#             params.add(N_param_name, value=species_data1['N'][i], min=0)
            
#             # Master v_rad parameter
#             v_rad_param_name = f'master_v_rad_{species_idx}_{i}'
#             params.add(v_rad_param_name, value=species_data1['v_rad'][i], vary=True)
            
#             # Store parameter names for later reference
#             master_params[(species_idx, i)] = {
#                 'b': b_param_name,
#                 'N': N_param_name,
#                 'v_rad': v_rad_param_name
#             }
    
#     # Add parameters for spectrum 1
#     for species_idx in species_params1:
#         species_data = species_params1[species_idx]
        
#         # Transition parameters
#         n_trans = species_data['lambda'].size
#         params.add(f'spec1_n_trans_{species_idx}', value=n_trans, vary=False)
        
#         for i in range(n_trans):
#             params.add(f'spec1_lambda_{species_idx}_{i}', value=species_data['lambda'][i], vary=False)
#             params.add(f'spec1_f_{species_idx}_{i}', value=species_data['f'][i], vary=False)
#             params.add(f'spec1_gamma_{species_idx}_{i}', value=species_data['gamma'][i], vary=False)
        
#         # Component parameters
#         n_component = species_data['v_rad'].size
#         params.add(f'spec1_n_component_{species_idx}', value=n_component, vary=False)
        
#         # Link component parameters to master parameters
#         for i in range(n_component):
#             params.add(f'spec1_b_{species_idx}_{i}', expr=master_params[(species_idx, i)]['b'])
#             params.add(f'spec1_N_{species_idx}_{i}', expr=master_params[(species_idx, i)]['N'])
#             params.add(f'spec1_v_rad_{species_idx}_{i}', expr=master_params[(species_idx, i)]['v_rad'])
    
#     # Add parameters for spectrum 2
#     for species_idx in species_params2:
#         species_data = species_params2[species_idx]
        
#         # Transition parameters
#         n_trans = species_data['lambda'].size
#         params.add(f'spec2_n_trans_{species_idx}', value=n_trans, vary=False)
        
#         for i in range(n_trans):
#             params.add(f'spec2_lambda_{species_idx}_{i}', value=species_data['lambda'][i], vary=False)
#             params.add(f'spec2_f_{species_idx}_{i}', value=species_data['f'][i], vary=False)
#             params.add(f'spec2_gamma_{species_idx}_{i}', value=species_data['gamma'][i], vary=False)
        
#         # Component parameters  
#         n_component = species_data['v_rad'].size
#         params.add(f'spec2_n_component_{species_idx}', value=n_component, vary=False)
        
#         # Link component parameters to master parameters (same as spectrum 1)
#         for i in range(n_component):
#             params.add(f'spec2_b_{species_idx}_{i}', expr=master_params[(species_idx, i)]['b'])
#             params.add(f'spec2_N_{species_idx}_{i}', expr=master_params[(species_idx, i)]['N'])
#             params.add(f'spec2_v_rad_{species_idx}_{i}', expr=master_params[(species_idx, i)]['v_rad'])
    
#     # Create a combined model function
#     def dual_voigt_model(**params_dict):
#         """Combined model that computes both spectra and concatenates results"""
        
#         # Extract parameters for spectrum 1
#         spec1_params = {}
#         spec2_params = {}
        
#         # Global parameters
#         spec1_params['n_species'] = params_dict['n_species']
#         spec1_params['v_resolution'] = params_dict['v_resolution'] 
#         spec1_params['n_step'] = params_dict['n_step']
#         spec1_params['wavegrid'] = wavegrid1
        
#         spec2_params['n_species'] = params_dict['n_species']
#         spec2_params['v_resolution'] = params_dict['v_resolution']
#         spec2_params['n_step'] = params_dict['n_step'] 
#         spec2_params['wavegrid'] = wavegrid2
        
#         # Copy spectrum-specific parameters
#         for key, value in params_dict.items():
#             if key.startswith('spec1_'):
#                 new_key = key.replace('spec1_', '')
#                 spec1_params[new_key] = value
#             elif key.startswith('spec2_'):
#                 new_key = key.replace('spec2_', '')
#                 spec2_params[new_key] = value
        
#         # Compute both models
#         model1 = Voigt_fit_wrapper(**spec1_params)
#         model2 = Voigt_fit_wrapper(**spec2_params)
        
#         # Return concatenated results
#         return np.concatenate([model1, model2])
    
#     # Concatenate the observed data
#     combined_ydata = np.concatenate([ydata1, ydata2])
    
#     # Create weights (inverse of standard deviations)
#     weights1 = np.full_like(ydata1, 1/std_dev1)
#     weights2 = np.full_like(ydata2, 1/std_dev2)
#     combined_weights = np.concatenate([weights1, weights2])
    
#     # Create and fit the model
#     dual_model = Model(dual_voigt_model)
#     result = dual_model.fit(combined_ydata, params, weights=combined_weights)
    
#     return result


# --- Simultaneous Continuum and Voigt Fitting ---

def continuum_voigt_wrapper(**params_list):
    """
    Wrapper model that multiplies a Cubic Spline continuum with the Voigt absorption model.
    
    Physics:
        Model(lambda) = Continuum(lambda) * Absorption(lambda)
        
        Continuum: Cubic Spline defined by knots (knot_x, knot_y).
        Absorption: Voigt profile sum (e^(-tau)).
        
    Expected additional parameters in params_list:
    - knot_x_array: Array of wavelength positions for the knots (passed as independent var).
    - knot_y_{i}: Flux values at each knot position.
    """
    # 1. Calculate Absorption Component (normalized flux)
    # We pass all params to Voigt_fit_wrapper. It only uses what it needs.
    absorption_model = Voigt_fit_wrapper(**params_list)
    
    # 2. Calculate Continuum Component
    wavegrid = params_list['wavegrid']
    knot_x_array = params_list['knot_x_array']
    
    # Extract knot_y values from parameters
    # The number of knots is determined by the length of knot_x_array
    n_knots = len(knot_x_array)
    knot_y_values = []
    
    # Get spline order (default to 3/Cubic)
    spline_order = params_list.get('spline_order', 3)
    # Ensure it's a scalar integer (sometimes passed as 0-d array by lmfit)
    if hasattr(spline_order, 'shape') and spline_order.shape == ():
        spline_order = int(spline_order)
    
    for i in range(n_knots):
        # We expect params named 'knot_y_0', 'knot_y_1', etc.
        knot_y_values.append(params_list[f'knot_y_{i}'])
        
    # Create Spline for continuum
    if spline_order == 3:
        # User standard CubicSpline for k=3 (preserves existing behavior)
        spline = CubicSpline(knot_x_array, knot_y_values)
    else:
        # Use make_interp_spline for Linear (k=1) or Quadratic (k=2)
        from scipy.interpolate import make_interp_spline
        spline = make_interp_spline(knot_x_array, knot_y_values, k=spline_order)
        
    continuum_model = spline(wavegrid)
    
    # 3. Combine
    return continuum_model * absorption_model


def _generate_smart_knots(wavegrid, species_params, n_knots, avoidance_width=0.5, absorption_ranges=None):
    """
    Helper to generate continuum knots that avoid absorption line centers.
    It accounts for the Radial Velocity (v_rad) shift using the initial guesses.
    
    Parameters:
    - avoidance_width: Minimum distance (in Angstroms) a knot should be from a line center.
    - absorption_ranges: List of tuples [(start, end), ...] defining regions to strictly avoid.
    """
    w_min, w_max = np.min(wavegrid), np.max(wavegrid)
    initial_knots = np.linspace(w_min, w_max, n_knots)
    
    # 1. Handle Explicit Absorption Ranges (User Provided)
    if absorption_ranges is not None:
        # Ensure it's a list of tuples
        if isinstance(absorption_ranges, tuple):
            absorption_ranges = [absorption_ranges]
            
        final_knots = []
        for k in initial_knots:
            is_bad = False
            for (start, end) in absorption_ranges:
                if start <= k <= end:
                    is_bad = True
                    break
            if not is_bad:
                final_knots.append(k)
        
        # Add anchor points just outside the absorption ranges
        for (start, end) in absorption_ranges:
            # Add knots slightly outside the range to pin the spline
            # Check if they are within grid bounds
            if start - avoidance_width >= w_min:
                final_knots.append(start - avoidance_width)
            if end + avoidance_width <= w_max:
                final_knots.append(end + avoidance_width)
                
        return np.sort(np.unique(final_knots))

    # 2. Fallback to Automatic Detection (if no ranges provided)
    c_light = 299792.458  # km/s
    
    # Collect all shifted line centers
    shifted_line_centers = []
    
    for species_idx in species_params:
        sp = species_params[species_idx]
        
        # Get Wavelengths (Transitions)
        lambdas = sp.get('lambda', [])
        if np.isscalar(lambdas): lambdas = [lambdas]
        else: lambdas = np.asarray(lambdas).flatten()
        
        # Get Radial Velocities (Components)
        v_rads = sp.get('v_rad', [])
        if np.isscalar(v_rads): v_rads = [v_rads]
        else: v_rads = np.asarray(v_rads).flatten()
                 
        # Generate all observed center wavelengths
        # Each transition produces a line at each component velocity
        for lam in lambdas:
            for v in v_rads:
                # Apply Doppler shift: obs = rest * (1 + v/c)
                shifted_lam = lam * (1 + v / c_light)
                
                # Only care if it's within our grid range (plus buffer)
                if w_min - 2.0 < shifted_lam < w_max + 2.0:
                    shifted_line_centers.append(shifted_lam)
            
    shifted_line_centers = np.array(shifted_line_centers)
    
    final_knots = []
    for k in initial_knots:
        # Check distance to all lines
        if len(shifted_line_centers) > 0:
            dists = np.abs(shifted_line_centers - k)
            collision_idx = np.where(dists < avoidance_width)[0]
            
            if len(collision_idx) > 0:
                # Collision detected!
                # Move knot to the right edge of the avoidance zone
                new_k = k + avoidance_width
                if new_k > w_max:
                    new_k = k - avoidance_width # Try left if right is OOB
                k = new_k
                
        final_knots.append(k)
        
    return np.sort(np.unique(final_knots))


def astro_simultaneous_fit(
    wavegrid, 
    ydata, 
    species_params,
    n_knots=10,
    knots_x_array=None,
    v_resolution=0.0, 
    n_step=25, 
    std_dev=0.002,
    avoidance_width=0,
    absorption_ranges=None,
    spline_order=3
):
    """
    Perform a simultaneous fit of the continuum (Spline) and absorption (Voigt).
    
    This avoids the error propagation issues of normalizing first and then fitting.
    
    Parameters:
    -----------
    ...
    spline_order : int
        Order of the spline for continuum fitting. 1=Linear, 2=Quadratic, 3=Cubic.
    wavegrid : array
        Wavelength grid
    ydata : array
        Observed flux (NOT normalized)
    species_params : dict
        Parameters for the Voigt species (same format as astro_voigt_fit)
    n_knots : int
        Number of spline knots to generate (equally spaced) if knots_x_array is None.
    knots_x_array : array, optional
        Explicit wavelength positions for the continuum knots.
    v_resolution : float
        Velocity resolution
    n_step : int
        Number of steps for Voigt calculation
    std_dev : float
        Standard deviation for weighting
    avoidance_width : float
        Distance (in A) to shift knots away from line centers.
    absorption_ranges : list of tuples, optional
        List of (start, end) tuples defining regions to strictly avoid placing knots.
        e.g. [(6707.5, 6708.2)]
        
    Returns:
    --------
    result : lmfit.model.ModelResult
        Fitting result
    """
    
    # 1. Setup Continuum Parameters
    if knots_x_array is None:
        # Smart generation
        knots_x_array = _generate_smart_knots(
            wavegrid, species_params, n_knots, avoidance_width, absorption_ranges
        )
        n_knots = len(knots_x_array)
    else:
        knots_x_array = np.asarray(knots_x_array)
        n_knots = len(knots_x_array)
        
    # 2. Setup Voigt Parameters (Reuse logic from astro_voigt_fit)
    n_species = len(species_params)
    
    # Ensure numpy arrays
    for species_idx in species_params:
        for key in ['lambda', 'f', 'gamma', 'b', 'N', 'v_rad']:
            species_params[species_idx][key] = np.asarray(species_params[species_idx][key])
            
    params = Parameters()
    params.add('n_species', value=n_species, vary=False)
    params.add('v_resolution', value=v_resolution, vary=False)
    params.add('n_step', value=n_step, vary=False)
    
    # Add Voigt Species Parameters
    for species_idx in range(n_species):
        species_data = species_params[species_idx]
        
        # Transition parameters
        n_trans = species_data['lambda'].size
        params.add(f'n_trans_{species_idx}', value=n_trans, vary=False)
        
        for i in range(n_trans):
            params.add(f'lambda_{species_idx}_{i}', value=species_data['lambda'][i], vary=False)
            params.add(f'f_{species_idx}_{i}', value=species_data['f'][i], vary=False)
            params.add(f'gamma_{species_idx}_{i}', value=species_data['gamma'][i], vary=False)
        
        # Component parameters
        n_component = species_data['v_rad'].size
        params.add(f'n_component_{species_idx}', value=n_component, vary=False)
        
        tie_b = species_data.get('tie_b', False) if species_idx > 0 else False
        tie_v_rad = species_data.get('tie_v_rad', False) if species_idx > 0 else False
        
        for i in range(n_component):
            # b parameter
            if tie_b:
                params.add(f'b_{species_idx}_{i}', expr=f'b_0_{i}')
            else:
                params.add(f'b_{species_idx}_{i}', value=species_data['b'][i], min=0.2, max=5.5, vary=True)

            # N parameter
            params.add(f'N_{species_idx}_{i}', value=species_data['N'][i], min=0, vary=True)

            # v_rad parameter
            if tie_v_rad:
                params.add(f'v_rad_{species_idx}_{i}', expr=f'v_rad_0_{i}')
            else:
                params.add(f'v_rad_{species_idx}_{i}', value=species_data['v_rad'][i], vary=True)

    # 3. Add Continuum Spline Parameters
    # Guess initial y-values from data
    # (Use interp to get values at knot positions)
    if len(ydata) == len(wavegrid):
        initial_knot_y = np.interp(knots_x_array, wavegrid, ydata)
    else:
        # Fallback if ydata/wavegrid mismatch (shouldn't happen in fit)
        initial_knot_y = np.ones(n_knots)
    
    for i in range(n_knots):
        # Allow knots to vary. min=0 usually safe for flux.
        params.add(f'knot_y_{i}', value=initial_knot_y[i], min=0, vary=True)

    # 4. Create Model
    # We declare knot_x_array and spline_order as independent so they are passed to the wrapper
    model = Model(continuum_voigt_wrapper, independent_vars=['wavegrid', 'knot_x_array', 'spline_order'])
    
    # 5. Fit
    result = model.fit(ydata, params, wavegrid=wavegrid, knot_x_array=knots_x_array, spline_order=spline_order, weights=1/std_dev)
    
    return result

