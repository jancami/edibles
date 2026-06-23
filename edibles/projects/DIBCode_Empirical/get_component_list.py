import numpy as np
def get_component_list(x, params):
    """
    Uses a dictionary of model parameters to build a list of separate profile arrays for the model.
    Args:
        x : Independent variable data points (wavelengths) as a NumPy array.
        params: Parameters from the model (Gaussian, Lorentizan, Chebyeshev polynomials)
    Returns:
        profiles: List of profiles.
        fids: Sorted list of parameters.

    """
    profiles = []

    relevant_keys = [k for k in params.keys() if not k.startswith('c')]
    fids = sorted(list(set(k.split('_')[0] for k in relevant_keys if '_amp' in k)))

    for fid in fids:
        a = params[f'{fid}_amp']
        c = params[f'{fid}_cen']
        w = params[f'{fid}_wid']

        if fid.startswith('g'):
            profiles.append(a * np.exp(-((x - c) / w) ** 2))
        elif fid.startswith('l'):
            # Using the same Lorentzian math from the composite_model
            profiles.append(a * w ** 2 / ((x - c) ** 2 + w ** 2))

    return profiles, fids
