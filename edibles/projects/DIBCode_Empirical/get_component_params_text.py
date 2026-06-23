def get_component_params_text(params, cand_type, cand_idx):
    """Extracts and formats specific model parameters for targeted display.

    Supports extracting Gaussian peaks, Lorentzian peaks, or Chebyshev
    coefficients.

    Args:
        params: An lmfit Parameters object containing the model coefficients.
        cand_type: The type of component ('gaussian', 'lorentzian', or
          'chebyshev').
        cand_idx: The index identifier of the candidate component (e.g., 0, 1).

    Returns:
        A formatted string summarizing the parameters of the chosen component.
    """
    if cand_type == 'gaussian':
        prefix = f'g{cand_idx}_'
        if f'{prefix}amp' in params:
            amp = params[f'{prefix}amp'].value
            cen = params[f'{prefix}cen'].value
            wid = params[f'{prefix}wid'].value
            return f"Amp={amp:.2f}, Cen={cen:.8f}, Wid={wid:.8f}"
    elif cand_type == 'lorentzian':
        prefix = f'l{cand_idx}_'
        if f'{prefix}amp' in params:
            amp = params[f'{prefix}amp'].value
            cen = params[f'{prefix}cen'].value
            wid = params[f'{prefix}wid'].value
            return f"Amp={amp:.2f}, Cen={cen:.8f}, Wid={wid:.8f}"
    elif cand_type == 'chebyshev':
        if f'c{cand_idx}' in params:
            val = params[f'c{cand_idx}'].value
            return f"c{cand_idx}={val:.4f}"
    return ""
