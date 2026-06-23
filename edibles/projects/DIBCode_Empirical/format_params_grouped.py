def format_params_grouped(params):
    """Formats parameters grouped by component type for structured display.

    Extracts and sequences parameter values to clearly separate the continuum
    (Chebyshev coefficients), Gaussian components, and Lorentzian components.

    Args:
        params: An lmfit Parameters object containing the model coefficients.

    Returns:
        A multiline formatted string separating the parameters by type.
    """
    lines = []

    # Continuum (Chebyshev)
    cheb_params = []
    i = 0
    while f'c{i}' in params:
        val = params[f'c{i}'].value
        cheb_params.append(f'c{i}={val:.4f}')
        i += 1
    if cheb_params:
        lines.append('Continuum: ' + ', '.join(cheb_params))

    # Gaussians
    gauss_list = []
    i = 0
    while f'g{i}_amp' in params:
        amp = params[f'g{i}_amp'].value
        cen = params[f'g{i}_cen'].value
        wid = params[f'g{i}_wid'].value
        gauss_list.append(f'g{i}[{amp:.8f}, {cen:.8f}, {wid:.8f}]')
        i += 1
    if gauss_list:
        lines.append('Gaussians: ' + ', '.join(gauss_list))

    # Lorentzians
    lorentz_list = []
    i = 0
    while f'l{i}_amp' in params:
        amp = params[f'l{i}_amp'].value
        cen = params[f'l{i}_cen'].value
        wid = params[f'l{i}_wid'].value
        lorentz_list.append(f'l{i}[{amp:.8f}, {cen:.8f}, {wid:.8f}]')
        i += 1
    if lorentz_list:
        lines.append('Lorentzians: ' + ', '.join(lorentz_list))

    return '\n'.join(lines)
