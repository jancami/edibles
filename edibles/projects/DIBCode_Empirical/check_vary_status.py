def check_vary_status(params):
    """Checks and reports which parameters are configured to vary during a fit.

    Scans through the parameters and separates them based on their boolean `vary` property.

    Args:
        params: An lmfit Parameters object containing the model coefficients.

    Returns:
        A tuple containing two lists:
            - A list of parameter names set to vary (True).
            - A list of parameter names set to fixed (False).
    """
    vary_true = []
    vary_false = []
    for name, param in params.items():
        if param.vary:
            vary_true.append(name)
        else:
            vary_false.append(name)
    return vary_true, vary_false
