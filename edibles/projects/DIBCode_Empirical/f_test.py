from scipy import stats
def f_test(chi2_simple, chi2_complex, df_simple, df_complex):
    """Performs an F-test for nested models to compare fit improvement.

    Null Hypothesis (H0): The complex model does not significantly improve the
    fit compared to the simple model.

    Args:
        chi2_simple: Chi-square value from the simpler (fewer parameters) model.
        chi2_complex: Chi-square value from the complex (more parameters)
          model.
        df_simple: Degrees of freedom for the simpler model.
        df_complex: Degrees of freedom for the complex model.

    Returns:
        A tuple containing:
            - The calculated F-statistic as a float (0.0 if no improvement).
            - The resulting p-value as a float (1.0 if no improvement).
    """
    delta_chi2 = chi2_simple - chi2_complex
    delta_df = df_simple - df_complex

    if delta_df <= 0 or delta_chi2 <= 0:
        return 0, 1.0  # No improvement

    F = (delta_chi2 / delta_df) / (chi2_complex / df_complex)
    p_value = 1 - stats.f.cdf(F, delta_df, df_complex)

    return F, p_value
