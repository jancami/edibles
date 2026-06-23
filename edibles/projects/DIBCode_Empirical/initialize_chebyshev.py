def initialize_chebyshev(order, prefix='c'):
    """Initializes a specific order Chebyshev coefficient to zero.

    Args:
        order: The order/index of the Chebyshev polynomial coefficient (e.g., 0, 1, 2).
        prefix: String identifier prefix to prepend to the coefficient name.

    Returns:
        A dictionary mapping the named coefficient to its initialized value of zero.
    """
    return {f'{prefix}{order}': {'value': 0.0}}
