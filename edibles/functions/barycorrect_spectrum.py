import astropy.constants as cst


def barycorrectSpectrum(wave_array, v_bary):
    """Barycentric wavelength correction tool.

    INPUT:
        wave_array: [ndarray]   The wavelength array of the spectrum
        v_bary:     [float]     The barycentric velocity
                                - probably from EdiblesSpectrum
    OUTPUT:
        bary_wave:  [ndarray] The corrected wavelength array
    """
    wave_array = wave_array + (v_bary / cst.c.to("km/s").value) * wave_array

    return wave_array
