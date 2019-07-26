def read_sky_transmission(transmission_dat, scale_factor=1.0):

    import vac2air_ciddor
    
    """ transmission input spectrum in nm
    """
    wt,ft = np.loadtxt(transmission_dat, unpack=True)
    wtc = vac2air_ciddor(wt*10.0)
    ft = ft**scale_factor
    return, wtc, ft
