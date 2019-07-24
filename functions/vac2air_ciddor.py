def vac2air_ciddor(vacw):
    """ Convert vacuum wavelengths in Angstroms to air wavelengths.

    This uses the relation from Ciddor 1996, Applied Optics LP,
    vol. 35, Issue 9, p.1566. Only valid for wavelengths > 2000 Ang.
    """
    k0 = 238.0185
    k1 = 1e-8 * 5792105.
    k2 = 57.362
    k3 = 1e-8 * 167917.
    s2 = (1e4 / vacw)**2
    n = 1 + k1/(k0 - s2) + k3/(k2 - s2)
    airw = vacw / n

    return airw

""" example:

wt,ft=np.loadtxt("transmission.dat", unpack=True)
wtc = vac2air_ciddor(wt*10.0) #make sure wave is in \AA.

"""
