import numpy as np
from pysynphot import observation
from pysynphot import spectrum

def rebin_spec(wave, specin, wavnew):

    import numpy as np
    from pysynphot import observation
    from pysynphot import spectrum 
   
    spec = spectrum.ArraySourceSpectrum(wave=wave, flux=specin)
    f = np.ones(len(wave))
    filt = spectrum.ArraySpectralElement(wave, f, waveunits='angstrom')
    obs = observation.Observation(spec, filt, binset=wavnew, force='taper')

    return obs.binflux
