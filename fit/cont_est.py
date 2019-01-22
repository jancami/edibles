import numpy as np

def cont_est(wave, spec, cw):

    try:
        if len(cw[0]) >= 1:
            xcl = cw[0][0]
            xcr = cw[0][-1]
    except TypeError:
        if len(cw) == 1:
            xcl = np.mean(cw[0])
            xcr = np.mean(cw[0])
        if len(cw) > 1:
            xcl = np.mean(cw[0])
            xcr = np.mean(cw[-1])

    idxl = np.nonzero((wave > xcl-2) & (wave < xcl-0.5))
    cnxl = wave[idxl]
    cnyl = spec[idxl]
    idxr = np.nonzero((wave > xcr+0.5) & (wave < xcr+2))
    cnxr = wave[idxr]
    cnyr = spec[idxr]

    # fit a straight line
    cont_x = np.concatenate((np.array(cnxl), np.array(cnxr)), axis=0)
    cont_y = np.concatenate((np.array(cnyl), np.array(cnyr)), axis=0)

    mean = np.mean(cont_y)
    weights = cont_y / mean

    pol = np.poly1d(np.polyfit(cont_x, cont_y, 2, w=weights))
    continuum_lim = np.mean(pol(cont_x))

    return continuum_lim
