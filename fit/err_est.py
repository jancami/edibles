import numpy as np


def errorEst(wave, flux, fit_flag=None):

    if fit_flag is None:
        fit_flag = ''

    wave = np.array(wave)
    flux = np.array(flux)

    # find a flat region based on heighest S/N
    first_idx = 0
    SNR_f = 0
    Nw = len(wave)
    if Nw > 50: Nd = 5
    else: Nd = 5

    for loop_sn in range(Nd):

        # split region into sub-regions
        if fit_flag != 'multi_trans':
            if loop_sn != Nd-1:
                sub_wave = wave[loop_sn*Nd:(loop_sn+1)*Nd]
                sub_flux = flux[loop_sn*Nd:(loop_sn+1)*Nd]
            else:
                sub_wave = wave[loop_sn*Nd:]
                sub_flux = flux[loop_sn*Nd:]
            first_idx += len(sub_wave)
        else:
            idx = np.nonzero(wave < min(wave)+4)
            sub_wave = wave[idx]
            sub_flux = flux[idx]


        # compute its S/N
        SNR = np.mean(sub_flux)/np.std(sub_flux)
        if SNR >= SNR_f:
            hist = [x - np.mean(sub_flux) for x in sub_flux]
            SNR_f = SNR

    # error
    SD = np.std(hist)
    SE = SD/np.sqrt(len(hist))

    return np.full(len(wave), SD)
