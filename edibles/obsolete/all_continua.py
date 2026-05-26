import numpy as np
import pandas as pd

from edibles.utils.edibles_spectrum import EdiblesSpectrum
from edibles.continuum import Continuum


obslog = pd.read_csv('edibles/data/DR4_ObsLog.csv')

for i in range(len(obslog)):

    obs = obslog.iloc[i]

    if obs['Order'] == 'ALL':
        print("Skipping merged file")
        continue

    centwave = obs['Setting']
    if centwave == 346:
        # trim
        tr_min = 0.20
        tr_max = 0.12
    elif centwave == 437:
        # trim
        tr_min = 0.08
        tr_max = 0.05
        if int(obs['Order']) - 1 > 20:
            tr_min = 0.04
            tr_max = 0.02
    elif centwave == 564:
        # trim
        if 'redl' in obs['Filename']:
            tr_min = 0.28
            tr_max = 0.12
        if 'redu' in obs['Filename']:
            tr_min = 0.2
            tr_max = 0.08
    elif centwave == 860:
        # trim
        if 'redl' in obs['Filename']:
            tr_min_list = [0.35, 0.1, 0.01, 0.005, 0.005, 0.005, 0.002, 0.002, 0.002, 0.002,
                           0.004, 0, 0, 0.03, 0.0, 0.0, 0.00, 0.0, 0, 0, 0.0]
            tr_min = tr_min_list[int(obs['Order']) - 1]
            tr_max_list = [0.15, 0.1, 0.01, 0.075, 0.05, 0.045, 0.002, 0.002, 0.002, 0.002,
                           0.002, 0, 0, 0.03, 0.1, 0.1, 0.12, 0.0, 0.15, 0.3, 0.1]
            tr_max = tr_max_list[int(obs['Order']) - 1]
        if 'redu' in obs['Filename']:
            tr_min_list = [0, 0.12, 0.130, 0.12, 0.10, 0.07, 0.045, 0.03,
                           0.040, 0.15, 0.58, 0.77, 0.87, 0.87]
            tr_min = tr_min_list[int(obs['Order']) - 1]
            tr_max_list = [0, 0.005, 0.025, 0.02, 0.02, 0.02, 0.02, 0.01, 0.025, 0.02, 0, 0, 0, 0]
            tr_max = tr_max_list[int(obs['Order']) - 1]

    # determine the order numbers
    if centwave == 346:
        order_number = [34]
    if centwave == 437:
        order_number = [32]
    if centwave == 564:
        order_number = [25, 17]
    if centwave == 860:
        order_number = [21, 14]

    wavemin = float(obs['WaveMin'])
    wavemax = float(obs['WaveMax'])
    spectrum_length = wavemax - wavemin
    min_cutoff = spectrum_length * tr_min
    max_cutoff = spectrum_length * tr_max
    wavemin = wavemin + min_cutoff
    wavemax = wavemax - max_cutoff

    sp = EdiblesSpectrum(obs['Filename'])

    idx = np.where(np.logical_and(sp.wave > wavemin, sp.wave < wavemax))
    sp.wave = sp.wave[idx]
    sp.flux = sp.flux[idx]

    # build a 4 anchor points spline
    cont = Continuum(sp, method="spline", n_anchors=4, plot=False, verbose=0)
    # Guess the model parameters
    params = cont.model.guess(sp.flux, x=sp.wave)
    # Fit the model
    result = cont.model.fit(data=sp.flux, params=params, x=sp.wave)
    # Get the output of the fit model
    out = result.eval(params=result.params, x=sp.wave)

    resid = sp.flux - out
    std = np.std(resid)

    idx2 = np.where(sp.flux >= (out - std))
    result2 = cont.model.fit(data=sp.flux[idx2], params=params, x=sp.wave[idx2])
    # Get the output of the fit model
    out2 = result2.eval(params=result2.params, x=sp.wave)

    cont.add_to_csv(
        user="Klay Kulik", comments="Initial fit of order"
    )

    if i % 10 == 0:
        print('Completed {} of {} splines'.format(i, len(obslog)))

print('Completed all {} splines'.format(len(obslog)))
