from edibles.projects.global_voigt_fit.voigt_fitting import voigt_fit_wrapper
from pathlib import Path
import pandas as pd
import numpy as np

#possible button for the gui?

# to use this, organize your fit results from the voigt_gui.py into folders that are named by target.
# e.g: HD 23180 (make sure to include a space in between 'HD' and the numbers)

# keep all of these folders in one main folder, and paste that path in the base_dir = Path(...) down below

# keep in mind that this will go through ALL of the target folders in your main folder at once.

base_dir = Path(r"C:\Users\User\Downloads\PYTHON PROJECTS\edibles\edibles\data\voigt_fitting_data") 

for star_dir in base_dir.iterdir():
    for sav_file in star_dir.glob('*.sav'):
        csv_file = sav_file.with_suffix('.csv')
        dat_file = sav_file.with_suffix('.dat')

        fit_df = pd.read_csv(csv_file)
        fit_spec = np.genfromtxt(dat_file, unpack=True)

        fit_df['v_rad_init'] = fit_df['v_rad_fit']
        fit_df['b_init'] = fit_df['b_fit']
        fit_df['n_init'] = fit_df['n_fit']
        fit_df['weight'] = 1.0
        #fit_df.loc[fit_df['n_fit'] == 1e9, 'n_init'] = 0.0 (only use this for targets where you mask out a line)

        result = voigt_fit_wrapper(fit_df, fit_spec)

        for name, par in result.params.items():
            parts = name.rsplit('_', 1)
            if len(parts) != 2:
                continue
            param_type, comp_idx = parts
            try:
                comp_idx = int(comp_idx)
            except ValueError:
                continue
            if param_type == 'v_rad':
                fit_df.loc[fit_df['v_comp'] == comp_idx, 'v_rad_fit'] = par.value
                fit_df.loc[fit_df['v_comp'] == comp_idx, 'v_rad_err'] = par.stderr
            elif param_type == 'b':
                fit_df.loc[fit_df['b_comp'] == comp_idx, 'b_fit'] = par.value
                fit_df.loc[fit_df['b_comp'] == comp_idx, 'b_err'] = par.stderr
            elif param_type == 'n':
                fit_df.loc[fit_df['n_comp'] == comp_idx, 'n_fit'] = par.value
                fit_df.loc[fit_df['n_comp'] == comp_idx, 'n_err'] = par.stderr

        fit_df.to_csv(star_dir / (sav_file.stem + '_errors.csv'), index=False)
        print(f'Done: {star_dir.name} / {sav_file.stem}')