# This script runs through all your fitting results for every sightline completed by voigt_gui and creates a table
# with the fit results. 


import os
from pathlib import Path
import pandas as pd
import numpy as np

base_dir = Path(r"C:\Users\User\Downloads\PYTHON PROJECTS\edibles\edibles\data\voigt_fitting_data\TiII fits")
output_csv_path = base_dir / "fit_results_TiII.csv"

all_rows = []

def find_csv(folder):
    """Find the fit-results CSV in a folder, regardless of naming convention."""
    candidates = list(folder.glob("*.csv"))
    candidates = [c for c in candidates if c.name != output_csv_path.name]
    if not candidates:
        return None
    for c in candidates:
        if 'error' in c.name.lower():
            return c
    return candidates[0]

def compute_chisqr_redchi(folder, csv_file, df_csv):
    stem = csv_file.stem
    base_stem = stem.replace('_errors', '')

    spec_path = folder / f"{base_stem}.dat"
    bestfit_path = folder / f"{base_stem}_bestfit.dat"

    if not spec_path.exists() or not bestfit_path.exists():
        return None, None

    try:
        spec = np.genfromtxt(spec_path, unpack=True)
        bestfit = np.genfromtxt(bestfit_path)

        flux = spec[1]
        err = spec[2]

        if len(bestfit) != len(flux):
            print(f"  Warning: length mismatch between spectrum and bestfit for {folder.name}")
            return None, None

        residuals = (flux - bestfit) / err
        valid = np.isfinite(residuals)
        residuals = residuals[valid]

        chisqr = np.sum(residuals**2)

        n_v_params = df_csv['v_comp'].nunique() * 2
        n_n_params = df_csv['n_comp'].nunique()
        n_free_params = n_v_params + n_n_params

        n_points = len(residuals)
        dof = n_points - n_free_params
        redchi = chisqr / dof if dof > 0 else None

        return chisqr, redchi

    except Exception as e:
        print(f"  Error computing chisqr for {folder.name}: {e}")
        return None, None

for folder in base_dir.iterdir():
    if not folder.is_dir():
        continue

    sightline_name = folder.name

    csv_file = find_csv(folder)
    if csv_file is None:
        print(f"Warning: No CSV file found in {folder.name}. Skipping...")
        continue

    try:
        df_csv = pd.read_csv(csv_file)

        sav_path = csv_file.with_suffix('.sav')
        base_sav_path = folder / f"{csv_file.stem.replace('_errors', '')}.sav"

        chisqr_val = None
        redchi_val = None

        if sav_path.exists():
            with open(sav_path, 'r') as f:
                for line in f:
                    if line.startswith('chisqr:'):
                        chisqr_val = float(line.split(':')[1].strip())
                    elif line.startswith('redchi:'):
                        redchi_val = float(line.split(':')[1].strip())
        elif base_sav_path.exists():
            with open(base_sav_path, 'r') as f:
                for line in f:
                    if line.startswith('chisqr:'):
                        chisqr_val = float(line.split(':')[1].strip())
                    elif line.startswith('redchi:'):
                        redchi_val = float(line.split(':')[1].strip())
        else:
            chisqr_val, redchi_val = compute_chisqr_redchi(folder, csv_file, df_csv)
            if chisqr_val is not None:
                new_sav_path = folder / f"{csv_file.stem.replace('_errors', '')}.sav"
                with open(new_sav_path, 'w') as f:
                    f.write(f'chisqr: {chisqr_val}\n')
                    f.write(f'redchi: {redchi_val}\n')
                print(f"  Created sav file: {new_sav_path.name}")
            else:
                print(f"  Warning: Could not compute chisqr/redchi for {folder.name} (missing .dat or bestfit.dat)")

        comps = df_csv[['v_comp']].drop_duplicates().sort_values('v_comp')
        n_clouds = len(comps)

        first_row_for_sightline = True

        for _, comp_row in comps.iterrows():
            v_comp = comp_row['v_comp']
            sub_df = df_csv[df_csv['v_comp'] == v_comp].iloc[0]

            b_fit = sub_df.get('b_fit')
            b_err = sub_df.get('b_err')
            n_fit = sub_df.get('n_fit')
            n_err = sub_df.get('n_err')
            v_fit = sub_df.get('v_rad_fit')
            v_err = sub_df.get('v_rad_err')

            b_pct = abs(b_err / b_fit * 100) if pd.notna(b_err) and b_fit not in (0, None) else None
            n_pct = abs(n_err / n_fit * 100) if pd.notna(n_err) and n_fit not in (0, None) else None
            v_pct = abs(v_err / v_fit * 100) if pd.notna(v_err) and v_fit not in (0, None) else None

            b_init = sub_df.get('b_init')
            n_init = sub_df.get('n_init')
            v_init = sub_df.get('v_rad_init')
            init_str = f"({b_init}, {n_init}, {v_init})"

            row = {
                'Sightline': sightline_name if first_row_for_sightline else '',
                'No. of cloud': n_clouds if first_row_for_sightline else '',
                'Width (b)': b_fit,
                'Width error (+/-)': b_err,
                'Width error (%)': b_pct,
                'Number Density (N)': n_fit,
                'Number Density error (+/-)': n_err,
                'Number Density (%)': n_pct,
                'Radial Velocity (V)': v_fit,
                'Radial Velocity error (+/-)': v_err,
                'Radial Velocity error (%)': v_pct,
                'Initial values (b, N, V)': init_str,
                'chisqr': chisqr_val if first_row_for_sightline else '',
                'redchi': redchi_val if first_row_for_sightline else '',
                'Comment': '',
            }
            all_rows.append(row)
            first_row_for_sightline = False

        print(f"Processed: {sightline_name}")

    except Exception as e:
        print(f"Error processing {folder.name}: {e}")

if all_rows:
    final_df = pd.DataFrame(all_rows)
    final_df.to_csv(output_csv_path, index=False)
    print(f"\nSuccess! Master file created at: {output_csv_path}")
else:
    print("\nNo data pairs found to compile.")