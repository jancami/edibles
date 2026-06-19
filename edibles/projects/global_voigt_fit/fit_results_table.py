import os
from pathlib import Path
import pandas as pd

base_dir = Path(r"C:\Users\User\Downloads\PYTHON PROJECTS\edibles\edibles\data\voigt_fitting_data\TiII fits")
output_csv_path = base_dir / "fit_results.csv"

all_rows = []

for folder in base_dir.iterdir():
    if not folder.is_dir():
        continue

    sightline_name = folder.name

    csv_files = list(folder.glob("*.csv"))
    if not csv_files:
        print(f"Warning: No CSV file found in {folder.name}. Skipping...")
        continue

    csv_file = csv_files[0]
    sav_path = csv_file.with_suffix('.sav')

    if not sav_path.exists():
        print(f"Warning: Missing SAV file for {folder.name}. Skipping...")
        continue

    try:
        df_csv = pd.read_csv(csv_file)

        # Get one row per unique velocity component
        comps = df_csv[['v_comp']].drop_duplicates().sort_values('v_comp')
        n_clouds = len(comps)

        chisqr_val = None
        redchi_val = None
        with open(sav_path, 'r') as f:
            for line in f:
                if line.startswith('chisqr:'):
                    chisqr_val = float(line.split(':')[1].strip())
                elif line.startswith('redchi:'):
                    redchi_val = float(line.split(':')[1].strip())

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