from edibles.utils.voigt_fitting import voigt_fit_wrapper
import tkinter as tk
from importlib.resources import files
import pandas as pd

atomic_line_file = files('edibles') / 'data/auxiliary_data/line_catalogs/edibles_linelist_atoms.csv'
atomic_line_list = pd.read_csv(atomic_line_file)
atomic_line_list = atomic_line_list.dropna(subset=['Gamma'])

def make_default_df(in_df: pd.DataFrame, v_comp: int):
    i_df = in_df.copy()
    i_df.loc[:, 'v_comp'] = v_comp
    i_df.loc[:, f'v_rad_init'] = 0
    i_df.loc[:, f'v_rad_min'] = -100
    i_df.loc[:, f'v_rad_max'] = 100
    i_df.loc[:, 'b_comp'] = v_comp
    i_df.loc[:, f'b_init'] = 2
    i_df.loc[:, f'b_min'] = 0
    i_df.loc[:, f'b_max'] = 20
    

    return i_df



def main():
    root = tk.Tk()

    # Setting some window properties
    root.title("Voigt fitter")
    # root.geometry("1400x1000")
    root.geometry("700x500")

    # Create two labels
    lbl = tk.Label(root, text="Select a species.")
    lbl.grid()

    # function to display text when
    # button is clicked
    def clicked(elem):
        lbl.configure(text = f"{elem} selected")
        elem_inds = atomic_line_list[atomic_line_list['Species'] == elem].index
        elem_df = atomic_line_list.loc[elem_inds]
        print(elem_df)

        ext_df = make_default_df(elem_df, 0)

        print(ext_df)

    # button widget with red color text
    # inside
    btn = tk.Button(root, text = "KI", fg = "red", command=lambda: clicked("KI"))
    # set Button grid
    btn.grid(column=1, row=0)
    root.mainloop()


if __name__ == "__main__":
    main()
