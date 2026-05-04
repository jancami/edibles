from edibles.utils.voigt_fitting import voigt_fit_wrapper
import tkinter as tk
from importlib.resources import files
import pandas as pd
from edibles.utils.edibles_oracle import EdiblesOracle
import numpy as np
from edibles.projects.edr5_integration import dr5_io, util_functions
from edibles import DATADIR
from PyAstronomy import pyasl
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)


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

    # make wavelength range +- 150 km/s around line center
    c = 299792.458 # speed of light in km/s
    i_df.loc[:, 'w_min'] = i_df['WavelengthAir'] * (1 - 150/c)
    i_df.loc[:, 'w_max'] = i_df['WavelengthAir'] * (1 + 150/c)

    # if wavelength ranges overlap, merge them
    i_df = i_df.sort_values(by='w_min').reset_index(drop=True)
    merged_ranges = []
    current_range = [i_df.loc[0, 'w_min'], i_df.loc[0, 'w_max']]
    for i in range(1, len(i_df)):
        w_min = i_df.loc[i, 'w_min']
        w_max = i_df.loc[i, 'w_max']
        if w_min <= current_range[1]: # if ranges overlap
            current_range[1] = max(current_range[1], w_max) # merge ranges
        else:
            merged_ranges.append(current_range)
            current_range = [w_min, w_max]
    merged_ranges.append(current_range) # add last range        
    

    return i_df



def main():
    root = tk.Tk()
    # root.state('iconic')

    global fit_df
    fit_df = pd.DataFrame()


    # Setting some window properties
    root.title("Voigt fitter")
    root.geometry("1600x1000")
    # root.geometry("700x500")

    elem_lbl = tk.Label(root, text="Select a species.")
    elem_lbl.grid(column=0, row=0)
    # elem_lbl.pack()

    # function to display text when
    # button is clicked
    def clicked(elem):
        global fit_df
        elem_lbl.configure(text = f"{elem} selected")
        elem_inds = atomic_line_list[atomic_line_list['Species'] == elem].index
        elem_df = atomic_line_list.loc[elem_inds]
        print(elem_df)

        # Make initial dataframe (include wavelength range)
        ext_df = make_default_df(elem_df, 0)

        # load spectra within wavelength range
        fit_df = pd.concat([fit_df, ext_df], ignore_index=True)
        print(fit_df)

    # Button for elements
    elem_btn = tk.Button(root, text = "KI", fg = "red", command=lambda: clicked("KI"))
    elem_btn.grid(column=0, row=1)
    # btn.pack()

    # Text box for star name
    star_lbl = tk.Label(root, text="Enter star name:")
    star_lbl.grid(column=0, row=2)
    # star_lbl.pack()

    star_entry = tk.Entry(root, width=10)
    star_entry.grid(column=0, row=3)
    # star_entry.pack()

        # the figure that will contain the plot
    fig = Figure(figsize = (15, 10), dpi = 100)

    # creating the Tkinter canvas
    canvas = FigureCanvasTkAgg(fig, master = root)  
    # containing the Matplotlib figure
    canvas.draw()

    # placing the canvas on the Tkinter window
    canvas.get_tk_widget().grid(column=1, row=1, rowspan=4)
    # canvas.get_tk_widget().pack()

    # creating the Matplotlib toolbar
    toolbar_frame = tk.Frame(master=root)
    toolbar_frame.grid(column=1, row=0)
    # toolbar = NavigationToolbar2TkAgg(canvas, toolbar_frame)
    toolbar = NavigationToolbar2Tk(canvas, toolbar_frame)
    toolbar.update()

    # placing the toolbar on the Tkinter window
    # canvas.get_tk_widget().grid(column=2, row=2)
    # canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)




    # make button to get star name and load spectrum
    def load_spectrum():
        star_name = star_entry.get()
        print(f"Loading spectrum for {star_name}...")
        # load spectrum for star_name
        # extract fitting ranges
        range_df = fit_df[['w_min', 'w_max']].drop_duplicates().reset_index(drop=True).sort_values(by=['w_min'])

        file_lists = []
        coadded_spectra = []
        for i, wave_range in range_df.iterrows():
            pythia = EdiblesOracle()
            file_list = pythia.getFilteredObsList(object=[star_name], MergedOnly=True, Wave=np.mean(wave_range))
            file_lists.append(file_list)

            # adding the subplot
            plot1 = fig.add_subplot(2, 2, i+1)
            spectra = []

            for file in file_list:
                print(f"Loaded file: {file}, wave range: {wave_range['w_min']:.2f} - {wave_range['w_max']:.2f}")
                spec = dr5_io.read_combined_spec(DATADIR / file, bary_corr=True)
                spec = util_functions.crop_spectrum(spec, *wave_range)
                my_order = np.nanmedian(spec[4])
                spec = spec[:, spec[4]==my_order]
                if len(spec) > 5:
                    if not np.isnan(spec[6]).all():
                        spec[1] = spec[6]

                x, y = pyasl.equidistantInterpolation(spec[0], spec[1], '2x')
                # plot the spectrum in the GUI using matplotlib

                # plotting the graph
                plot1.plot(x, y)

                spectra.append(np.array([x, y]))

            coadded_spec = pyasl.coaddSpec(spectra, method='mean')
            coadded_spectra.append(coadded_spec)

        canvas.draw()
        print('file lists:', file_lists)


    load_btn = tk.Button(root, text="Load Spectrum", command=load_spectrum)
    load_btn.grid(column=0, row=4)
    # load_btn.pack()


    root.mainloop()


if __name__ == "__main__":
    main()
