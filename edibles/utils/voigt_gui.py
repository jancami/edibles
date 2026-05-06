from edibles.utils.voigt_fitting import voigt_fit_wrapper
import tkinter as tk
from importlib.resources import files
import pandas as pd
from edibles.utils.edibles_oracle import EdiblesOracle
import numpy as np
from edibles.projects.edr5_integration import dr5_io, util_functions
from edibles import DATADIR
from PyAstronomy import pyasl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from pathlib import Path
from scipy.interpolate import interp1d
from matplotlib.widgets import SpanSelector
from edibles.utils import transformations


atomic_line_file = files('edibles') / 'data/auxiliary_data/line_catalogs/edibles_linelist_atoms.csv'
atomic_line_list = pd.read_csv(atomic_line_file)
atomic_line_list = atomic_line_list.dropna(subset=['Gamma'])

# remove lines which are contaminated by telluric lines
atomic_line_list = atomic_line_list.loc[~atomic_line_list['WavelengthAir'].between(7664, 7666)]

fitting_dir = files('edibles') / 'data/voigt_fitting_data'

def make_default_df(in_df: pd.DataFrame, v_comp: int) -> pd.DataFrame:
    """
    Create a default DataFrame for Voigt fitting.

    Parameters
    ----------
    in_df : pd.DataFrame
        DataFrame loaded from a line list.
    v_comp : int
        velocity component number (used to link lines that belong to the same component and set the same initial v_rad and b values)

    Returns
    -------
    pd.DataFrame
        DataFrame with default values for fitting parameters and wavelength ranges for fitting.
    """

    i_df = in_df.copy()
    i_df.loc[:, 'v_comp'] = v_comp
    i_df.loc[:, f'v_rad_init'] = 0
    i_df.loc[:, f'v_rad_min'] = -100
    i_df.loc[:, f'v_rad_max'] = 100
    i_df.loc[:, 'b_comp'] = v_comp
    i_df.loc[:, f'b_init'] = 0.001
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

def resample(spectrum: np.array, wave_new: np.array, assume_sorted=True) -> np.array:
    """
    Resample a spectrum to given wavelength points.

    Parameters
    ----------
    spectrum : np.array
        Input spectrum (*np.array([wave, flux])*)

    wave_new : numpy array of wavelength points for the new spectrum

    assume_sorted : bool

    Returns
    -------
    np_array:
        resampled spectrum
    """
    new_cols = []
    for col in spectrum[1:]:
        f = interp1d(spectrum[0], col, assume_sorted=assume_sorted, bounds_error=False, fill_value="extrapolate")
        col_new = f(wave_new)
        new_cols.append(col_new)
    return np.array([wave_new, *new_cols])


def coadd_spectra(spectra: list, ref_spec_num=0) -> np.array:
    """
    Coadds a list of spectra.
    Before coadding, the spectra get resampled to the same wavelength points using the spectrum with index ref_spec_num as reference.

    Parameters
    ----------
    spectra : list
        List of spectra to coadd.
    ref_spec_num : int, optional
        Index of the reference spectrum, by default 0

    Returns
    -------
    np.array
        Coadded spectrum
    """
    ref_spec = spectra[ref_spec_num]
    x, _ = pyasl.equidistantInterpolation(ref_spec[0], ref_spec[1], '2x')

    flux_list = []
    weight_list = []
    for spec in spectra:
        res_spec = resample(spec, x)
        flux_list.append(res_spec[1])
        weight_list.append(1/res_spec[2])

    
    flux_list = np.array(flux_list)
    weight_list = np.array(weight_list)

    masked_flux_list = np.ma.masked_array(flux_list, mask=np.isnan(flux_list))
    masked_weight_list = np.ma.masked_array(weight_list, mask=np.isnan(weight_list))

    coadd_flux = np.ma.average(masked_flux_list, weights=masked_weight_list, axis=0)

    return np.array([x, coadd_flux])


def results_to_df(result, fit_df):
    """
    Merges the fitting DataFrame and the results to one DataFrame containing all the information.

    Parameters
    ----------
    result : lmfit.result
        lmfit result object containing the best fit values for all parameters and the best fit model.
    fit_df : pd.DataFrame
        DataFrame with default values for fitting parameters and wavelength ranges for fitting.

    Returns
    -------
    pd.DataFrame
        Merged DataFrame containing all fitting information.
    """
    best_values = result.best_values

    print('results to file')

    range_df = fit_df[['w_min', 'w_max']].drop_duplicates().reset_index(drop=True).sort_values(by=['w_min'])

    c_comps = fit_df[['v_comp', 'Species']].drop_duplicates().reset_index(drop=True)


    # continuum parameters
    for i, wave_range in range_df.iterrows():
        for k, row in fit_df.iterrows():
            if row['w_min'] == wave_range['w_min'] and row['w_max'] == wave_range['w_max']:
                fit_df.loc[k, f'cont'] = best_values[f'cont_{i}']
                fit_df.loc[k, f'slope'] = best_values[f'slope_{i}']

    # v_rad, b and N for each component
    for i, v_comp in c_comps.iterrows():
        for k, row in fit_df.iterrows():
            if row['v_comp'] == v_comp['v_comp'] and row['Species'] == v_comp['Species']:
                fit_df.loc[k, f'v_rad_fit'] = best_values[f'v_rad_{v_comp["v_comp"]}']
                fit_df.loc[k, f'b_fit'] = best_values[f'b_{v_comp["v_comp"]}']
                fit_df.loc[k, f'n_fit'] = best_values[f'n_{v_comp["v_comp"]}_{v_comp["Species"]}']

    print(fit_df)


    return fit_df

def main():
    root = tk.Tk()

    # defining global variables for the GUI
    root.fit_df = pd.DataFrame()
    root.fit_spec = None
    root.cid = None
    root.vlines = {}

    # Setting some window properties
    root.title("Voigt fitter")
    root.geometry("1600x1000")

    # make label for important messages
    msg_lbl = tk.Label(root, text="")
    msg_lbl.grid(column=1, row=0)

    def print_msg(msg: str):
        msg_lbl.config(text=msg)
        print(msg)

    # Button for selecting element to fit. Only KI for now, but can be easily extended to include more elements.
    elem_lbl = tk.Label(root, text="Select a species.")
    elem_lbl.grid(column=0, row=0)

    # function to display text when
    # button is clicked
    def add_elem(elem):
        elem_lbl.configure(text = f"{elem} selected")
        elem_inds = atomic_line_list[atomic_line_list['Species'] == elem].index
        elem_df = atomic_line_list.loc[elem_inds]

        v_comp = root.fit_df['v_comp'].max() + 1 if len(root.fit_df) > 0 else 0

        # Make initial dataframe (include wavelength range)
        ext_df = make_default_df(elem_df, v_comp)

        # load spectra within wavelength range
        root.fit_df = pd.concat([root.fit_df, ext_df], ignore_index=True)
        # plot_fit_info()

    # Button for elements
    elem_btn = tk.Button(root, text = "KI", fg = "red", command=lambda: add_elem("KI"))
    elem_btn.grid(column=0, row=1)

    # Text box for star name
    star_lbl = tk.Label(root, text="Enter star name:")
    star_lbl.grid(column=0, row=2)

    star_entry = tk.Entry(root, width=10)
    star_entry.grid(column=0, row=3)

    # the figure that will contain the plot ==============================================================
    nrows = 1
    ncols = 2
    root.fig, root.axs = plt.subplots(figsize = (15, 8), dpi = 100, nrows=nrows, ncols=ncols)
    # creating the Tkinter canvas
    root.canvas = FigureCanvasTkAgg(root.fig, master = root)  
    # containing the Matplotlib figure
    root.canvas.draw()
    # placing the canvas on the Tkinter window
    root.canvas.get_tk_widget().grid(column=1, row=1, rowspan=20)
    # creating the Matplotlib toolbar
    toolbar_frame = tk.Frame(master=root)
    toolbar_frame.grid(column=1, row=1)
    toolbar = NavigationToolbar2Tk(root.canvas, toolbar_frame)
    toolbar.update()

    def plot_fit_info():
        """
        Plotting relevant fitting info in the spectrum plots.
        """
        print('Plotting fit infooooooooooooooooooooooooooooooooooooooo')
        # c_comps = root.fit_df['v_comp'].drop_duplicates().reset_index(drop=True)

        # iterate through plot windows
        range_df = root.fit_df[['w_min', 'w_max']].drop_duplicates().reset_index(drop=True).sort_values(by=['w_min'])


        for k, row in root.fit_df.iterrows():
            for i, wave_range in range_df.iterrows():
                # getting the subplot axis
                j = i // ncols
                if nrows == 1:
                    plot1 = root.axs[i]
                else:
                    plot1 = root.axs[j, i % ncols]

                print('wave range: ',i)

                if wave_range['w_min'] < row['WavelengthAir'] < wave_range['w_max']:
                    x = transformations.doppler_shift_wl(row['WavelengthAir'], row['v_rad_init'])
                    print(x)
                    key = (i, int(row['v_comp']))
                    print(root.vlines)
                    print(root.vlines.get(key))
                    if root.vlines.get(key) is None:
                        root.vlines[key] = plot1.axvline(x, color='red', linestyle='--')
                        print('xdataaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa')
                        print(root.vlines[key].get_xdata())
                    else:
                        line = root.vlines[key]
                        line.set_xdata([x, x])
        
                    root.canvas.draw()



    # make button to get star name and load spectrum
    # TODO: add option to exclude a spectrum from coaddidtion if it is bad (e.g. by plotting the spectra and letting the user click on the bad spectra)
    def load_spectrum():
        root.star_name = star_entry.get()
        msg_lbl.config(text=f"Loading spectrum for {root.star_name}...")
        print(f"Loading spectrum for {root.star_name}...")
        # load spectrum for star_name
        # extract fitting ranges
        if root.fit_df.empty:
            print_msg("No element selected. Spectra cannot be loaded as no fitting ranges are defined.")
            return
        
        range_df = root.fit_df[['w_min', 'w_max']].drop_duplicates().reset_index(drop=True).sort_values(by=['w_min'])

        file_lists = []
        coadded_spectra = []
        for i, wave_range in range_df.iterrows():
            pythia = EdiblesOracle()
            file_list = pythia.getFilteredObsList(object=[root.star_name], MergedOnly=True, Wave=np.mean(wave_range))
            file_lists.append(file_list)

            # getting the subplot
            j = i // ncols
            if nrows == 1:
                plot1 = root.axs[i]
            else:
                plot1 = root.axs[j, i % ncols]
            plot1.clear()
            spectra = []
            if len(file_list) == 0:
                print_msg(f"No spectra found for {root.star_name} in wavelength range {wave_range['w_min']:.2f} - {wave_range['w_max']:.2f}.")

            for file in file_list:
                file = Path(file)
                print(f"Loaded file: {file}, wave range: {wave_range['w_min']:.2f} - {wave_range['w_max']:.2f}")
                spec = dr5_io.read_combined_spec(DATADIR / file, bary_corr=True)
                spec = util_functions.crop_spectrum(spec, *wave_range)
                my_order = np.nanmedian(spec[4])
                spec = spec[:, spec[4]==my_order]
                if len(spec) > 5:
                    if not np.isnan(spec[6]).all():
                        spec[1] = spec[6]

                # plot the spectrum in the GUI using matplotlib
                # plotting the graph
                plot1.plot(spec[0], spec[1], label = file.name, alpha=0.5)

                spectra.append(np.array([spec[0], spec[1], spec[2]]))  # wavelength, flux, error

            # coad spectra
            coadded_spec = coadd_spectra(spectra)
            plot1.plot(coadded_spec[0], coadded_spec[1], label='Coadd', color='k')
            coadded_spectra.append(coadded_spec)
            plt.legend()

        root.canvas.draw()
        print('file lists:', file_lists)
        root.fit_spec = np.concatenate(coadded_spectra, axis=1)


    load_btn = tk.Button(root, text="Load Spectrum", command=load_spectrum)
    load_btn.grid(column=0, row=4)

    # set starting values for fit using plotted spectrum
    # for each component v_comp
    v_rad_comp_entry = tk.Entry(root)
    v_rad_comp_entry.grid(column=0, row=7)


    def set_v_rad_init():
        # select initial wavelength from plot
        def onclick(event):
            v_comp = v_rad_comp_entry.get()
            v_comp_sub_df = root.fit_df.loc[root.fit_df['v_comp'] == int(v_comp)]

            ix= event.xdata
            print(f'Clicked at x = {ix}')

            # find correspinding wavelength in df
            line_idx = (root.fit_df['WavelengthAir'] - ix).abs().idxmin()
            print(f'Selected line: {root.fit_df.loc[line_idx, "WavelengthAir"]}')

            # calculate doppler shift between selected wavelength and line center
            c = 299792.458 # speed of light in km/s
            line_center = root.fit_df.loc[line_idx, 'WavelengthAir']
            v_rad_init = (ix - line_center) / line_center * c
            print(f'Calculated radial velocity: {v_rad_init:.2f} km/s')

            # update fit_df with new v_rad_init for selected component
            print(v_comp_sub_df.index)
            for i, row in root.fit_df.iterrows():
                if row['v_comp'] == int(v_comp):
                    root.fit_df.loc[i, f'v_rad_init'] = v_rad_init

            print(root.fit_df)
            plot_fit_info()


        root.cid = root.canvas.mpl_connect('button_press_event', onclick)

        root.canvas.draw()

    v_rad_btn = tk.Button(root, text='v_rad init', command=set_v_rad_init)
    v_rad_btn.grid(column=0, row=6)

    root.span = []
    # change wavelength range
    def range_function():
        if root.cid is not None:
            root.canvas.mpl_disconnect(root.cid)
        print_msg("Select wavelength range by clicking and dragging on the plot.")
        # make range selection tool using matplotlib span selector
        # apply it on canvas
        def onselect(xmin, xmax):
            print(f'Selected wavelength range: {xmin:.2f} - {xmax:.2f}')
            # update fit_df with new wavelength range for all lines
            for i in range(len(root.fit_df)):
                # change wavelength range to selected range if the prior range overlaps with the selected range
                if root.fit_df.loc[i, 'w_min'] < xmax and root.fit_df.loc[i, 'w_max'] > xmin:
                    root.fit_df.loc[i, 'w_min'] = xmin
                    root.fit_df.loc[i, 'w_max'] = xmax
            print_msg(f"Updated wavelength range for {len(root.fit_df)} lines.")
            print(root.fit_df)    

        root.span.clear()
        for _, ax in enumerate(root.axs.flatten()):
            selector = SpanSelector(ax, onselect, 'horizontal', useblit=True, props=dict(alpha=0.5, facecolor='red'))
            root.span.append(selector)

        root.canvas.draw_idle()

    range_btn = tk.Button(root, text='Set wavelength range', command=range_function)
    range_btn.grid(column=0, row=8)


    # fitting
    def fit_spectrum():
        """
        Fitting the spectrum with the voigt model.
        Assigns result to root.result and prints best fit values to console.
        """
        if root.fit_spec is None:
            print_msg("No spectrum loaded.")
            return
        
        result = voigt_fit_wrapper(root.fit_df, root.fit_spec)

        range_df = root.fit_df[['w_min', 'w_max']].drop_duplicates().reset_index(drop=True).sort_values(by=['w_min'])

        for i, wave_range in range_df.iterrows():
            # adding the subplot
            j = i // ncols
            if nrows == 1:
                plot1 = root.axs[i]
            else:
                plot1 = root.axs[j, i % ncols]

            plot1.clear()

            plot_spec = root.fit_spec[:, (root.fit_spec[0] >= wave_range['w_min']) & (root.fit_spec[0] <= wave_range['w_max'])]
            plot1.plot(plot_spec[0], plot_spec[1], label = 'Data', alpha=0.5)
            plot1.plot(plot_spec[0], result.best_fit[(root.fit_spec[0] >= wave_range['w_min']) & (root.fit_spec[0] <= wave_range['w_max'])], label='Fit', color='k')
            plot1.legend()
            root.canvas.draw()

        print(result.best_values)
        root.result = result

    fit_btn = tk.Button(root, text="Fit Spectrum", command=fit_spectrum)
    fit_btn.grid(column=0, row=5)

    save_btn = tk.Button(root, text="Save fit results", command=lambda: results_to_df(root.result, root.fit_df).to_csv(fitting_dir / f'{root.star_name}_KI.csv', index=False))
    save_btn.grid(column=0, row=9)

    def load_fit_result():
        """
        Loads saved fit results from csv file and assigns it to root.fit_df. 
        The fit results can then be plotted by clicking the "Fit Spectrum" button after loading a spectrum.
        """
        star_name = star_entry.get()
        print_msg(f'Loading fit results: {fitting_dir / f"{star_name}_KI.csv"}')
        if star_name is None:
            print("No star name entered.")
            return
        try:
            fit_df = pd.read_csv(fitting_dir / f'{star_name}_KI.csv')
            print(fit_df)
            fit_df['v_rad_init'] = fit_df['v_rad_fit']
            fit_df['b_init'] = fit_df['b_fit']
            fit_df['n_init'] = fit_df['n_fit']
            root.fit_df = fit_df
            # plot_fit_info()
        except FileNotFoundError:
            print(f"No fit results found for {star_name}.")

    load_res_btn = tk.Button(root, text="Load fit results", command=load_fit_result)
    load_res_btn.grid(column=0, row=10)

    def clear_df_func():
        root.fit_df = pd.DataFrame()
        print_msg('Clearing the present fit_df DataFrame. A new fit can be started.')
    
    clear_df_btn = tk.Button(root, text="Clear fit DataFrame", command=clear_df_func)
    clear_df_btn.grid(column=0, row=11)

    root.mainloop()


if __name__ == "__main__":
    main()
