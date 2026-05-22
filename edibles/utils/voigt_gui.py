from edibles.utils.voigt_fitting import voigt_fit_wrapper
import tkinter as tk
from importlib.resources import files
import pandas as pd
from edibles.utils.edibles_oracle import EdiblesOracle
from edibles.utils.edibles_spectrum import EdiblesSpectrum
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
from lmfit.model import save_modelresult, load_modelresult


atomic_line_file = files('edibles') / 'data/auxiliary_data/line_catalogs/edibles_linelist_atoms.csv'
atomic_line_list = pd.read_csv(atomic_line_file)
atomic_line_list = atomic_line_list.dropna(subset=['Gamma'])

# remove lines which are contaminated by telluric lines
atomic_line_list = atomic_line_list.loc[~atomic_line_list['WavelengthAir'].between(7664, 7666)]
atomic_line_list = atomic_line_list.loc[~atomic_line_list['WavelengthAir'].between(4044, 4045)]

molecular_line_file = files('edibles') / 'data/auxiliary_data/line_catalogs/edibles_linelist_molecules.csv'
molecular_line_list = pd.read_csv(molecular_line_file)

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
    i_df.loc[:, f'v_rad_init'] = 0.0
    # i_df.loc[:, f'v_rad_min'] = -100.0
    # i_df.loc[:, f'v_rad_max'] = 100.0
    i_df.loc[:, 'b_comp'] = v_comp
    i_df.loc[:, f'b_init'] = 0.1
    i_df.loc[:, f'b_min'] = 0.0
    i_df.loc[:, f'b_max'] = 6

    # make wavelength range +- 150 km/s around line center
    c = 299792.458 # speed of light in km/s
    i_df.loc[:, 'w_min'] = i_df.loc[:, 'WavelengthAir'] * (1 - 50/c)
    i_df.loc[:, 'w_max'] = i_df.loc[:, 'WavelengthAir'] * (1 + 50/c)

    # if wavelength ranges overlap, merge them
    i_df = i_df.sort_values(by='w_min').reset_index(drop=True)
    for i, row in i_df.iterrows():
        w_min = row['w_min']
        w_max = row['w_max']
        if i + 1 < len(i_df):
            if i_df.loc[i+1, 'w_min'] < w_max:
                i_df.loc[i+1, 'w_min'] = w_min
                i_df.loc[i, 'w_max'] = i_df.loc[i+1, 'w_max']
    
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


def coadd_spectra(spectra: list, ref_spec_num=0, return_error=True) -> np.array:
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
    masked_error_list = 1/masked_weight_list

    coadd_flux = np.ma.sum(masked_flux_list, axis=0)

    coadd_error = np.sqrt(np.ma.sum(masked_error_list**2, axis=0))

    if return_error:
        return np.array([x, coadd_flux, coadd_error])
    else:
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
    root.v_rad_active = False
    root.w_range_active = False
    root.result = None
    root.errorbar = True

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
        print(ext_df)

        # if w_min, w_max is aready changed in fit_df, copy the values to ext_df
        if not root.fit_df.empty:
            for i, i_row in ext_df.iterrows():
                for _, k_row in root.fit_df.iterrows():
                    if i_row['Species'] == k_row['Species'] and i_row['WavelengthAir'] == k_row['WavelengthAir']:
                        ext_df.loc[i, 'w_min'] = k_row['w_min']
                        ext_df.loc[i, 'w_max'] = k_row['w_max']
            

        # load spectra within wavelength range
        root.fit_df = pd.concat([root.fit_df, ext_df], ignore_index=True)
        print_msg(f'Adding a species {elem} to fit_df.')
        print(root.fit_df)
        # plot_fit_info()

    def add_ch_plus():
        molec = 'CH+'
        elem_lbl.configure(text = f"{molec} selected")
        print(molecular_line_list)
        elem_inds = molecular_line_list[molecular_line_list['Species'] == molec].index
        elem_df = molecular_line_list.loc[elem_inds]
        print(elem_df)
        elem_df = elem_df.loc[((elem_df.loc[:, 'WavelengthAir'] > 4229) & (elem_df.loc[:, 'WavelengthAir'] < 4233)) |
                              ((elem_df.loc[:, 'WavelengthAir'] > 3957) & (elem_df.loc[:, 'WavelengthAir'] < 3958))]
        print(elem_df)

        elem_df.replace('CH+', 'CHplus', inplace=True)

        v_comp = root.fit_df['v_comp'].max() + 1 if len(root.fit_df) > 0 else 0

        # Make initial dataframe (include wavelength range)
        ext_df = make_default_df(elem_df, v_comp)
        print("ext_df", ext_df)

        # if w_min, w_max is aready changed in fit_df, copy the values to ext_df
        if not root.fit_df.empty:
            for i, i_row in ext_df.iterrows():
                for _, k_row in root.fit_df.iterrows():
                    if i_row['Species'] == k_row['Species'] and i_row['WavelengthAir'] == k_row['WavelengthAir']:
                        ext_df.loc[i, 'w_min'] = k_row['w_min']
                        ext_df.loc[i, 'w_max'] = k_row['w_max']
            

        # load spectra within wavelength range
        root.fit_df = pd.concat([root.fit_df, ext_df], ignore_index=True)
        print_msg(f'Adding a species {molec} to fit_df.')
        print(root.fit_df)


    # Button for elements
    ki_btn = tk.Button(root, text = "KI", fg = "red", command=lambda: add_elem("KI"))
    ki_btn.grid(column=0, row=1)

    nai_btn = tk.Button(root, text = "NaI", fg = "green", command=lambda: add_elem("NaI"))
    nai_btn.grid(column=0, row=2)

    caI_btn = tk.Button(root, text = "CaI", fg = "blue", command=lambda: add_elem("CaI"))
    caI_btn.grid(column=0, row=3)

    fei_btn = tk.Button(root, text = "FeI", fg = "blue", command=lambda: add_elem("FeI"))
    fei_btn.grid(column=0, row=4)

    tiii_btn = tk.Button(root, text = "TiII", fg = "purple", command=lambda: add_elem("TiII"))
    tiii_btn.grid(column=0, row=5)

    ch_plus_btn = tk.Button(root, text = "CH$^+$", fg = "purple", command=lambda: add_ch_plus())
    ch_plus_btn.grid(column=0, row=6)

    elnum = 6

    # Text box for star name
    star_lbl = tk.Label(root, text="Enter star name:")
    star_lbl.grid(column=0, row=elnum+1)

    star_entry = tk.Entry(root, width=10)
    star_entry.grid(column=0, row=elnum+2)

    # the figure that will contain the plot ==============================================================
    nrows = 2
    ncols = 3
    root.fig, root.axs = plt.subplots(nrows=nrows, ncols=ncols)
    # creating the Tkinter canvas
    root.canvas = FigureCanvasTkAgg(root.fig, master = root)  
    # containing the Matplotlib figure
    root.canvas.draw()
    # placing the canvas on the Tkinter window
    root.canvas.get_tk_widget().grid(column=1, row=1, rowspan=20, sticky='nesw')
    # creating the Matplotlib toolbar
    toolbar_frame = tk.Frame(master=root)
    toolbar_frame.grid(column=1, row=1)
    toolbar = NavigationToolbar2Tk(root.canvas, toolbar_frame)
    toolbar.update()

    def plot_fit_info():
        """
        Plotting relevant fitting info in the spectrum plots.
        """
        print('Plotting fit info')
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
        for plot1 in root.axs.flatten():
            plot1.clear()
        for i, wave_range in range_df.iterrows():
            pythia = EdiblesOracle()
            file_list = pythia.getFilteredObsList(object=[root.star_name], OrdersOnly=True, Wave=np.mean(wave_range), closest_order=True)
            file_lists.append(file_list)

            # getting the subplot
            j = i // ncols
            if nrows == 1:
                plot1 = root.axs[i]
            else:
                plot1 = root.axs[j, i % ncols]
            plot1.clear()
            root.vlines = {}
            spectra = []
            if len(file_list) == 0:
                print_msg(f"No spectra found for {root.star_name} in wavelength range {wave_range['w_min']:.2f} - {wave_range['w_max']:.2f}.")

            for file in file_list:
                print(f"Loaded file: {file}, wave range: {wave_range['w_min']:.2f} - {wave_range['w_max']:.2f}")
                spec = EdiblesSpectrum(file)
                file = Path(file)

                if spec.c_flux is None:
                    spec = np.array([spec.bary_wave, spec.flux, spec.flux_err])
                else:
                    spec = np.array([spec.bary_wave, spec.c_flux, spec.flux_err])

                # crop spectrum
                spec = util_functions.crop_spectrum(spec, *wave_range)
                # my_order = np.nanmedian(spec[4])
                # spec = spec[:, spec[4]==my_order]
                # if len(spec) > 5:
                #     if not np.isnan(spec[6]).all():
                #         spec[1] = spec[6]

                if 3300 < np.mean(wave_range) < 3305:
                    spec[0] = transformations.doppler_shift_wl(spec[0], -1)
                    # spec[2] /= 10
                # plot the spectrum in the GUI using matplotlib
                # plotting the graph
                if root.errorbar:
                    plot1.errorbar(spec[0], spec[1], yerr=spec[2], label = file.name, alpha=0.5)
                else:
                    plot1.plot(spec[0], spec[1], label = file.name, alpha=0.5)

                spectra.append(np.array([spec[0], spec[1], spec[2]]))  # wavelength, flux, error

            # coad spectra
            coadded_spec = coadd_spectra(spectra)
            if root.errorbar:
                plot1.errorbar(coadded_spec[0], coadded_spec[1], yerr=coadded_spec[2], label='Coadd', color='k')
            else:
                plot1.plot(coadded_spec[0], coadded_spec[1], label='Coadd', color='k')
            coadded_spectra.append(coadded_spec)
            plot1.legend()

        root.canvas.draw()
        print('file lists:', file_lists)
        root.fit_spec = np.concatenate(coadded_spectra, axis=1)


    load_btn = tk.Button(root, text="Load Spectrum", command=load_spectrum)
    load_btn.grid(column=0, row=elnum+3)

    # set starting values for fit using plotted spectrum
    # for each component v_comp
    v_rad_comp_entry = tk.Entry(root)
    v_rad_comp_entry.grid(column=0, row=elnum+6)

    root.span = []

    def set_v_rad_init():
        for i, _ in enumerate(root.span):
            root.span[i].set_visible(False)


        root.v_rad_active = True
        # select initial wavelength from plot
        def onclick(event):
            if root.v_rad_active:
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
                print('v_comp_sub_df.index', v_comp_sub_df.index)
                print(root.fit_df)
                for i, row in root.fit_df.iterrows():
                    if row['v_comp'] == int(v_comp):
                        root.fit_df.loc[i, f'v_rad_init'] = v_rad_init

                print(root.fit_df)
                plot_fit_info()
                root.v_rad_active = False


        root.cid = root.canvas.mpl_connect('button_press_event', onclick)

        root.canvas.draw()

    v_rad_btn = tk.Button(root, text='v_rad init', command=set_v_rad_init)
    v_rad_btn.grid(column=0, row=elnum+5)

    # change wavelength range
    def range_function():
        if root.cid is not None:
            root.canvas.mpl_disconnect(root.cid)
        root.w_range_active = True
        print_msg("Select wavelength range by clicking and dragging on the plot.")
        # make range selection tool using matplotlib span selector
        # apply it on canvas
        def onselect(xmin, xmax):
            if root.w_range_active:
                print(f'Selected wavelength range: {xmin:.2f} - {xmax:.2f}')
                # update fit_df with new wavelength range for all lines
                for i in range(len(root.fit_df)):
                    # change wavelength range to selected range if the prior range overlaps with the selected range
                    if root.fit_df.loc[i, 'w_min'] < xmax and root.fit_df.loc[i, 'w_max'] > xmin:
                        root.fit_df.loc[i, 'w_min'] = xmin
                        root.fit_df.loc[i, 'w_max'] = xmax

                print_msg(f"Updated wavelength range for {len(root.fit_df)} lines.")
                print(root.fit_df)    
                root.w_range_active = False
                range_df = root.fit_df[['w_min', 'w_max']].drop_duplicates().reset_index(drop=True).sort_values(by=['w_min'])

                new_spec_list = []
                for i, wave_range in range_df.iterrows():
                    spec = util_functions.crop_spectrum(root.fit_spec, *wave_range)
                    new_spec_list.append(spec)
                    # getting the subplot
                    j = i // ncols
                    if nrows == 1:
                        plot1 = root.axs[i]
                    else:
                        plot1 = root.axs[j, i % ncols]
                    plot1.clear()
                    if root.errorbar:
                        plot1.errorbar(spec[0], spec[1], yerr=spec[2], color='k', label = 'Data')
                    else:
                        plot1.plot(spec[0], spec[1], color='k', label = 'Data')

                    if root.result is not None:
                        if len(root.fit_spec[0]) == len(root.result.best_fit):
                            plot1.plot(spec[0], root.result.best_fit[(root.fit_spec[0] >= wave_range['w_min']) & (root.fit_spec[0] <= wave_range['w_max'])], label='Fit', color='r')
                    plot1.legend()
                    root.canvas.draw()



                root.fit_spec = np.concatenate(new_spec_list, axis=1)



        root.span.clear()
        for _, ax in enumerate(root.axs.flatten()):
            selector = SpanSelector(ax, onselect, 'horizontal', useblit=True, props=dict(alpha=0.5, facecolor='red'))
            root.span.append(selector)

        root.canvas.draw_idle()

    range_btn = tk.Button(root, text='Set wavelength range', command=range_function)
    range_btn.grid(column=0, row=elnum+7)


    # fitting
    def fit_spectrum():
        """
        Fitting the spectrum with the voigt model.
        Assigns result to root.result and prints best fit values to console.
        """
        if root.fit_spec is None:
            print_msg("No spectrum loaded.")
            return
        
        print_msg('Fitting voigt model.')
        
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
            root.vlines = {}

            plot_spec = root.fit_spec[:, (root.fit_spec[0] >= wave_range['w_min']) & (root.fit_spec[0] <= wave_range['w_max'])]
            # plot1.plot(plot_spec[0], plot_spec[1], 'k', label = 'Data')
            if root.errorbar:
                plot1.errorbar(plot_spec[0], plot_spec[1], yerr=plot_spec[2], color='k', label = 'Data')
            else:
                plot1.plot(plot_spec[0], plot_spec[1], color='k', label = 'Data')
            plot1.plot(plot_spec[0], result.best_fit[(root.fit_spec[0] >= wave_range['w_min']) & (root.fit_spec[0] <= wave_range['w_max'])], label='Fit', color='r')
            plot1.legend()
            root.canvas.draw()

        print(result.best_values)
        root.result = result

    fit_btn = tk.Button(root, text="Fit Spectrum", command=fit_spectrum)
    fit_btn.grid(column=0, row=elnum+4)

    def save_function():
        # Save fitting results in csv file
        res_df = results_to_df(root.result, root.fit_df)
        elem_list = root.fit_df.loc[:, 'Species'].unique()
        elem_str = '_'.join(elem_list)
        print(elem_str)
        res_df.to_csv(fitting_dir / f'{root.star_name}_{elem_str}.csv', index=False)

        # Save spectrum
        np.savetxt(fitting_dir / f'{root.star_name}_{elem_str}.dat', root.fit_spec.T)

        # Save model result 
        save_modelresult(root.result, fitting_dir / f'{root.star_name}_{elem_str}.sav')


    save_btn = tk.Button(root, text="Save fit results", command=save_function)
    save_btn.grid(column=0, row=elnum+8)

    def load_fit_result():
        """
        Loads saved fit results from csv file and assigns it to root.fit_df. 
        The fit results can then be plotted by clicking the "Fit Spectrum" button after loading a spectrum.
        """
        star_name = star_entry.get()
        elem_list = root.fit_df.loc[:, 'Species'].unique()
        elem_str = '_'.join(elem_list)
        print_msg(f'Loading fit results: {fitting_dir / f"{star_name}_{elem_str}.csv"}')
        if star_name is None:
            print("No star name entered.")
            return
        try:
            # load fit df
            fit_df = pd.read_csv(fitting_dir / f'{star_name}_{elem_str}.csv')
            print(fit_df)
            fit_df['v_rad_init'] = fit_df['v_rad_fit']
            fit_df['b_init'] = fit_df['b_fit']
            fit_df['n_init'] = fit_df['n_fit']

            root.fit_df = fit_df

            # load spectrum
            root.fit_spec = np.genfromtxt(fitting_dir / f'{star_name}_{elem_str}.dat', unpack=True)

            # load model result
            root.result = load_modelresult(fitting_dir / f'{star_name}_{elem_str}.sav')

            range_df = root.fit_df[['w_min', 'w_max']].drop_duplicates().reset_index(drop=True).sort_values(by=['w_min'])

            # update fit_df with new wavelength range for all lines
            for i, wave_range in range_df.iterrows():
                spec = util_functions.crop_spectrum(root.fit_spec, *wave_range)
                # getting the subplot
                j = i // ncols
                if nrows == 1:
                    plot1 = root.axs[i]
                else:
                    plot1 = root.axs[j, i % ncols]
                plot1.clear()
                if root.errorbar:
                    plot1.errorbar(spec[0], spec[1], yerr=spec[2], color='k', label = 'Data')
                else:
                    plot1.plot(spec[0], spec[1], color='k', label = 'Data')

                if root.result is not None:
                    if len(root.fit_spec[0]) == len(root.result.best_fit):
                        plot1.plot(spec[0], root.result.best_fit[(root.fit_spec[0] >= wave_range['w_min']) & (root.fit_spec[0] <= wave_range['w_max'])], label='Fit', color='r')
                plot1.legend()
                root.canvas.draw()

        except FileNotFoundError:
            print(f"No fit results found for {star_name}.")

    load_res_btn = tk.Button(root, text="Load fit results", command=load_fit_result)
    load_res_btn.grid(column=0, row=elnum+9)

    def clear_df_func():
        root.fit_df = pd.DataFrame()
        root.fit_spec = None
        root.result = None
        print_msg('Clearing the present fit_df DataFrame. A new fit can be started.')
    
    clear_df_btn = tk.Button(root, text="Clear fit DataFrame", command=clear_df_func)
    clear_df_btn.grid(column=0, row=elnum+10)

    # set starting values for fit using plotted spectrum
    # for each component v_comp

    # v_rad_window_comp_entry = tk.Entry(root)
    # v_rad_window_comp_entry.grid(column=0, row=elnum+12)

    # def set_v_rad_shift():
    #     """
    #     Setting a radial valocity shift for a specific fitting window.
    #     """
    #     for i, _ in enumerate(root.span):
    #         root.span[i].set_visible(False)


    #     root.v_rad_active = True

    #     v_comp = v_rad_window_comp_entry.get()
    #     v_comp_sub_df = root.fit_df.loc[root.fit_df['v_comp'] == int(v_comp)]

    #     # select initial wavelength from plot
    #     def onclick(event):
    #         if root.v_rad_active:
    #             ix= event.xdata
    #             print(f'Clicked at x = {ix}')

    #             # find correspinding wavelength in df
    #             line_idx = (root.fit_df['WavelengthAir'] - ix).abs().idxmin()
    #             print(f'Selected line: {root.fit_df.loc[line_idx, "WavelengthAir"]}')

    #             # calculate doppler shift between selected wavelength and line center
    #             c = 299792.458 # speed of light in km/s
    #             line_center = root.fit_df.loc[line_idx, 'WavelengthAir']
    #             v_rad_init = (ix - line_center) / line_center * c
    #             print(f'Calculated radial velocity: {v_rad_init:.2f} km/s')

    #             # update fit_df with new v_rad_init for selected component
    #             print('v_comp_sub_df.index', v_comp_sub_df.index)
    #             print(root.fit_df)
    #             for i, row in root.fit_df.iterrows():
    #                 if row['v_comp'] == int(v_comp):
    #                     root.fit_df.loc[i, f'v_rad_init'] = v_rad_init

    #             print(root.fit_df)
    #             plot_fit_info()
    #             root.v_rad_active = False


    #     root.cid = root.canvas.mpl_connect('button_press_event', onclick)

    #     root.canvas.draw()

    # window_shift_btn = tk.Button(root, text="Apply shift to spectrum", command=set_v_rad_shift)
    # window_shift_btn.grid(column=0, row=elnum+11)

    # change wavelength range
    def cont_range_function():
        if root.cid is not None:
            root.canvas.mpl_disconnect(root.cid)
        root.w_range_active = True
        root.range_counter = 0
        range_list = []
        print_msg("Select wavelength range by clicking and dragging on the plot.")
        # make range selection tool using matplotlib span selector
        # apply it on canvas
        def onselect(xmin, xmax):
            range_df = root.fit_df[['w_min', 'w_max']].drop_duplicates().reset_index(drop=True).sort_values(by=['w_min'])
            if root.range_counter == 0:
                range_list.append([xmin, xmax])
                root.range_counter += 1
            elif root.range_counter == 1:
                range_list.append([xmin, xmax])
                print(f'Selected wavelength range: {xmin:.2f} - {xmax:.2f}')
                # update fit_df with new wavelength range for all lines
                for i, wave_range in range_df.iterrows():
                    # adding the subplot
                    j = i // ncols
                    if nrows == 1:
                        plot1 = root.axs[i]
                    else:
                        plot1 = root.axs[j, i % ncols]
                    if (wave_range['w_min'] < xmin) & (wave_range['w_max'] > xmax):
                        cut_spec = root.fit_spec[:, (root.fit_spec[0] >= wave_range['w_min']) & (root.fit_spec[0] <= wave_range['w_max'])]
                        cont_anchors = []
                        for range in range_list:
                            co = util_functions.crop_spectrum(cut_spec, range[0], range[1])
                            point = np.nanmean(co, axis=1)
                            cont_anchors.append(point)
                        print(cont_anchors)
                        cut_spec_norm = util_functions.normalize_spectrum_linear(cut_spec, cont_anchors[0], cont_anchors[1], additional_normalized_columns=[2])
                        root.fit_spec[:, (root.fit_spec[0] >= wave_range['w_min']) & (root.fit_spec[0] <= wave_range['w_max'])] = cut_spec_norm
                        plot1.clear()
                        root.vlines = {}

                        plot_spec = root.fit_spec[:, (root.fit_spec[0] >= wave_range['w_min']) & (root.fit_spec[0] <= wave_range['w_max'])]
                        if root.errorbar:
                            plot1.errorbar(plot_spec[0], plot_spec[1], yerr=plot_spec[2], color='k', label = 'Data')
                        else:
                            plot1.plot(plot_spec[0], plot_spec[1], color='k', label = 'Data')

                        if root.result is not None:
                            if len(root.fit_spec[0]) == len(root.result.best_fit):
                                plot1.plot(plot_spec[0], root.result.best_fit[(root.fit_spec[0] >= wave_range['w_min']) & (root.fit_spec[0] <= wave_range['w_max'])], label='Fit', color='r')
                        plot1.legend()
                        root.canvas.draw()


                print_msg(f"Updated wavelength range for {len(root.fit_df)} lines.")
                print(root.fit_df)    
                root.range_counter += 1

        root.span.clear()
        for _, ax in enumerate(root.axs.flatten()):
            selector = SpanSelector(ax, onselect, 'horizontal', useblit=True, props=dict(alpha=0.5, facecolor='red'))
            root.span.append(selector)

        root.canvas.draw_idle()

    range_btn = tk.Button(root, text='Select continuum ranges', command=cont_range_function)
    range_btn.grid(column=0, row=elnum+12)

    def mask_function():
        if root.cid is not None:
            root.canvas.mpl_disconnect(root.cid)
        root.w_range_active = True
        print_msg("Select wavelength range by clicking and dragging on the plot.")
        # make range selection tool using matplotlib span selector
        # apply it on canvas
        def onselect(xmin, xmax):
            print(f'Selected wavelength range: {xmin:.2f} - {xmax:.2f}')
            spec_inds = root.fit_spec[0].searchsorted([xmin, xmax])
            print(spec_inds)
            root.fit_spec[2, spec_inds[0]:spec_inds[1]] = np.inf
            range_df = root.fit_df[['w_min', 'w_max']].drop_duplicates().reset_index(drop=True).sort_values(by=['w_min'])

            for i, wave_range in range_df.iterrows():
                # adding the subplot
                j = i // ncols
                if nrows == 1:
                    plot1 = root.axs[i]
                else:
                    plot1 = root.axs[j, i % ncols]
                if (wave_range['w_min'] < xmin) & (wave_range['w_max'] > xmax):
                    plot1.clear()
                    root.vlines = {}

                    plot_spec = root.fit_spec[:, (root.fit_spec[0] >= wave_range['w_min']) & (root.fit_spec[0] <= wave_range['w_max'])]
                    if root.errorbar:
                        plot1.errorbar(plot_spec[0], plot_spec[1], yerr=plot_spec[2], color='k', label = 'Data')
                    else:
                        plot1.plot(plot_spec[0], plot_spec[1], color='k', label = 'Data')

                    if root.result is not None:
                        if len(root.fit_spec[0]) == len(root.result.best_fit):
                            plot1.plot(plot_spec[0], root.result.best_fit[(root.fit_spec[0] >= wave_range['w_min']) & (root.fit_spec[0] <= wave_range['w_max'])], label='Fit', color='r')
                    plot1.legend()
                    root.canvas.draw()

                    print_msg(f"Added mask to fitting spectrum.")
                root.w_range_active = False

        root.span.clear()
        for _, ax in enumerate(root.axs.flatten()):
            selector = SpanSelector(ax, onselect, 'horizontal', useblit=True, props=dict(alpha=0.5, facecolor='red'))
            root.span.append(selector)

        root.canvas.draw_idle()

    mask_btn = tk.Button(root, text='Select mask ranges', command=mask_function)
    mask_btn.grid(column=0, row=elnum+13)

    root.grid_columnconfigure(1, weight=1)
    root.grid_rowconfigure(elnum+14, weight=1)

    root.mainloop()


if __name__ == "__main__":
    main()
