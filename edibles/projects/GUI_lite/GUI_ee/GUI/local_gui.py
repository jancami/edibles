"""EDIBLES GUI for spectrum co-adding + Voigt fitting.

Keep this file just about the GUI. The actual science lives elsewhere:
    flux_wave_find.wave_flux_data    -> grab one spectrum from a file
    co_adding_flux.coadd_spectra     -> co-add the spectra you picked
    main_run.get_species_data        -> read species.txt
    main_run.astrovoigtfit_run       -> the joint continuum + Voigt fit
    main.master_function             -> per-species absorption (we also call
                                        this directly to draw the individual
                                        and component-wise overlay curves)

So if you want to change how the fit or co-add behaves, go to those files.
This one is for layout, widgets, plotting style, and wiring things together.

EdiblesApp at the bottom is the app window. AnalysisTab is one analysis
session — you can open multiple tabs.
"""

import customtkinter as ctk
import tkinter as tk
from tkinter import messagebox, filedialog
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from co_adding_flux import coadd_spectra
from matplotlib.ticker import AutoMinorLocator
import numpy as np
import json
from main_run import get_species_data, astrovoigtfit_run
from edibles.utils.edibles_oracle import EdiblesOracle
from flux_wave_find import wave_flux_data
from main import master_function

# Speed of light in km/s
c_light = 299792.458

# Set appearance mode and default color theme
ctk.set_appearance_mode("Dark")
ctk.set_default_color_theme("blue")

def style_axis(ax, theme="Dark", grid=True, label_color=None, remove_minor_grid=False):
    """Applies the requested astronomy plotting style to an axis."""
    
    # Determine colors based on theme
    if theme == "Light":
        bg_color = 'white'
        fg_color = 'black'
    else:
        bg_color = '#2b2b2b' # Dark gray for plot background
        fg_color = 'white'
        
    if label_color:
        fg_color = label_color

    # Major + minor ticks on all sides
    ax.tick_params(
        which="both",
        direction="in",
        top=True,
        right=True,
        colors=fg_color,
        labelcolor=fg_color
    )
    
    # Minor tick locators
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    
    # Tick lengths
    ax.tick_params(which="major", length=6)
    ax.tick_params(which="minor", length=3)
    
    # Grid
    if grid:
        ax.grid(True, which="major", linestyle="--", alpha=0.6, color=fg_color)
        if not remove_minor_grid:
            ax.grid(True, which="minor", linestyle=":", alpha=0.3, color=fg_color)
    else:
        ax.grid(False)
    
    # Spines color
    for spine in ax.spines.values():
        spine.set_color(fg_color)
        
    # Labels and Title Color
    ax.xaxis.label.set_color(fg_color)
    ax.yaxis.label.set_color(fg_color)
    ax.title.set_color(fg_color)

def format_key_results(fitresult, molecules, star_name=None):
    """Formats key fit parameters (b, N, v_rad) for display/copy."""
    output = []
    if star_name:
        output.append(f"Star: {star_name}")
        output.append("-" * 20)
        
    params = fitresult.params
    
    for i, mol in enumerate(molecules):
        output.append(f"Species {i+1}: {mol}")
        
        comp_idx = 0
        while True:
            b_key = f"b_{i}_{comp_idx}"
            if b_key not in params:
                break
                
            b = params[f"b_{i}_{comp_idx}"]
            N = params[f"N_{i}_{comp_idx}"]
            v = params[f"v_rad_{i}_{comp_idx}"]
            
            def fmt_param(p):
                val = p.value
                err = p.stderr
                if err is None:
                    return f"{val:.4g} (fixed)"
                percent = abs(err / val) * 100 if val != 0 else 0
                return f"{val:.4g} +/- {err:.4g} ({percent:.2f}%)"

            output.append(f"  Component {comp_idx+1}:")
            output.append(f"    b:     {fmt_param(b)}")
            output.append(f"    N:     {fmt_param(N)}")
            output.append(f"    v_rad: {fmt_param(v)}")
            
            comp_idx += 1
        output.append("") 

    output.append("-" * 30)
    output.append(f"Chi-Square:         {fitresult.chisqr:.4f}")
    output.append(f"Reduced Chi-Square: {fitresult.redchi:.4f}")
    
    return "\n".join(output)

class AnalysisTab(ctk.CTkFrame):
    def __init__(self, parent):
        super().__init__(parent)
        
        # Load Species Data
        self.species_data = get_species_data()
        self.species_list = list(self.species_data.keys())
        
        # Data storage for plots
        self.plot_data = {}
        self.key_results_text = ""
        self.available_files = [] 
        self.general_data = None
        self.general_file_path = ""

        # --- Paned Window for Resizable Sidebar ---
        self.paned_window = tk.PanedWindow(self, orient=tk.HORIZONTAL, sashwidth=4, bg="#2b2b2b", bd=0)
        self.paned_window.pack(fill="both", expand=True)

        # --- Left Sidebar (Tabview) ---
        self.sidebar_frame = ctk.CTkFrame(self.paned_window, width=340, corner_radius=0)
        self.sidebar_frame.grid_rowconfigure(0, weight=1)
        self.sidebar_frame.grid_columnconfigure(0, weight=1)
        
        self.sidebar_tabs = ctk.CTkTabview(self.sidebar_frame)
        self.sidebar_tabs.grid(row=0, column=0, padx=10, pady=10, sticky="nsew")
        self.sidebar_tabs.add("Inputs")
        self.sidebar_tabs.add("Plot Settings")
        
        self.setup_inputs_tab()
        self.setup_settings_tab()
        
        # Add sidebar to paned window
        self.paned_window.add(self.sidebar_frame, minsize=250)

        # --- Right Main Area ---
        self.main_frame = ctk.CTkFrame(self.paned_window, corner_radius=0, fg_color="transparent")
        self.main_frame.grid_rowconfigure(1, weight=1)
        self.main_frame.grid_columnconfigure(0, weight=1)
        
        # Add main frame to paned window
        self.paned_window.add(self.main_frame, stretch="always")

        # 1. Top: File Selection
        self.setup_file_selection(self.main_frame)

        # 2. Middle: Plots
        self.tabview = ctk.CTkTabview(self.main_frame)
        self.tabview.grid(row=1, column=0, padx=20, pady=(0, 10), sticky="nsew")
        self.tabview.add("Continuum Fit")
        self.tabview.add("Voigt Fit")
        self.tabview.add("Fit Results")
        
        self.canvas1 = None
        self.toolbar1 = None
        self.canvas2 = None
        self.toolbar2 = None
        self.setup_results_tab()

        # 3. Bottom: Fit Parameters
        self.setup_fit_params(self.main_frame)

    def setup_inputs_tab(self):
        tab = self.sidebar_tabs.tab("Inputs")
        tab.grid_columnconfigure(0, weight=1)
        tab.grid_rowconfigure(6, weight=1) 

        # Session Management
        session_frame = ctk.CTkFrame(tab, fg_color="transparent")
        session_frame.grid(row=0, column=0, padx=10, pady=(10, 5), sticky="ew")
        ctk.CTkButton(session_frame, text="💾 Save", width=60, command=self.save_session).pack(side="left", padx=(0, 5), expand=True, fill="x")
        ctk.CTkButton(session_frame, text="📂 Load", width=60, command=self.load_session).pack(side="right", expand=True, fill="x")

        ctk.CTkLabel(tab, text="Configuration", font=ctk.CTkFont(size=16, weight="bold")).grid(row=1, column=0, pady=(10, 5))
        
        # Mode dropdown. Right now only Local Mode actually works — Global Mode
        # is here as a stub so whoever adds it next has an obvious place to plug in.
        # When you add it: handle it in on_mode_change below and add whatever
        # extra widgets that mode needs.
        self.mode_var = ctk.StringVar(value="Local Mode")
        self.mode_selector = ctk.CTkOptionMenu(
            tab,
            values=["Local Mode", "Global Mode"],
            variable=self.mode_var,
            command=self.on_mode_change,
        )
        self.mode_selector.grid(row=2, column=0, padx=10, pady=(0, 10), sticky="ew")

        # Star Name
        self.star_label = ctk.CTkLabel(tab, text="Star Name:", anchor="w")
        self.star_label.grid(row=3, column=0, padx=10, sticky="w")
        self.star_entry = ctk.CTkEntry(tab)
        self.star_entry.grid(row=4, column=0, padx=10, pady=(0, 5), sticky="ew")
        self.star_entry.insert(0, "HD 183143")

        # Molecules List
        ctk.CTkLabel(tab, text="Species:", anchor="w").grid(row=5, column=0, padx=10, pady=(10, 0), sticky="w")
        self.mol_scroll = ctk.CTkScrollableFrame(tab, height=300)
        self.mol_scroll.grid(row=6, column=0, padx=10, pady=(0, 10), sticky="nsew")
        
        self.molecule_entries = []
        # Initial molecules
        self.add_molecule_inputs("13CH+_4032", "2, 2", f"{1e13/70}, {1e13/70}", "-11, 4")
        self.add_molecule_inputs("12CH+_4032", "2, 2", "1e13, 1e13", "-11, 4")
        
        # Local Mode Inputs (Knots)
        self.n_knots_label = ctk.CTkLabel(self.mol_scroll, text="N Knots:", anchor="w")
        self.n_knots_entry = ctk.CTkEntry(self.mol_scroll, placeholder_text="10")
        self.n_knots_entry.insert(0, "10")
        self.n_knots_label.pack(pady=(10,0), anchor="w", padx=5)
        self.n_knots_entry.pack(pady=(0,5), fill="x", padx=5)
        
        # Continuum Order
        self.cont_order_label = ctk.CTkLabel(self.mol_scroll, text="Cont. Order:", anchor="w")
        self.cont_order_combo = ctk.CTkComboBox(self.mol_scroll, values=["Linear", "Quadratic", "Cubic"])
        self.cont_order_combo.set("Cubic")
        self.cont_order_label.pack(pady=(5,0), anchor="w", padx=5)
        self.cont_order_combo.pack(pady=(0,5), fill="x", padx=5)
        
        self.knots_x_label = ctk.CTkLabel(self.mol_scroll, text="Knots X Array (optional):", anchor="w")
        self.knots_x_entry = ctk.CTkEntry(self.mol_scroll, placeholder_text="e.g. 6707.5, 6708.0")
        self.knots_x_label.pack(pady=(5,0), anchor="w", padx=5)
        self.knots_x_entry.pack(pady=(0,5), fill="x", padx=5)
        
        self.show_knots_switch = ctk.CTkSwitch(self.mol_scroll, text="Show Knots")
        self.show_knots_switch.select()
        self.show_knots_switch.pack(pady=(5,5), anchor="w", padx=5)
        
        self.show_continuum_switch = ctk.CTkSwitch(self.mol_scroll, text="Show Continuum")
        self.show_continuum_switch.select()
        self.show_continuum_switch.pack(pady=(0,5), anchor="w", padx=5)
        
        self.normalize_plot_switch = ctk.CTkSwitch(self.mol_scroll, text="Plot Normalized")
        self.normalize_plot_switch.deselect()
        self.normalize_plot_switch.pack(pady=(0,10), anchor="w", padx=5)

        # Buttons
        btn_frame = ctk.CTkFrame(tab, fg_color="transparent")
        btn_frame.grid(row=7, column=0, padx=10, pady=5, sticky="ew")
        ctk.CTkButton(btn_frame, text="+ Add Species", width=80, command=lambda: self.add_molecule_inputs("", "", "", "")).pack(side="left", padx=(0, 5), expand=True, fill="x")
        ctk.CTkButton(btn_frame, text="- Remove Species", width=80, fg_color="#D32F2F", hover_color="#B71C1C", command=self.remove_last_molecule).pack(side="right", expand=True, fill="x")

        ctk.CTkButton(tab, text="Run Analysis", command=self.run_analysis, height=40, font=ctk.CTkFont(size=14, weight="bold")).grid(row=8, column=0, padx=10, pady=20, sticky="ew")

    def setup_file_selection(self, parent):
        self.file_frame = ctk.CTkFrame(parent)
        self.file_frame.grid(row=0, column=0, padx=20, pady=10, sticky="ew")
        
        # Local Mode Widgets
        self.local_fetch_btn = ctk.CTkButton(self.file_frame, text="🔍 Fetch Files", width=100, command=self.fetch_files_local)
        self.local_fetch_btn.pack(side="left", padx=10)
        
        self.local_select_btn = ctk.CTkButton(self.file_frame, text="Select Files ⬇", width=100, command=self.open_file_selector)
        self.local_select_btn.pack(side="left", padx=10)
        
        self.local_coadd_btn = ctk.CTkButton(self.file_frame, text="Co-add Selected", width=100, fg_color="green", command=self.coadd_spectra_action)
        self.local_coadd_btn.pack(side="left", padx=10)
        
        self.local_status_label = ctk.CTkLabel(self.file_frame, text="No files fetched", font=("Arial", 12), text_color="gray")
        self.local_status_label.pack(side="left", padx=10)
        
        self.local_sigma_frame = ctk.CTkFrame(self.file_frame, fg_color="transparent")
        self.local_sigma_frame.pack(side="left", padx=10)
        
        # Two sigma fields — they look similar but feed different things:
        # "Sigma" goes into coadd_spectra as the per-spectrum noise (used to
        # weight each spectrum during co-addition).
        ctk.CTkLabel(self.local_sigma_frame, text="Sigma:", width=50).pack(side="left")
        self.local_sigma_entry = ctk.CTkEntry(self.local_sigma_frame, width=100)
        self.local_sigma_entry.pack(side="left")
        self.local_sigma_entry.insert(0, "0.002")

        # "Fit σ" is the noise of the *co-added* spectrum and gets passed to
        # the fit as std_dev (lmfit uses weights = 1/std_dev). Default matches
        # Sigma; bump it if your co-added spectrum is noisier than that.
        ctk.CTkLabel(self.local_sigma_frame, text="Fit σ:", width=50).pack(side="left", padx=(10, 0))
        self.local_std_dev_entry = ctk.CTkEntry(self.local_sigma_frame, width=100)
        self.local_std_dev_entry.pack(side="left")
        self.local_std_dev_entry.insert(0, "0.002")

        self.local_files = []
        self.local_selected_vars = []

    def setup_fit_params(self, parent):
        # This frame now only contains Wavelength and Absorption Range, Degree/Knots are in setup_inputs_tab
        frame = ctk.CTkFrame(parent)
        frame.grid(row=2, column=0, padx=20, pady=10, sticky="ew")
        
        # Wavelength Range
        ctk.CTkLabel(frame, text="Fit Range (Min, Max):").pack(side="left", padx=(10, 5))
        self.wave_min = ctk.CTkEntry(frame, width=80)
        self.wave_min.pack(side="left", padx=5)
        self.wave_min.insert(0, "4231.5")
        self.wave_max = ctk.CTkEntry(frame, width=80)
        self.wave_max.pack(side="left", padx=5)
        self.wave_max.insert(0, "4233.5")
        
        # Absorption Range
        ctk.CTkLabel(frame, text="Abs. Range (Min, Max):").pack(side="left", padx=(20, 5))
        self.abs_min = ctk.CTkEntry(frame, width=80)
        self.abs_min.pack(side="left", padx=5)
        self.abs_min.insert(0, "4232.05")
        self.abs_max = ctk.CTkEntry(frame, width=80)
        self.abs_max.pack(side="left", padx=5)
        self.abs_max.insert(0, "4232.8")
        
        # Degree (This is now handled in setup_inputs_tab, so removing from here)
        # ctk.CTkLabel(frame, text="Poly. Degree:").pack(side="left", padx=(20, 5))
        # self.degree_combo = ctk.CTkComboBox(frame, values=["0", "1", "2", "3", "4"], width=60)
        # self.degree_combo.pack(side="left", padx=5)
        # self.degree_combo.set("3")

    def setup_settings_tab(self):
        tab = self.sidebar_tabs.tab("Plot Settings")
        
        # Create Scrollable Frame
        scroll_frame = ctk.CTkScrollableFrame(tab, fg_color="transparent")
        scroll_frame.pack(fill="both", expand=True)
        scroll_frame.grid_columnconfigure(0, weight=1)
        
        ctk.CTkLabel(scroll_frame, text="Appearance", font=ctk.CTkFont(size=14, weight="bold")).grid(row=0, column=0, pady=(10, 5))
        
        # Theme
        ctk.CTkLabel(scroll_frame, text="Theme:", anchor="w").grid(row=1, column=0, padx=10, sticky="w")
        self.theme_var = ctk.StringVar(value="Light")
        self.theme_switch = ctk.CTkSegmentedButton(scroll_frame, values=["Dark", "Light"], variable=self.theme_var, command=self.on_theme_change)
        self.theme_switch.grid(row=2, column=0, padx=10, pady=(0, 10), sticky="ew")
        
        # Show Grid Toggle
        self.grid_var = ctk.BooleanVar(value=True)
        ctk.CTkSwitch(scroll_frame, text="Show Grid", variable=self.grid_var, command=self.update_plots).grid(row=4, column=0, padx=10, pady=5, sticky="w")
        
        # Remove Minor Grid Toggle
        self.remove_minor_grid_var = ctk.BooleanVar(value=False)
        ctk.CTkSwitch(scroll_frame, text="Remove Minor Grid", variable=self.remove_minor_grid_var, command=self.update_plots).grid(row=5, column=0, padx=10, pady=5, sticky="w")
        
        # Plot Style
        ctk.CTkLabel(scroll_frame, text="Plot Style:", anchor="w").grid(row=6, column=0, padx=10, sticky="w")
        self.plot_style_combo = ctk.CTkComboBox(scroll_frame, values=["Line", "Boxy", "Scatter"], command=self.update_plots)
        self.plot_style_combo.set("Line")
        self.plot_style_combo.grid(row=7, column=0, padx=10, pady=(0, 5), sticky="ew")
        
        # Individual Lines Toggle
        self.show_individual_lines_var = ctk.BooleanVar(value=False)
        ctk.CTkSwitch(scroll_frame, text="Show Individual Lines", variable=self.show_individual_lines_var, command=self.update_plots).grid(row=8, column=0, padx=10, pady=5, sticky="w")
        
        # Component-wise Lines Toggle
        self.show_comp_lines_var = ctk.BooleanVar(value=False)
        ctk.CTkSwitch(scroll_frame, text="Show Component Lines", variable=self.show_comp_lines_var, command=self.update_plots).grid(row=9, column=0, padx=10, pady=5, sticky="w")
        
        # Show Total Fit Toggle
        self.show_total_fit_switch = ctk.BooleanVar(value=True)
        ctk.CTkSwitch(scroll_frame, text="Show Total Fit", variable=self.show_total_fit_switch, command=self.update_plots).grid(row=10, column=0, padx=10, pady=5, sticky="w")
        
        # Show Velocity Axis Toggle
        self.show_velocity_axis_var = ctk.BooleanVar(value=False)
        ctk.CTkSwitch(scroll_frame, text="Show Velocity Axis (km/s)", variable=self.show_velocity_axis_var, command=self.update_plots).grid(row=11, column=0, padx=10, pady=5, sticky="w")
        
        # Show Sigma Toggle
        self.show_sigma_var = ctk.BooleanVar(value=False)
        ctk.CTkSwitch(scroll_frame, text="Show Residual Sigma", variable=self.show_sigma_var, command=self.update_plots).grid(row=12, column=0, padx=10, pady=5, sticky="w")
        
        # Sigma Level Slider
        ctk.CTkLabel(scroll_frame, text="Sigma Level (Multiplier):", anchor="w").grid(row=13, column=0, padx=10, pady=(5, 0), sticky="w")
        self.sigma_level_slider = ctk.CTkSlider(scroll_frame, from_=1.0, to=3.0, number_of_steps=20, command=self.update_plots)
        self.sigma_level_slider.set(1.0) 
        self.sigma_level_slider.grid(row=14, column=0, padx=10, pady=(0, 10), sticky="ew")
        
        # Individual Spectra Alpha (for GM3)
        ctk.CTkLabel(scroll_frame, text="Individual Spectra Alpha:", anchor="w").grid(row=15, column=0, padx=10, pady=(10, 0), sticky="w")
        self.spectra_alpha_slider = ctk.CTkSlider(scroll_frame, from_=0.0, to=1.0, number_of_steps=100, command=self.update_plots)
        self.spectra_alpha_slider.set(0.8) # Default alpha (High visibility)
        self.spectra_alpha_slider.grid(row=16, column=0, padx=10, pady=(0, 10), sticky="ew")

        ctk.CTkLabel(scroll_frame, text="Labels", font=ctk.CTkFont(size=14, weight="bold")).grid(row=17, column=0, pady=(15, 5))
        
        # Custom Labels
        ctk.CTkLabel(scroll_frame, text="Plot Title:", anchor="w").grid(row=18, column=0, padx=10, sticky="w")
        self.title_entry = ctk.CTkEntry(scroll_frame, placeholder_text="Auto")
        self.title_entry.grid(row=19, column=0, padx=10, pady=(0, 5), sticky="ew")
        
        ctk.CTkLabel(scroll_frame, text="X-Axis Label:", anchor="w").grid(row=20, column=0, padx=10, sticky="w")
        self.xlabel_entry = ctk.CTkEntry(scroll_frame, placeholder_text="Wavelength (Å)")
        self.xlabel_entry.grid(row=21, column=0, padx=10, pady=(0, 5), sticky="ew")
        
        ctk.CTkLabel(scroll_frame, text="Y-Axis Label:", anchor="w").grid(row=22, column=0, padx=10, sticky="w")
        self.ylabel_entry = ctk.CTkEntry(scroll_frame, placeholder_text="Normalized Flux")
        self.ylabel_entry.grid(row=23, column=0, padx=10, pady=(0, 5), sticky="ew")
        
        # Colors
        colors = ['white', 'black', 'gray', 'red', 'blue', 'green', 'cyan', 'magenta', 'yellow', '#7B2CBF']
        
        ctk.CTkLabel(scroll_frame, text="Label Color:", anchor="w").grid(row=24, column=0, padx=10, sticky="w")
        self.color_entry = ctk.CTkComboBox(scroll_frame, values=colors)
        self.color_entry.grid(row=25, column=0, padx=10, pady=(0, 5), sticky="ew")
        self.color_entry.set("black") # Changed default to black for Light theme
        
        ctk.CTkLabel(scroll_frame, text="Continuum Fit Color:", anchor="w").grid(row=26, column=0, padx=10, sticky="w")
        self.cont_color_combo = ctk.CTkComboBox(scroll_frame, values=colors)
        self.cont_color_combo.grid(row=27, column=0, padx=10, pady=(0, 5), sticky="ew")
        self.cont_color_combo.set("red")
        
        ctk.CTkLabel(scroll_frame, text="Voigt Fit Color:", anchor="w").grid(row=28, column=0, padx=10, sticky="w")
        self.voigt_color_combo = ctk.CTkComboBox(scroll_frame, values=colors)
        self.voigt_color_combo.grid(row=29, column=0, padx=10, pady=(0, 10), sticky="ew")
        self.voigt_color_combo.set("green")
        
        ctk.CTkButton(scroll_frame, text="🔄 Update Plots", command=self.update_plots).grid(row=30, column=0, padx=10, pady=20, sticky="ew")

    def on_theme_change(self, value):
        """Automatically updates label color based on theme."""
        if value == "Light":
            self.color_entry.set("black")
        else:
            self.color_entry.set("white")

    def on_mode_change(self, value):
        """Fires when the user picks something from the mode dropdown.

        Only Local Mode does anything for now. If someone picks Global Mode we
        just pop a message and flip back. When Global Mode gets implemented,
        you'll probably want to:
          - add its widgets somewhere in setup_inputs_tab / setup_file_selection
          - branch on self.mode_var.get() in run_analysis and coadd_spectra_action
          - tag self.plot_data['mode'] so update_plots can do the right thing
        """
        if value == "Global Mode":
            messagebox.showinfo(
                "Not implemented",
                "Global Mode isn't wired up yet — sticking with Local Mode."
            )
            self.mode_var.set("Local Mode")

    def setup_results_tab(self):
        tab = self.tabview.tab("Fit Results")
        
        # Copy Button Frame
        btn_frame = ctk.CTkFrame(tab, fg_color="transparent")
        btn_frame.pack(fill="x", padx=10, pady=5)
        
        ctk.CTkButton(btn_frame, text="💾 Save Key Results", command=self.save_key_results_to_file, width=150).pack(side="left")
        ctk.CTkButton(btn_frame, text="📋 Copy Key Results", command=self.copy_key_results, width=150).pack(side="right")
        
        # Scrollable Frame for Results
        self.results_scroll = ctk.CTkScrollableFrame(tab)
        self.results_scroll.pack(fill="both", expand=True, padx=10, pady=(0, 10))

    def create_result_card(self, parent, molecule_name, components_data):
        """Creates a styled card for a molecule's results."""
        card = ctk.CTkFrame(parent, border_width=1, border_color="gray")
        card.pack(fill="x", pady=10, padx=5)
        
        # Header
        ctk.CTkLabel(card, text=molecule_name, font=("Arial", 16, "bold")).pack(anchor="w", padx=10, pady=(10, 5))
        
        # Grid Frame
        grid = ctk.CTkFrame(card, fg_color="transparent")
        grid.pack(fill="x", padx=10, pady=5)
        
        # Headers
        headers = ["Param", "Initial", "Best Fit", "Error", "% Uncert"]
        for col, h in enumerate(headers):
            ctk.CTkLabel(grid, text=h, font=("Arial", 12, "bold"), text_color="gray").grid(row=0, column=col, padx=10, pady=5, sticky="w")
            
        # Data Rows
        row = 1
        for comp in components_data:
            ctk.CTkLabel(grid, text=f"Comp {comp['id']}", font=("Arial", 12, "bold")).grid(row=row, column=0, columnspan=5, sticky="w", pady=(5,0))
            row += 1
            
            for param in ['b', 'N', 'v_rad']:
                p_data = comp[param]
                ctk.CTkLabel(grid, text=param).grid(row=row, column=0, padx=10, sticky="w")
                ctk.CTkLabel(grid, text=f"{p_data['init']:.4g}").grid(row=row, column=1, padx=10, sticky="w")
                ctk.CTkLabel(grid, text=f"{p_data['value']:.4g}", text_color="#4CAF50" if p_data['stderr'] else "white").grid(row=row, column=2, padx=10, sticky="w")
                
                err_text = f"{p_data['stderr']:.4g}" if p_data['stderr'] else "fixed"
                ctk.CTkLabel(grid, text=err_text).grid(row=row, column=3, padx=10, sticky="w")
                
                pct_text = f"{p_data['percent']:.2f}%" if p_data['stderr'] else "-"
                ctk.CTkLabel(grid, text=pct_text).grid(row=row, column=4, padx=10, sticky="w")
                
                row += 1

    def display_results(self, fitresult, molecules):
        # Clear previous results
        for widget in self.results_scroll.winfo_children():
            widget.destroy()
            
        params = fitresult.params
        
        # 1. Statistics Card
        stats_card = ctk.CTkFrame(self.results_scroll, border_width=1, border_color="gray")
        stats_card.pack(fill="x", pady=10, padx=5)
        ctk.CTkLabel(stats_card, text="Fit Statistics", font=("Arial", 14, "bold")).pack(anchor="w", padx=10, pady=(10, 5))
        
        stats_grid = ctk.CTkFrame(stats_card, fg_color="transparent")
        stats_grid.pack(fill="x", padx=10, pady=5)
        ctk.CTkLabel(stats_grid, text=f"Chi-Square: {fitresult.chisqr:.4f}").pack(side="left", padx=20)
        ctk.CTkLabel(stats_grid, text=f"Reduced Chi-Square: {fitresult.redchi:.4f}").pack(side="left", padx=20)

        # 2. Molecule Cards
        for i, mol in enumerate(molecules):
            components_data = []
            comp_idx = 0
            while True:
                b_key = f"b_{i}_{comp_idx}"
                if b_key not in params:
                    break
                
                comp_data = {'id': comp_idx + 1}
                for p_name in ['b', 'N', 'v_rad']:
                    p = params[f"{p_name}_{i}_{comp_idx}"]
                    comp_data[p_name] = {
                        'init': p.init_value,
                        'value': p.value,
                        'stderr': p.stderr,
                        'percent': abs(p.stderr / p.value) * 100 if p.stderr and p.value != 0 else 0
                    }
                components_data.append(comp_data)
                comp_idx += 1
            
            self.create_result_card(self.results_scroll, f"Species {i+1}: {mol}", components_data)

        # 3. Full Report (Expandable)
        ctk.CTkLabel(self.results_scroll, text="Full Report", font=("Arial", 14, "bold")).pack(anchor="w", padx=10, pady=(20, 5))
        report_box = ctk.CTkTextbox(self.results_scroll, height=200, font=("Courier", 12))
        report_box.pack(fill="x", padx=5, pady=5)
        report_box.insert("1.0", fitresult.fit_report())

    def copy_key_results(self):
        if self.key_results_text:
            self.clipboard_clear()
            self.clipboard_append(self.key_results_text)
            messagebox.showinfo("Copied", "Key results copied to clipboard!")
        else:
            messagebox.showwarning("Empty", "No results to copy yet.")

    def save_key_results_to_file(self):
        if not self.key_results_text:
            messagebox.showwarning("Empty", "No results to save yet.")
            return
            
        file_path = filedialog.asksaveasfilename(
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")],
            title="Save Fit Results"
        )
        if file_path:
            try:
                with open(file_path, 'w') as f:
                    f.write(self.key_results_text)
                messagebox.showinfo("Saved", "Results saved successfully!")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to save results: {e}")

    def add_molecule_inputs(self, name_val, b_val, N_val, v_val, tie_b_val="Custom", tie_v_val="Custom"):
        index = len(self.molecule_entries)
        frame = ctk.CTkFrame(self.mol_scroll)
        frame.pack(fill="x", pady=5)
        
        ctk.CTkLabel(frame, text=f"Species {index+1}", font=("Arial", 12, "bold")).pack(anchor="w", padx=5)
        
        # Species Dropdown
        name = ctk.CTkComboBox(frame, values=self.species_list, command=lambda choice: self.update_fit_range(choice))
        name.pack(fill="x", padx=5, pady=2)
        if name_val in self.species_list:
            name.set(name_val)
        else:
            name.set(name_val if name_val else self.species_list[0] if self.species_list else "")
        
        # --- b parameter ---
        b_frame = ctk.CTkFrame(frame, fg_color="transparent")
        b_frame.pack(fill="x", padx=5, pady=2)
        
        b = ctk.CTkEntry(b_frame, placeholder_text="b-value (comma sep)")
        b.pack(side="left", fill="x", expand=True)
        b.insert(0, b_val)
        
        b_tie_combo = None
        if index > 0:
            b_tie_combo = ctk.CTkComboBox(b_frame, values=["Custom", "Tie to Species 1"], width=110,
                                          command=lambda choice, entry=b: self.toggle_entry(entry, choice))
            b_tie_combo.pack(side="right", padx=(5, 0))
            b_tie_combo.set(tie_b_val)
            self.toggle_entry(b, tie_b_val) # Set initial state

        # --- N parameter ---
        N = ctk.CTkEntry(frame, placeholder_text="column density (comma sep)")
        N.pack(fill="x", padx=5, pady=2)
        N.insert(0, N_val)
        
        # --- v_rad parameter ---
        v_frame = ctk.CTkFrame(frame, fg_color="transparent")
        v_frame.pack(fill="x", padx=5, pady=2)
        
        v = ctk.CTkEntry(v_frame, placeholder_text="radial velocity (comma sep)")
        v.pack(side="left", fill="x", expand=True)
        v.insert(0, v_val)
        
        v_tie_combo = None
        if index > 0:
            v_tie_combo = ctk.CTkComboBox(v_frame, values=["Custom", "Tie to Species 1"], width=110,
                                          command=lambda choice, entry=v: self.toggle_entry(entry, choice))
            v_tie_combo.pack(side="right", padx=(5, 0))
            v_tie_combo.set(tie_v_val)
            self.toggle_entry(v, tie_v_val) # Set initial state
        
        self.molecule_entries.append({
            'frame': frame, 
            'name': name, 
            'b': b, 
            'N': N, 
            'v': v,
            'b_tie': b_tie_combo,
            'v_tie': v_tie_combo
        })

    def toggle_entry(self, entry_widget, choice):
        if choice == "Tie to Species 1":
            entry_widget.configure(state="disabled", fg_color="gray25")
        else:
            entry_widget.configure(state="normal", fg_color=["#F9F9FA", "#343638"]) # Default colors

    def remove_last_molecule(self):
        if not self.molecule_entries:
            return
        last_entry = self.molecule_entries.pop()
        last_entry['frame'].destroy()

    def update_fit_range(self, species_name):
        """Updates the fit range and absorption range based on the selected species."""
        if species_name in self.species_data:
            wrange = self.species_data[species_name]['wrange']
            
            # Update Fit Range
            self.wave_min.delete(0, tk.END)
            self.wave_min.insert(0, str(wrange[0]))
            self.wave_max.delete(0, tk.END)
            self.wave_max.insert(0, str(wrange[1]))
            
            # Update Absorption Range (wrange_min + 0.5, wrange_max - 0.5)
            abs_min_val = wrange[0] + 0.5
            abs_max_val = wrange[1] - 0.5
            
            self.abs_min.delete(0, tk.END)
            self.abs_min.insert(0, str(abs_min_val))
            self.abs_max.delete(0, tk.END)
            self.abs_max.insert(0, str(abs_max_val))

    def fetch_files_local(self):
        """Fetches files for GM3 and initializes selection variables."""
        star_name = self.star_entry.get()
        if not star_name:
            messagebox.showerror("Error", "Please enter a star name.")
            return

        try:
            # Use EdiblesOracle to find files
            # Determine wavelength from first molecule or default to 6708
            target_wave = 6708
            wrange = [6707, 6709] # Default
            
            if self.molecule_entries:
                mol_name = self.molecule_entries[0]['name'].get()
                if mol_name in self.species_data:
                    # Use center of wrange
                    wrange = self.species_data[mol_name]['wrange']
                    target_wave = float(self.species_data[mol_name]['line']) # Use the 'line' column value
            
            # Store for co-add
            self.local_target_wave = target_wave
            self.local_wrange = wrange
            
            oracle = EdiblesOracle()
            files = oracle.getFilteredObsList([star_name], Wave=target_wave, MergedOnly=False)
            
            # Convert to list if it's a pandas Series/DataFrame or numpy array
            if hasattr(files, 'tolist'):
                files = files.tolist()
            
            self.local_files = []
            self.local_selected_vars = []
            
            if files and len(files) > 0:
                self.local_files = files
                # Initialize BooleanVars for each file
                self.local_selected_vars = [ctk.BooleanVar(value=False) for _ in files]
                
                self.local_status_label.configure(text=f"Found {len(files)} files. Click 'Select Files'.", text_color="black")
            else:
                self.local_status_label.configure(text="No files found.", text_color="red")
                
        except Exception as e:
            messagebox.showerror("Error", f"Failed to fetch files: {e}")

    def open_file_selector(self):
        """Opens a popup to select files."""
        if not self.local_files:
            messagebox.showwarning("Warning", "Please fetch files first.")
            return
            
        # Create Toplevel Window
        top = ctk.CTkToplevel(self)
        top.title("Select Files to Co-add")
        top.geometry("400x400")
        top.attributes("-topmost", True)
        
        ctk.CTkLabel(top, text=f"Select Files for {self.star_entry.get()}", font=("Arial", 14, "bold")).pack(pady=10)
        
        # Scrollable Frame
        scroll = ctk.CTkScrollableFrame(top)
        scroll.pack(fill="both", expand=True, padx=10, pady=5)
        
        for i, f in enumerate(self.local_files):
            cb = ctk.CTkCheckBox(scroll, text=f"{i}: {f}", variable=self.local_selected_vars[i])
            cb.pack(anchor="w", padx=5, pady=2)
            
        ctk.CTkButton(top, text="Done", command=top.destroy).pack(pady=10)

    def coadd_spectra_action(self):
        """Co-adds selected spectra for GM3."""
        if not self.local_files:
             messagebox.showwarning("Warning", "No files available.")
             return

        selected_indices = [i for i, var in enumerate(self.local_selected_vars) if var.get()]
        selected_files = [self.local_files[i] for i in selected_indices]
        
        if not selected_files:
            messagebox.showwarning("Warning", "Please select at least one file to co-add.")
            return
            
        try:
            wave_list = []
            flux_list = []
            
            # Get wavelength range from inputs with buffer
            try:
                w_min = float(self.wave_min.get())
                w_max = float(self.wave_max.get())
                wrange_input = [w_min, w_max]
            except ValueError:
                wrange_input = self.local_wrange if hasattr(self, 'local_wrange') else [6700, 6715]

            star_name = self.star_entry.get()
            target_wave = self.local_target_wave if hasattr(self, 'local_target_wave') else 6708
            
            for i in selected_indices:
                filename = self.local_files[i]
                try:
                    # Use wave_flux_data as requested
                    # wave0, flux0 = wave_flux_data(star_name="HD 147889", wave = 6707,wrange=[6705, 6710],file_number=0)
                    wave, flux = wave_flux_data(star_name=star_name, wave=target_wave, wrange=wrange_input, file_name=filename)
                    
                    if wave is None or len(wave) == 0:
                        print(f"Warning: No data found for {filename} in range {wrange_input}")
                        continue
                        
                    wave_list.append(wave)
                    flux_list.append(flux)
                    
                except Exception as e:
                    print(f"Error getting data for {filename}: {e}")
                    continue
                
            if not wave_list:
                messagebox.showerror("Error", f"No valid spectra found in the selected files for the {w_min}-{w_max} A range.")
                return
                
            # Co-add
            # Parse Sigma Input
            sigma_input = self.local_sigma_entry.get().strip()
            sigma_list = 0.002 # Default
            
            if sigma_input:
                try:
                    # Check for comma-separated list
                    if "," in sigma_input:
                        parts = [float(x.strip()) for x in sigma_input.split(",")]
                        if len(parts) == 1:
                            sigma_list = parts[0]
                        elif len(parts) == len(selected_indices):
                            sigma_list = parts
                        else:
                            messagebox.showerror("Error", f"Number of sigma values ({len(parts)}) must match selected files ({len(selected_indices)}) or be a single value.")
                            return
                    else:
                        sigma_list = float(sigma_input)
                except ValueError:
                    messagebox.showerror("Error", "Invalid Sigma input. Must be a float or comma-separated floats.")
                    return
            
            wave_coadd, flux_coadd, sigma_coadd = coadd_spectra(wave_list, flux_list, sigma_list=sigma_list)
            
            # Retrieve rest wavelength (lambda_0) for the first species
            lambda_0 = None
            if self.molecule_entries:
                mol_name = self.molecule_entries[0]['name'].get()
                if mol_name:
                    all_species = get_species_data('species.txt')
                    if mol_name in all_species:
                        # 'line' is the rest wavelength (string or float in dict)
                        try:
                            lambda_0 = float(all_species[mol_name]['line'])
                        except (ValueError, TypeError):
                            pass

            # Initial Plotting
            self.plot_data = { # Initialize
                'wave': wave_coadd,
                'flux': flux_coadd,
                'continuum': np.ones_like(wave_coadd), # Initial dummy continuum
                'norm_flux': flux_coadd, # Initially same
                'model': np.zeros_like(wave_coadd),
                'residuals': np.zeros_like(wave_coadd),
                'mode': 'Local Mode', # Mark as GM3
                'sigma': sigma_coadd,
                'individual_spectra': list(zip(wave_list, flux_list)), # Store individual spectra
                'lambda_0_ref': lambda_0 # Store lambda_0 for velocity conversion
            }
            
            # Plot in Continuum Panel
            self.update_plots(plot_title=f"Co-added Spectrum ({len(selected_files)} files)")
            self.local_status_label.configure(text=f"Co-added {len(selected_files)} files.", text_color="green")
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to co-add spectra: {e}")

    def save_session(self):
        """Saves current inputs to a JSON file."""
        data = {
            "star": self.star_entry.get(),
            "wave_min": self.wave_min.get(),
            "wave_max": self.wave_max.get(),
            "abs_min": self.abs_min.get(),
            "abs_max": self.abs_max.get(),
            "n_knots": self.n_knots_entry.get(),
            "cont_order": self.cont_order_combo.get(),
            "knots_x": self.knots_x_entry.get(),
            "show_knots": self.show_knots_switch.get(),
            "molecules": [
                {
                    "name": entry['name'].get(),
                    "b": entry['b'].get(),
                    "N": entry['N'].get(),
                    "v": entry['v'].get(),
                    "tie_b": entry['b_tie'].get() if entry.get('b_tie') else "Custom",
                    "tie_v": entry['v_tie'].get() if entry.get('v_tie') else "Custom"
                }
                for entry in self.molecule_entries
            ]
        }
        
        data["mode"] = "Local Mode"
        
        # Save Selected Files
        if hasattr(self, 'local_files'):
            data["local_files"] = self.local_files
            data["local_selected_indices"] = [i for i, var in enumerate(self.local_selected_vars) if var.get()]
            if hasattr(self, 'local_target_wave'):
                data["local_target_wave"] = self.local_target_wave
            if hasattr(self, 'local_wrange'):
                data["local_wrange"] = self.local_wrange
            if hasattr(self, 'plot_data') and 'lambda_0_ref' in self.plot_data:
                data["local_lambda_0_ref"] = self.plot_data['lambda_0_ref']
        
        # Add Plotting Style Settings
        data["theme"] = self.theme_var.get()
        data["show_grid"] = self.grid_var.get()
        data["remove_minor_grid"] = self.remove_minor_grid_var.get()
        data["plot_style"] = self.plot_style_combo.get()
        data["show_individual_lines"] = self.show_individual_lines_var.get()
        data["show_comp_lines"] = self.show_comp_lines_var.get()
        data["show_total_fit"] = self.show_total_fit_switch.get()
        data["plot_title"] = self.title_entry.get()
        data["xlabel"] = self.xlabel_entry.get()
        data["ylabel"] = self.ylabel_entry.get()
        data["label_color"] = self.color_entry.get()
        data["cont_color"] = self.cont_color_combo.get()
        data["voigt_color"] = self.voigt_color_combo.get()
        data["show_velocity_axis"] = self.show_velocity_axis_var.get()
        data["show_sigma"] = self.show_sigma_var.get()
        data["sigma_level"] = self.sigma_level_slider.get()
        
        file_path = filedialog.asksaveasfilename(defaultextension=".json", filetypes=[("JSON files", "*.json")])
        if file_path:
            try:
                with open(file_path, 'w') as f:
                    json.dump(data, f, indent=4)
                messagebox.showinfo("Saved", "Session saved successfully!")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to save session: {e}")

    def load_session(self):
        """Loads inputs from a JSON file."""
        file_path = filedialog.askopenfilename(filetypes=[("JSON files", "*.json")])
        if not file_path:
            return
            
        try:
            with open(file_path, 'r') as f:
                data = json.load(f)
            
            # Restore Inputs
            self.star_entry.delete(0, tk.END)
            self.star_entry.insert(0, data.get("star", ""))
            
            self.mode_var.set("Local Mode")
            
            # Restore Molecules first (needed for fetch)
            while self.molecule_entries:
                self.remove_last_molecule()
            for mol in data.get("molecules", []):
                self.add_molecule_inputs(
                    mol["name"], 
                    mol["b"], 
                    mol["N"], 
                    mol["v"],
                    mol.get("tie_b", "Custom"),
                    mol.get("tie_v", "Custom")
                )

            # Restore Files and Selection
            local_files = data.get("local_files", data.get("gm3_files", []))
            local_indices = data.get("local_selected_indices", data.get("gm3_selected_indices", []))
            
            if local_files:
                self.local_files = local_files
                self.local_selected_vars = [ctk.BooleanVar(value=False) for _ in local_files]
                
                for i in local_indices:
                    if 0 <= i < len(self.local_selected_vars):
                        self.local_selected_vars[i].set(True)
                        
                self.local_status_label.configure(text=f"{len(local_files)} files fetched ({len(local_indices)} selected)", text_color="green")
            
            if "local_target_wave" in data or "gm3_target_wave" in data:
                self.local_target_wave = data.get("local_target_wave", data.get("gm3_target_wave"))
            if "local_wrange" in data or "gm3_wrange" in data:
                self.local_wrange = data.get("local_wrange", data.get("gm3_wrange"))
            if "local_lambda_0_ref" in data or "gm3_lambda_0_ref" in data:
                if not hasattr(self, 'plot_data'):
                    self.plot_data = {}
                self.plot_data['lambda_0_ref'] = data.get("local_lambda_0_ref", data.get("gm3_lambda_0_ref"))
            
            self.wave_min.delete(0, tk.END)
            self.wave_min.insert(0, data.get("wave_min", ""))
            self.wave_max.delete(0, tk.END)
            self.wave_max.insert(0, data.get("wave_max", ""))
            
            self.abs_min.delete(0, tk.END)
            self.abs_min.insert(0, data.get("abs_min", ""))
            self.abs_max.delete(0, tk.END)
            self.abs_max.insert(0, data.get("abs_max", ""))

            self.n_knots_entry.delete(0, tk.END)
            self.n_knots_entry.insert(0, data.get("n_knots", "10"))
            self.cont_order_combo.set(data.get("cont_order", "Cubic"))
            self.knots_x_entry.delete(0, tk.END)
            self.knots_x_entry.insert(0, data.get("knots_x", ""))
            if data.get("show_knots", True):
                self.show_knots_switch.select()
            else:
                self.show_knots_switch.deselect()
            
            # Restore Plotting Style Settings
            self.theme_var.set(data.get("theme", "Dark"))
            self.grid_var.set(data.get("show_grid", True))
            if "remove_minor_grid" in data:
                self.remove_minor_grid_var.set(data["remove_minor_grid"])
            
            if "show_individual_lines" in data:
                self.show_individual_lines_var.set(data["show_individual_lines"])
                
            if "show_comp_lines" in data:
                self.show_comp_lines_var.set(data["show_comp_lines"])

            if "show_total_fit" in data:
                self.show_total_fit_switch.set(data["show_total_fit"])
            if "plot_style" in data:
                self.plot_style_combo.set(data["plot_style"])
            elif "step_plot" in data:
                self.plot_style_combo.set("Boxy" if data["step_plot"] else "Line")
            else:
                self.plot_style_combo.set("Line")
            
            self.title_entry.delete(0, tk.END)
            self.title_entry.insert(0, data.get("plot_title", ""))
            
            self.xlabel_entry.delete(0, tk.END)
            self.xlabel_entry.insert(0, data.get("xlabel", ""))
            
            self.ylabel_entry.delete(0, tk.END)
            self.ylabel_entry.insert(0, data.get("ylabel", ""))
            
            self.color_entry.set(data.get("label_color", "white"))
            self.cont_color_combo.set(data.get("cont_color", "red"))
            self.voigt_color_combo.set(data.get("voigt_color", "green"))

            if data.get("show_velocity_axis", False):
                self.show_velocity_axis_var.set(True)
            else:
                self.show_velocity_axis_var.set(False)

            if data.get("show_sigma", False):
                self.show_sigma_var.set(True)
            else:
                self.show_sigma_var.set(False)

            self.sigma_level_slider.set(data.get("sigma_level", 1.0))
            
            messagebox.showinfo("Loaded", "Session loaded successfully!")
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load session: {e}")

    def run_analysis(self):
        try:
            # Get Inputs
            star = self.star_entry.get()
            wave_range = [float(self.wave_min.get()), float(self.wave_max.get())]
            absorption_range = (float(self.abs_min.get()), float(self.abs_max.get()))
            
            molecules = []
            species_params = {}
            
            for i, entry in enumerate(self.molecule_entries):
                mol_name = entry['name'].get()
                if not mol_name: continue
                molecules.append(mol_name)
                
                # Check if tying is enabled for this species
                tie_b_enabled = (i > 0 and entry.get('b_tie') and entry['b_tie'].get() == "Tie to Species 1")
                tie_v_enabled = (i > 0 and entry.get('v_tie') and entry['v_tie'].get() == "Tie to Species 1")
                
                # Parse b values
                b_text = entry['b'].get().strip()
                if tie_b_enabled and not b_text:
                    b_list = species_params[0]['b']
                else:
                    b_list = [float(x.strip()) for x in b_text.split(',') if x.strip()]
                
                # Parse N values
                N_list = [float(x.strip()) for x in entry['N'].get().split(',') if x.strip()]
                
                # Parse v_rad values
                v_text = entry['v'].get().strip()
                if tie_v_enabled and not v_text:
                    v_list = species_params[0]['v_rad']
                else:
                    v_list = [float(x.strip()) for x in v_text.split(',') if x.strip()]
                
                params_dict = {'b': b_list, 'N': N_list, 'v_rad': v_list}
                
                if i > 0:
                    if tie_b_enabled:
                        params_dict['tie_b'] = True
                    if tie_v_enabled:
                        params_dict['tie_v_rad'] = True
                
                species_params[i] = params_dict
            
            if not molecules:
                messagebox.showwarning("Warning", "Please add at least one molecule.")
                return

            # Use co-added data from self.plot_data
            if not self.plot_data:
                messagebox.showwarning("Warning", "Please co-add spectra first.")
                return
            
            # Filter by wave_range
            full_wave = self.plot_data['wave']
            full_flux = self.plot_data['flux']
            
            mask = (full_wave >= wave_range[0]) & (full_wave <= wave_range[1])
            wave = full_wave[mask]
            flux = full_flux[mask]
            
            if len(wave) == 0:
                raise ValueError("No data points in the specified wavelength range.")

            n_knots = int(self.n_knots_entry.get())
            knots_x_str = self.knots_x_entry.get()
            knots_x_array = [float(x.strip()) for x in knots_x_str.split(',')] if knots_x_str.strip() else None
            
            species_file = 'species.txt'
            
            # Get Spline Order
            order_map = {"Linear": 1, "Quadratic": 2, "Cubic": 3}
            spline_order = order_map.get(self.cont_order_combo.get(), 3)

            # Read the "Fit σ" box. If the user typed garbage, just fall
            # back to 0.002 — same as the old hardcoded default.
            try:
                fit_std_dev = float(self.local_std_dev_entry.get().strip())
            except ValueError:
                fit_std_dev = 0.002

            fitresult, continuum_normalized_flux, continuum, std_dev = astrovoigtfit_run(
                wave, flux, molecules, species_params, absorption_range, n_knots, knots_x_array,
                species_file=species_file, spline_order=spline_order, std_dev=fit_std_dev,
            )
            
            lambda_0 = None
            if molecules:
                all_species = get_species_data('species.txt')
                if molecules[0] in all_species:
                    try:
                        lambda_0 = float(all_species[molecules[0]]['line'])
                    except (ValueError, TypeError):
                        pass

            individual_spectra = self.plot_data.get('individual_spectra')
            
            self.plot_data = {
                'wave': wave,
                'flux': flux,
                'norm_flux': continuum_normalized_flux,
                'continuum': continuum,
                'model': fitresult.best_fit,
                'residuals': flux - fitresult.best_fit,
                'mode': 'Local Mode',
                'fitresult': fitresult,
                'molecules': molecules,
                'lambda_0_ref': lambda_0
            }
            
            if individual_spectra:
                self.plot_data['individual_spectra'] = individual_spectra
            
            self.key_results_text = format_key_results(fitresult, molecules, star)
            self.update_plots()
            self.display_results(fitresult, molecules)
            self.tabview.set("Voigt Fit")
            
        except Exception as e:
            messagebox.showerror("Error", f"Analysis failed: {e}")
            import traceback
            traceback.print_exc()

    def update_plots(self, plot_title=None):
        plt.close('all') # Clear previous figures to prevent memory leak
        if not self.plot_data:
            return
            
        # Get Settings
        theme = self.theme_var.get()
        show_grid = self.grid_var.get()
        if not plot_title:
            plot_title = self.title_entry.get()
        xlabel = self.xlabel_entry.get() or "Wavelength (Å)"
        ylabel = self.ylabel_entry.get() or "Normalized Flux"
        label_color = self.color_entry.get()

        # Determine colors based on theme
        data_color = 'white' if theme == 'Dark' else 'black'
        cont_color = self.cont_color_combo.get() # Use selected continuum color
        
        # Plot Style (Line, Boxy, Scatter)
        plot_style = self.plot_style_combo.get()
        
        # Helper to plot data based on style
        def plot_data(ax, x, y, color, label=None, alpha=0.8):
            if plot_style == "Boxy":
                ax.step(x, y, where='mid', color=color, alpha=alpha, label=label)
            elif plot_style == "Scatter":
                ax.scatter(x, y, color=color, s=10, alpha=alpha, label=label, marker='o') # s is size
            else: # Line (Default)
                ax.plot(x, y, color=color, alpha=alpha, label=label)
        
        # Absorption Range Shading
        abs_min = float(self.abs_min.get())
        abs_max = float(self.abs_max.get())
        
        # Shading color: Light/Cyan for Dark mode, Gray for Light mode
        shade_color = 'cyan' if theme == 'Dark' else 'gray'
        shade_alpha = 0.2 if theme == 'Dark' else 0.3

        # 1. Continuum Fit Plot
        # For GM2, we now show the continuum fit in this tab too, as requested.
        # Logic is same for all modes now: Data + Continuum.
        
        if self.canvas1: 
            self.canvas1.get_tk_widget().destroy()
        if self.toolbar1:
            self.toolbar1.destroy()
            
        # Create 2 subplots: Top (Data + Continuum), Bottom (Normalized)
        # Equal height ratios [1, 1]
        # Create 2 subplots: Top (Data + Continuum), Bottom (Normalized)
        # Equal height ratios [1, 1]
        fig1, (ax1_top, ax1_bot) = plt.subplots(2, 1, figsize=(6, 6), dpi=100, gridspec_kw={'height_ratios': [1, 1]}, sharex=True)
        
        remove_minor = self.remove_minor_grid_var.get() if hasattr(self, 'remove_minor_grid_var') else False
        style_axis(ax1_top, theme, show_grid, label_color, remove_minor_grid=remove_minor)
        style_axis(ax1_bot, theme, show_grid, label_color, remove_minor_grid=remove_minor)
        
        # Prepare Data Arrays for Continuum Plot (Handle Velocity Conversion if enabled)
        plot_wave_cont = self.plot_data['wave']
        xlabel_cont = self.xlabel_entry.get() or "Wavelength (Å)"
        
        # Speed of light in km/s
        c_light = 299792.458 

        if self.show_velocity_axis_var.get() and self.plot_data.get('lambda_0_ref'):
            lambda_0_ref = self.plot_data['lambda_0_ref']
            # velocity = c * (wave - lambda_0) / lambda_0
            plot_wave_cont = c_light * (plot_wave_cont - lambda_0_ref) / lambda_0_ref
            xlabel_cont = "Velocity (km/s)" # Override label
            
            # Convert absorption range to velocity for shading
            abs_min_plot = c_light * (abs_min - lambda_0_ref) / lambda_0_ref
            abs_max_plot = c_light * (abs_max - lambda_0_ref) / lambda_0_ref
        else:
            abs_min_plot = abs_min
            abs_max_plot = abs_max

        # Top Plot: Flux + Continuum
        label_text = 'Data (Co-added)'
        plot_data(ax1_top, plot_wave_cont, self.plot_data['flux'], data_color, label=label_text)
            
        # GM3: Plot Individual Spectra
        if 'individual_spectra' in self.plot_data:
            alpha_val = self.spectra_alpha_slider.get()
            for i, (w, f) in enumerate(self.plot_data['individual_spectra']):
                # Convert individual spectra to velocity if enabled
                if self.show_velocity_axis_var.get() and self.plot_data.get('lambda_0_ref'):
                    w_plot = c_light * (w - lambda_0_ref) / lambda_0_ref
                else:
                    w_plot = w
                ax1_top.plot(w_plot, f, alpha=alpha_val, linewidth=0.8, label=f'Spec {i+1}' if i < 3 else None) # Label only first few to avoid clutter
            
        ax1_top.plot(plot_wave_cont, self.plot_data['continuum'], color=cont_color, linewidth=2, label='Continuum')
        ax1_top.axvspan(abs_min_plot, abs_max_plot, color=shade_color, alpha=shade_alpha, label='Absorption Range')
        
        # Plot Knots in Continuum Panel (GM2/GM3 only)
        if self.show_knots_switch.get():
            fitresult = self.plot_data.get('fitresult')
            if fitresult:
                used_knots_x = fitresult.userkws.get('knot_x_array')
                if used_knots_x is not None:
                    # Convert knots to velocity if enabled
                    if self.show_velocity_axis_var.get() and self.plot_data.get('lambda_0_ref'):
                        used_knots_x_plot = c_light * (used_knots_x - lambda_0_ref) / lambda_0_ref
                    else:
                        used_knots_x_plot = used_knots_x

                    knot_y_values = []
                    i = 0
                    while True:
                        name = f'knot_y_{i}'
                        if name in fitresult.params:
                            knot_y_values.append(fitresult.params[name].value)
                            i += 1
                        else:
                            break
                    if len(knot_y_values) == len(used_knots_x_plot):
                        ax1_top.plot(used_knots_x_plot, knot_y_values, 'o', color='orange', markersize=6, label='Knots')

        ax1_top.set_title(plot_title if plot_title else "Continuum Fit")
        ax1_top.set_ylabel("Flux")
        ax1_top.legend()
        
        # Bottom Plot: Normalized Flux
        plot_data(ax1_bot, plot_wave_cont, self.plot_data['norm_flux'], data_color)
            
        ax1_bot.axhline(1.0, color=cont_color, linestyle='--', linewidth=1.5)
        ax1_bot.axvspan(abs_min_plot, abs_max_plot, color=shade_color, alpha=shade_alpha)
        
        ax1_bot.set_xlabel(xlabel_cont)
        ax1_bot.set_ylabel("Norm. Flux")
        
        if theme == 'Dark':
            fig1.patch.set_facecolor('#2b2b2b')
            ax1_top.set_facecolor('#2b2b2b')
            ax1_bot.set_facecolor('#2b2b2b')
            
        plt.subplots_adjust(hspace=0.05)
        
        self.canvas1 = FigureCanvasTkAgg(fig1, master=self.tabview.tab("Continuum Fit"))
        self.canvas1.draw()
        self.canvas1.get_tk_widget().pack(fill="both", expand=True)
        
        self.toolbar1 = NavigationToolbar2Tk(self.canvas1, self.tabview.tab("Continuum Fit"))
        self.toolbar1.update()
        self.toolbar1.pack(side="bottom", fill="x")

        # 2. Voigt Fit Plot
        if self.canvas2:
            self.canvas2.get_tk_widget().destroy()
        if self.toolbar2:
            self.toolbar2.destroy()
            
        fig2, (ax2, ax3) = plt.subplots(2, 1, figsize=(6, 6), dpi=100, gridspec_kw={'height_ratios': [3, 1]}, sharex=True)
        
        remove_minor = self.remove_minor_grid_var.get() if hasattr(self, 'remove_minor_grid_var') else False
        style_axis(ax2, theme, show_grid, label_color, remove_minor_grid=remove_minor)
        style_axis(ax3, theme, show_grid, label_color, remove_minor_grid=remove_minor)
        
        fit_color = self.voigt_color_combo.get() # Use selected voigt color
        
        xlabel = self.xlabel_entry.get() or "Wavelength (Å)"
        plot_title = self.title_entry.get()
        
        # Prepare Data Arrays (Handle Velocity Conversion if enabled)
        plot_wave = self.plot_data['wave']
        
        if self.show_velocity_axis_var.get() and self.plot_data.get('lambda_0_ref'):
            lambda_0_ref = self.plot_data['lambda_0_ref']
            # velocity = c * (wave - lambda_0) / lambda_0
            plot_wave = c_light * (plot_wave - lambda_0_ref) / lambda_0_ref
            xlabel = "Velocity (km/s)" # Override label
            
        fitresult = self.plot_data.get('fitresult')
        
        # GM2: Check Normalize Toggle
        if self.normalize_plot_switch.get():
            # Plot Re-normalized (Flux / Continuum)
            # Norm Flux is already in plot_data['norm_flux']
            # Norm Model = Total Model / Continuum. 
            # Note: Total Model = Continuum * Absorption. So Norm Model = Absorption.
            # We can calculate it on the fly.
            
            norm_model = self.plot_data['model'] / self.plot_data['continuum']
            
            # 1. Normalized Flux
            plot_data(ax2, plot_wave, self.plot_data['norm_flux'], data_color, label='Data (Norm)', alpha=0.7)
                
            # 2. Normalized Model (Absorption only)
            if self.show_total_fit_switch.get():
            	ax2.plot(plot_wave, norm_model, color=fit_color, linewidth=2, label='Fit (Norm)')
            
            # 3. Flat Continuum Line (if enabled)
            if self.show_continuum_switch.get():
                cont_color = self.cont_color_combo.get()
                ax2.axhline(1.0, color=cont_color, linestyle='--', linewidth=1.5, label='Continuum')

            # 4. Individual Lines (if enabled)
            if self.show_individual_lines_var.get():
                fitresult = self.plot_data.get('fitresult')
                if fitresult:
                    # Extract params
                    n_species = int(fitresult.params['n_species'].value)
                    v_resolution = fitresult.params['v_resolution'].value
                    
                    # Loop through species
                    molecules = self.plot_data.get('keys', []) # We need molecule names. Using keys from fitresult or stored plot_data
                    # Reconstruct molecules list if not stored directly. 
                    # Actually 'molecules' list is local to run_analysis. We should store it in plot_data.
                    # For now, let's try to get it from keys/params or self.mol_scroll checkboxes if they match.
                    # Better approach: Extract from params by prefix logic or assume generic names if missing.
                    
                    # Let's extract per-species params from fitresult for master_function
                    
                    # We need to build kwargs_one for each species.
                    # The fitresult params have flattened N, b, v_rad, etc.
                    
                    for i in range(n_species):
                        # Suffix i corresponds to species index
                        suffix = str(i)
                        
                        # Build kwargs for THIS species only
                        kwargs_species = {}
                        
                        # Get transition params (fixed)
                        n_trans = int(fitresult.params[f'n_trans_{i}'].value)
                        lambdas = []
                        fs = []
                        gammas = []
                        for t in range(n_trans):
                            lambdas.append(fitresult.params[f'lambda_{i}_{t}'].value)
                            fs.append(fitresult.params[f'f_{i}_{t}'].value)
                            gammas.append(fitresult.params[f'gamma_{i}_{t}'].value)
                        
                        kwargs_species[f'lambda_{i}'] = np.array(lambdas)
                        kwargs_species[f'f_{i}'] = np.array(fs)
                        kwargs_species[f'gamma_{i}'] = np.array(gammas)
                        
                        # Get component params (fitted)
                        n_comp = int(fitresult.params[f'n_component_{i}'].value)
                        Ns = []
                        bs = []
                        vs = []
                        for c in range(n_comp):
                            Ns.append(fitresult.params[f'N_{i}_{c}'].value)
                            bs.append(fitresult.params[f'b_{i}_{c}'].value)
                            vs.append(fitresult.params[f'v_rad_{i}_{c}'].value)
                            
                        kwargs_species[f'N_{i}'] = np.array(Ns)
                        kwargs_species[f'b_{i}'] = np.array(bs)
                        kwargs_species[f'v_rad_{i}'] = np.array(vs)
                        
                        # Generate Model for this species
                        # We need to call master_function with ONLY this species params.
                        # BUT master_function expects suffixes to match. 
                        # Since we only pass one group with suffix 'i', it should work.
                        
                        with np.errstate(divide='ignore', invalid='ignore'):
                            model_species = master_function(
                                self.plot_data['wave'], # Always pass original wave to master_function
                                v_resolution=v_resolution,
                                **kwargs_species
                            )
                        
                        # Handle NaNs (master_function returns transmission, so 1.0 is continuum)
                        model_species = np.nan_to_num(model_species, nan=1.0)
                        
                        # Plot
                        # Try to get species name from stored list if possible, else generic
                        label_name = f'Species {i+1}'
                        # If we stored molecules in plot_data, use it.
                        if 'molecules' in self.plot_data and i < len(self.plot_data['molecules']):
                            label_name = self.plot_data['molecules'][i]
                        elif hasattr(self, 'molecule_vars') and i < len(self.molecule_vars): # Fallback to checkbox labels? risky if order changed.
                            pass
                            
                        ax2.plot(plot_wave, model_species, linestyle='--', linewidth=1.2, alpha=0.9, label=label_name)

            # 6. Component-wise Lines (if enabled)
            if self.show_comp_lines_var.get():
                fitresult = self.plot_data.get('fitresult')
                if fitresult:
                    n_species = int(fitresult.params['n_species'].value)
                    v_resolution = fitresult.params['v_resolution'].value
                    
                    for i in range(n_species):
                        # Extract basic species params (lambda, f, gamma)
                        n_trans = int(fitresult.params[f'n_trans_{i}'].value)
                        lambdas = [fitresult.params[f'lambda_{i}_{t}'].value for t in range(n_trans)]
                        fs = [fitresult.params[f'f_{i}_{t}'].value for t in range(n_trans)]
                        gammas = [fitresult.params[f'gamma_{i}_{t}'].value for t in range(n_trans)]
                        
                        n_comp = int(fitresult.params[f'n_component_{i}'].value)
                        
                        # Loop over components
                        for c in range(n_comp):
                            kwargs_comp = {}
                            # Transitions (apply to all components)
                            kwargs_comp[f'lambda_{i}'] = np.array(lambdas)
                            kwargs_comp[f'f_{i}'] = np.array(fs)
                            kwargs_comp[f'gamma_{i}'] = np.array(gammas)
                            
                            # Single component params (as single-element arrays)
                            kwargs_comp[f'N_{i}'] = np.array([fitresult.params[f'N_{i}_{c}'].value])
                            kwargs_comp[f'b_{i}'] = np.array([fitresult.params[f'b_{i}_{c}'].value])
                            kwargs_comp[f'v_rad_{i}'] = np.array([fitresult.params[f'v_rad_{i}_{c}'].value])
                            
                            with np.errstate(divide='ignore', invalid='ignore'):
                                model_comp = master_function(
                                    self.plot_data['wave'],
                                    v_resolution=v_resolution,
                                    **kwargs_comp
                                )
                            model_comp = np.nan_to_num(model_comp, nan=1.0)
                            
                            # Label
                            label_name = f'Species {i+1} (Cloud {c+1})'
                            if 'molecules' in self.plot_data and i < len(self.plot_data['molecules']):
                                label_name = f"{self.plot_data['molecules'][i]} (Cloud {c+1})"
                            
                            # Use dotted line for component-wise
                            ax2.plot(plot_wave, model_comp, linestyle=':', linewidth=1.5, alpha=0.9, label=label_name)

            ax2.set_ylabel("Norm. Flux")
            
            # Update residuals for normalized view
            # Residuals = Norm Flux - Norm Model
            norm_model = self.plot_data['model'] / self.plot_data['continuum']
            residuals = self.plot_data['norm_flux'] - norm_model
            
        else:
            # Plot Raw (Curved)
            # 1. Raw Flux
            plot_data(ax2, plot_wave, self.plot_data['flux'], data_color, label='Data', alpha=0.7)

            # 2. Total Model (Best Fit)
            if self.show_total_fit_switch.get():
                ax2.plot(plot_wave, self.plot_data['model'], color=fit_color, linewidth=2, label='Total Fit')
            
            # 3. Continuum (if enabled)
            if self.show_continuum_switch.get():
                cont_color = self.cont_color_combo.get()
                ax2.plot(plot_wave, self.plot_data['continuum'], color=cont_color, linestyle='--', linewidth=1.5, label='Continuum')
            
            # 4. Knots (if enabled)
            if self.show_knots_switch.get():
                fitresult = self.plot_data.get('fitresult')
                if fitresult:
                    # Extract knots
                    used_knots_x = fitresult.userkws.get('knot_x_array')
                    knot_y_values = []
                    i = 0
                    while True:
                        name = f'knot_y_{i}'
                        if name in fitresult.params:
                            knot_y_values.append(fitresult.params[name].value)
                            i += 1
                        else:
                            break
                    
                    if used_knots_x is not None and len(knot_y_values) == len(used_knots_x):
                         # Knots X must ALSO be converted if Velocity Axis is ON
                        plot_knots_x2 = used_knots_x
                        if self.show_velocity_axis_var.get() and self.plot_data.get('lambda_0_ref'):
                            lambda_0 = self.plot_data['lambda_0_ref']
                            plot_knots_x2 = [c_light * (kx - lambda_0) / lambda_0 for kx in used_knots_x]
                            
                        ax2.plot(plot_knots_x2, knot_y_values, 'o', color='orange', markersize=6, label='Knots')

            # 5. Individual Lines (if enabled)
            if self.show_individual_lines_var.get():
                fitresult = self.plot_data.get('fitresult')
                if fitresult:
                    n_species = int(fitresult.params['n_species'].value)
                    v_resolution = fitresult.params['v_resolution'].value
                    
                    for i in range(n_species):
                        suffix = str(i)
                        kwargs_species = {}
                        
                        n_trans = int(fitresult.params[f'n_trans_{i}'].value)
                        lambdas = []
                        fs = []
                        gammas = []
                        for t in range(n_trans):
                            lambdas.append(fitresult.params[f'lambda_{i}_{t}'].value)
                            fs.append(fitresult.params[f'f_{i}_{t}'].value)
                            gammas.append(fitresult.params[f'gamma_{i}_{t}'].value)
                        
                        kwargs_species[f'lambda_{i}'] = np.array(lambdas)
                        kwargs_species[f'f_{i}'] = np.array(fs)
                        kwargs_species[f'gamma_{i}'] = np.array(gammas)
                        
                        n_comp = int(fitresult.params[f'n_component_{i}'].value)
                        Ns = []
                        bs = []
                        vs = []
                        for c in range(n_comp):
                            Ns.append(fitresult.params[f'N_{i}_{c}'].value)
                            bs.append(fitresult.params[f'b_{i}_{c}'].value)
                            vs.append(fitresult.params[f'v_rad_{i}_{c}'].value)
                            
                        kwargs_species[f'N_{i}'] = np.array(Ns)
                        kwargs_species[f'b_{i}'] = np.array(bs)
                        kwargs_species[f'v_rad_{i}'] = np.array(vs)
                        
                        with np.errstate(divide='ignore', invalid='ignore'):
                            model_species = master_function(
                                self.plot_data['wave'],
                                v_resolution=v_resolution,
                                **kwargs_species
                            )
                            
                        model_species = np.nan_to_num(model_species, nan=1.0)
                        
                        # Scale by continuum for Raw Flux view!
                        model_species_flux = model_species * self.plot_data['continuum']
                        
                        label_name = f'Species {i+1}'
                        if 'molecules' in self.plot_data and i < len(self.plot_data['molecules']):
                            label_name = self.plot_data['molecules'][i]
                        elif hasattr(self, 'molecule_vars') and i < len(self.molecule_vars):
                            pass
                            
                        ax2.plot(plot_wave, model_species_flux, linestyle='--', linewidth=1.2, alpha=0.9, label=label_name)

            # 6. Component-wise Lines (if enabled)
            if self.show_comp_lines_var.get():
                fitresult = self.plot_data.get('fitresult')
                if fitresult:
                    n_species = int(fitresult.params['n_species'].value)
                    v_resolution = fitresult.params['v_resolution'].value
                    
                    for i in range(n_species):
                        n_trans = int(fitresult.params[f'n_trans_{i}'].value)
                        lambdas = [fitresult.params[f'lambda_{i}_{t}'].value for t in range(n_trans)]
                        fs = [fitresult.params[f'f_{i}_{t}'].value for t in range(n_trans)]
                        gammas = [fitresult.params[f'gamma_{i}_{t}'].value for t in range(n_trans)]
                        
                        n_comp = int(fitresult.params[f'n_component_{i}'].value)
                        
                        for c in range(n_comp):
                            kwargs_comp = {}
                            kwargs_comp[f'lambda_{i}'] = np.array(lambdas)
                            kwargs_comp[f'f_{i}'] = np.array(fs)
                            kwargs_comp[f'gamma_{i}'] = np.array(gammas)
                            
                            kwargs_comp[f'N_{i}'] = np.array([fitresult.params[f'N_{i}_{c}'].value])
                            kwargs_comp[f'b_{i}'] = np.array([fitresult.params[f'b_{i}_{c}'].value])
                            kwargs_comp[f'v_rad_{i}'] = np.array([fitresult.params[f'v_rad_{i}_{c}'].value])
                            
                            with np.errstate(divide='ignore', invalid='ignore'):
                                model_comp = master_function(
                                    self.plot_data['wave'],
                                    v_resolution=v_resolution,
                                    **kwargs_comp
                                )
                            model_comp = np.nan_to_num(model_comp, nan=1.0)
                            
                            # Scale by continuum for Raw Flux view!
                            model_comp_flux = model_comp * self.plot_data['continuum']

                            label_name = f'Species {i+1} (Cloud {c+1})'
                            if 'molecules' in self.plot_data and i < len(self.plot_data['molecules']):
                                label_name = f"{self.plot_data['molecules'][i]} (Cloud {c+1})"
                            
                            ax2.plot(plot_wave, model_comp_flux, linestyle=':', linewidth=1.5, alpha=0.9, label=label_name)

            ax2.set_ylabel("Flux")
            residuals = self.plot_data['residuals'] # Already calculated as Raw - Model
        
        ax2.legend()
        ax2.set_title(plot_title if plot_title else "Voigt Profile Fit")
        
        
        # Residuals
        plot_data(ax3, plot_wave, residuals, data_color, alpha=0.7)
            
        ax3.axhline(0, color=cont_color, linestyle='--', alpha=0.5)
        
        # Show Sigma Shading (Residuals)
        if hasattr(self, 'show_sigma_var') and self.show_sigma_var.get():
            # Calculate sigma with ddof=p (number of fitted parameters)
            p = 0
            fitresult = self.plot_data.get('fitresult')
            if fitresult and hasattr(fitresult, 'nvarys'):
                 p = fitresult.nvarys
            
            sigma_val = np.std(residuals, ddof=p)
            
            # Get Multiplier
            multiplier = self.sigma_level_slider.get()
            boundary = sigma_val * multiplier
            
            # Shade +sigma and -sigma
            ax3.fill_between(plot_wave, -boundary, boundary, color='gray', alpha=0.3, label=f'±{multiplier:.1f}$\sigma$ ({boundary:.3f})')
            ax3.legend(fontsize='small', loc='upper right')

        ax3.set_xlabel(xlabel)
        ax3.set_ylabel("Residuals")
        
        if theme == 'Dark':
            fig2.patch.set_facecolor('#2b2b2b')
            ax2.set_facecolor('#2b2b2b')
            ax3.set_facecolor('#2b2b2b')
            
        plt.subplots_adjust(hspace=0.05)
        
        self.canvas2 = FigureCanvasTkAgg(fig2, master=self.tabview.tab("Voigt Fit"))
        self.canvas2.draw()
        self.canvas2.get_tk_widget().pack(fill="both", expand=True)
        
        self.toolbar2 = NavigationToolbar2Tk(self.canvas2, self.tabview.tab("Voigt Fit"))
        self.toolbar2.update()
        self.toolbar2.pack(side="bottom", fill="x")

class EdiblesApp(ctk.CTk):
    def __init__(self):
        super().__init__()

        # Window configuration
        self.title("EDIBLES GUI - Astro Voigt Fit")
        self.geometry("1200x900")
        
        # Configure grid
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)

        # Tab Counter for unique naming
        self.tab_counter = 0

        # Main Tabview for Multiple Analysis Tabs
        self.main_tabs = ctk.CTkTabview(self)
        self.main_tabs.grid(row=0, column=0, padx=5, pady=5, sticky="nsew")
        
        # Add initial tab
        self.add_new_tab()
        
        # Add Tab Button (Overlay or separate frame)
        # Since CTkTabview doesn't have a built-in "+" button, we'll add a button below/above
        # or we can use a trick: add a tab named "+" and catch the click event.
        # For simplicity, let's add a button in a top bar.
        
        # Re-layout to make room for top bar
        self.main_tabs.grid(row=1, column=0, padx=5, pady=5, sticky="nsew")
        self.grid_rowconfigure(1, weight=1)
        self.grid_rowconfigure(0, weight=0)
        
        self.top_bar = ctk.CTkFrame(self, height=40, fg_color="transparent")
        self.top_bar.grid(row=0, column=0, sticky="ew", padx=5, pady=5)
        
        ctk.CTkButton(self.top_bar, text="+ New Tab", width=100, command=self.add_new_tab).pack(side="left", padx=5)
        ctk.CTkButton(self.top_bar, text="× Close Tab", width=100, fg_color="#D32F2F", hover_color="#B71C1C", command=self.close_current_tab).pack(side="left", padx=5)

        # Handle Window Close Event
        self.protocol("WM_DELETE_WINDOW", self.on_closing)

    def add_new_tab(self):
        self.tab_counter += 1
        tab_name = f"Analysis {self.tab_counter}"
        
        # Ensure name is unique (just in case)
        while tab_name in self.main_tabs._tab_dict:
            self.tab_counter += 1
            tab_name = f"Analysis {self.tab_counter}"
            
        self.main_tabs.add(tab_name)
        
        # Create AnalysisTab instance inside the new tab
        analysis_frame = AnalysisTab(self.main_tabs.tab(tab_name))
        analysis_frame.pack(fill="both", expand=True)
        
        self.main_tabs.set(tab_name)

    def close_current_tab(self):
        # Close the CURRENTLY SELECTED tab
        current_tab = self.main_tabs.get()
        if not current_tab:
            return

        # Don't close the last tab
        if len(self.main_tabs._tab_dict) <= 1:
            messagebox.showinfo("Info", "Cannot close the last tab.")
            return

        # Smart Focus Logic:
        # If we close a tab, we want to focus the one to the LEFT (previous index).
        # Unless we are closing the first tab, then we focus the new first tab (next index).
        
        tab_names = list(self.main_tabs._tab_dict.keys())
        try:
            index = tab_names.index(current_tab)
            
            # Determine which tab to select after closing
            if index > 0:
                next_tab_name = tab_names[index - 1]
            else:
                # We are closing the first tab (index 0), so select the one that was at index 1
                next_tab_name = tab_names[index + 1]
                
            self.main_tabs.delete(current_tab)
            self.main_tabs.set(next_tab_name)
            
        except ValueError:
            # Fallback if something goes wrong with index finding
            self.main_tabs.delete(current_tab)

    def on_closing(self):
        """Confirm before closing the application."""
        if messagebox.askokcancel("Quit", "Are you sure you want to close all windows?\nMake sure to save all analysis."):
            self.destroy()
            import sys
            sys.exit(0)

if __name__ == "__main__":
    app = EdiblesApp()
    app.mainloop()
