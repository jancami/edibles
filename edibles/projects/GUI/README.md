# EDIBLES GUI

A small Tkinter app for co-adding EDIBLES spectra and fitting Voigt profiles to interstellar absorption lines. This README is for whoever picks up the code next.

## What it does, end to end

You give it a star name and a species (say Li 6708). The GUI talks to EDIBLES, lists the observation files that cover that wavelength, lets you tick the ones you want, and co-adds them into a single cleaner spectrum. Then you set initial guesses for column density, b, and v_rad for each velocity component, and click Run Analysis. Behind the scenes it fits a spline continuum and a Voigt absorption model at the same time (not sequentially), and shows you the fit and residuals in two plot tabs.

## Run it

```bash
python local_gui.py
```

You need `customtkinter`, `matplotlib`, `numpy`, `scipy`, `lmfit`, and the `edibles` package.

## The files

The GUI itself is in one file, but the science is spread across several. The split is deliberate: if you only need to tweak the look, you stay inside `local_gui.py`. If you need to change how the fit works, you go elsewhere.

**`local_gui.py`** is the whole GUI. Buttons, plots, layout, session save and load. No science in here.

**`main_run.py`** is the thin glue between the GUI and the fitter. Reads `species.txt`, packages parameters, calls the fit, returns the continuum and normalized flux.

**`astrovoigtfit.py`** is the fitting machinery. Builds the lmfit Parameters object, places spline knots, runs the joint fit.

**`main.py`** is the Voigt absorption model. Given physical parameters, it returns a transmission spectrum.

**`flux_wave_find.py`** pulls one spectrum from the EDIBLES archive (wavelength cropped, median normalized).

**`co_adding_flux.py`** co adds a list of spectra into one. Plain numpy, no I/O.

**`species.txt`** has the atomic data per species (rest wavelengths, oscillator strengths, gammas, fit ranges).

The low-level Voigt math (the line profile, optical depth, FWHM helpers) comes from `edibles.utils.voigt_profile`, not from this project.

## How the pieces talk to each other

```
local_gui.py
    co_adding_flux.coadd_spectra
    flux_wave_find.wave_flux_data  ->  edibles.EdiblesSpectrum
    edibles.EdiblesOracle (direct, for the file list)
    main.master_function (for the per species overlay curves)
    main_run.get_species_data (reads species.txt)
    main_run.astrovoigtfit_run
        astrovoigtfit.astro_simultaneous_fit
            _generate_smart_knots
            continuum_voigt_wrapper      (the model lmfit fits)
                Voigt_fit_wrapper        (calls main.master_function)
                spline(knot_x, knot_y)
```

## What each button does

**Fetch Files.** Looks up every EDIBLES file for that star covering the target wavelength. Goes through `EdiblesOracle.getFilteredObsList`. Nothing is loaded yet, just listed.

**Select Files.** Opens a popup with checkboxes so you can pick which files to actually co add.

**Co add Selected.** For each ticked file, calls `flux_wave_find.wave_flux_data` to get the (wave, flux) chunk, then hands the lot to `co_adding_flux.coadd_spectra`. Stores the result in `self.plot_data` and plots it.

**Run Analysis.** Reads the species inputs, then calls `main_run.astrovoigtfit_run` which delegates to `astrovoigtfit.astro_simultaneous_fit`. The fit returns the joint best fit (continuum * absorption). Back in `main_run`, the continuum gets reconstructed from the fitted spline knot values so the GUI can plot it separately and compute a normalized flux.

## The Sigma and Fit σ fields

They look similar but feed different things.

* **Sigma** is the per-spectrum noise level you tell the co-add to use. It weights each input spectrum before averaging.
* **Fit σ** is the noise of the *co-added* spectrum. It goes into the fit as `std_dev`, which lmfit turns into `weights = 1 / std_dev`. If your co-added spectrum is noisier than 0.002, bump this up.

## Where to change things

If you want to change the layout, colors, fonts, or widgets, open `local_gui.py`.

If you want to change which plotting options are available, that's also `local_gui.py` (the Plot Settings tab is built in `setup_settings_tab`).

If you want to change how the fit is set up (bounds, ties, weights), go to `astrovoigtfit.py`.

If you want to change the Voigt model itself, that's `main.py`, and the underlying math is in `edibles.utils.voigt_profile`.

If you want to change how files are loaded from EDIBLES, that's `flux_wave_find.py`.

If you want to change how spectra are co added, that's `co_adding_flux.py`.

If you want to add or edit a species or its atomic data, that's `species.txt`.

## The Mode dropdown

Right now only **Local Mode** does anything. **Global Mode** is in the dropdown as a placeholder so the next person has an obvious place to plug it in. If you're that person:

1. Add whatever extra widgets Global Mode needs in `setup_inputs_tab` and `setup_file_selection`.
2. In `coadd_spectra_action` and `run_analysis`, branch on `self.mode_var.get()`.
3. Tag `self.plot_data['mode']` so `update_plots` can do mode specific things if you need.
4. Take out the warning popup in `on_mode_change`.

## Sessions

The 💾 Save and 📂 Load buttons just dump everything in the input forms to a JSON file and reload it later. They don't save the actual spectra or the fit, only the inputs and the selected file list. Useful for getting back to the same setup after restarting.

## A few gotchas

* `local_gui.py` imports `EdiblesOracle` directly in `fetch_files_local`. The science isn't really happening there, just a catalog lookup, so I left it as is. If you want a cleaner separation you can move that into `flux_wave_find.py`.
* `main_run.astrovoigtfit_run` returns four things: `fitresult`, `normalized_flux`, `continuum`, and `residual_std`. The last one is the std of `(flux - best_fit)`, which is NOT the same as the `std_dev` argument you pass in. The name shadowing is unfortunate but it was easier to keep the return tuple stable.
* `other_functions.py` is no longer imported by anything (the Voigt math comes from the `edibles` package now). It's safe to delete if you want a cleaner tree.
* The `mode_var` flag is still stamped into `self.plot_data['mode']` in a couple of places even though there's only one mode. That's intentional so that adding Global Mode later doesn't require touching the plotting code in the obvious spots.

## Where the data flows through `self.plot_data`

`plot_data` is the dict the GUI passes around between the co-add step, the fit step, and the plotting code. The keys it cares about:

* `wave`, `flux` (the co-added spectrum)
* `norm_flux`, `continuum` (filled in after the fit)
* `model`, `residuals` (best fit and the leftovers)
* `individual_spectra` (the list of (wave, flux) tuples that went into the co-add, used to draw the thin lines under the main curve)
* `lambda_0_ref` (rest wavelength of the first species, used for the velocity axis option)
* `fitresult` (the full lmfit ModelResult, used to extract per species and per component overlay curves)
* `mode` (currently always `'Local Mode'`)

If you add new plots, this is where the data is.

## Adding a new species

Open `species.txt` and copy the format of one of the existing rows. You need the species name, the line wavelength, the lambda / f / gamma triplets in brackets (one entry per transition), and the wavelength range you want the GUI to default to. The GUI picks it up automatically the next time you start it.
