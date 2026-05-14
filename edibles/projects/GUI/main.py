"""The Voigt absorption model — this is where the physics lives.

Two functions:
    mother_function   does the real work: sums up Voigt optical depths over
                      all transitions and velocity components for a single
                      species, then convolves with the instrumental resolution.
    master_function   is a thin wrapper around mother_function. It takes
                      per-species kwargs (lambda_1st, N_2nd, ...) — which is
                      the format astrovoigtfit.py's lmfit wrappers produce —
                      flattens them, and hands one big batch to mother_function.

The low-level Voigt math (FWHM, line profile, optical depth) is imported from
edibles.utils.voigt_profile; we just orchestrate it here.
"""

import numpy as np
from scipy.interpolate import interp1d
import astropy.constants as cst
from scipy.special import wofz
from edibles.utils.voigt_profile import voigt_profile, fwhm2sigma, VoigtFWHM, getVGrid
from scipy.ndimage import gaussian_filter
from edibles.utils.voigt_profile import voigt_optical_depth

# mother_function is used to model the spectrum
def mother_function(wavegrid, lambda0=0.0, f=0.0, gamma=0.0, b=0.0,
                    N=0.0, v_rad=0.0, v_resolution=0.0, n_step=25):
    
    # Pre-compute all Voigt FWHMs at once
    Voigt_FWHM = VoigtFWHM(lambda0, gamma, b)
    FWHM2use = np.min(np.append(Voigt_FWHM, v_resolution))
    
    # Optimize grid creation
    xgrid_test = np.asarray(wavegrid)
    dv_xgrid = np.median(np.diff(xgrid_test)) / np.mean(xgrid_test) * cst.c.to("km/s").value
    n_step = max(7, np.ceil(FWHM2use / dv_xgrid))
    
    # Vectorized wavelength range calculation
    v_stepsize = FWHM2use / n_step
    bluewaves = lambda0 * (1.0 + (v_rad - 8.5 * Voigt_FWHM) / cst.c.to("km/s").value)
    redwaves = lambda0 * (1.0 + (v_rad + 8.5 * Voigt_FWHM) / cst.c.to("km/s").value)
    
    minwave = max(np.min(bluewaves), wavegrid.min())
    maxwave = min(np.max(redwaves), wavegrid.max())
    
    # Check if the absorption line is outside the observed wavelength range
    if maxwave <= minwave:
        # No overlap - return flat spectrum (no absorption)
        return np.ones_like(wavegrid)
    
    # Create reference grid
    n_v = int(np.ceil((maxwave - minwave) / minwave * cst.c.to("km/s").value / v_stepsize))
    refgrid = minwave * (1.0 + np.arange(n_v) * v_stepsize / cst.c.to("km/s").value)
    
    # Process all lines at once if possible, or use parallel processing
    allcomponents = np.zeros((n_v, len(lambda0)))
    
    for lineloop in range(len(lambda0)):
        dv = getVGrid(
            lambda0[lineloop],
            gamma[lineloop],
            b[lineloop],
            v_resolution,
            n_step)
        thiswavegrid = lambda0[lineloop] * (1.0 + dv / cst.c.to("km/s").value)
        tau = voigt_optical_depth(
            thiswavegrid,
            lambda0=lambda0[lineloop],
            b=b[lineloop],
            N=N[lineloop],
            f=f[lineloop],
            gamma=gamma[lineloop],
            v_rad=v_rad[lineloop],
        )
        # if debug:
        #     print("Max tau:", tau.max())
        # Shift to the proper wavelength given the radial velocity
        vel = dv + v_rad[lineloop]
        thiswavegrid = lambda0[lineloop] * (
                1.0 + vel / cst.c.to("km/s").value
        )
        # Interpolate to reference grid
        interpolationfunction = interp1d(
            thiswavegrid, tau, kind="cubic", fill_value="extrapolate"
        )
        tau_grid = interpolationfunction(refgrid)
        tau_grid[np.where(refgrid > np.max(thiswavegrid))] = 0
        tau_grid[np.where(refgrid < np.min(thiswavegrid))] = 0
        # plt.plot(thiswavegrid,tau,marker="+")
        # plt.plot(refgrid,tau_grid, color='red')
        # plt.show()

        allcomponents[:, lineloop] = tau_grid
    
    # Vectorized operations for combining components
    tau = np.sum(allcomponents, axis=1)
    AbsorptionLine = np.exp(-tau)
    
    # Optimized smoothing
    smooth_sigma = fwhm2sigma(v_resolution) / v_stepsize
    gauss_smooth = gaussian_filter(AbsorptionLine, sigma=smooth_sigma)
    
    # Use faster interpolation for final step
    interpolated_model = np.interp(wavegrid, refgrid, gauss_smooth, left=1, right=1)
    
    return interpolated_model



# wrapper function of mother_function to properly distribute the parameter values.
def master_function(wavegrid, v_resolution=0.0, n_step=25, **kwargs):
    import re

    def process_species(lambdas, f, gamma, b, N, v_rad):
        lambdas = np.atleast_1d(lambdas)
        f = np.atleast_1d(f)
        gamma = np.atleast_1d(gamma)
        b = np.atleast_1d(b)
        N = np.atleast_1d(N)
        v_rad = np.atleast_1d(v_rad)

        n_lines = len(lambdas)
        n_clouds = len(N)
        total = n_lines * n_clouds

        lambdas_use = np.empty(total)
        N_use = np.empty(total)
        v_rad_use = np.empty(total)
        b_use = np.empty(total)
        f_use = np.empty(total)
        gamma_use = np.empty(total)

        for i in range(n_clouds):
            start = i * n_lines
            end = start + n_lines
            lambdas_use[start:end] = lambdas
            N_use[start:end] = N[i]
            v_rad_use[start:end] = v_rad[i]
            b_use[start:end] = b[i]
            f_use[start:end] = f
            gamma_use[start:end] = gamma

        return lambdas_use, N_use, v_rad_use, b_use, f_use, gamma_use

    # Collect species suffixes based on known parameter names ---
    expected_params = ["lambda", "f", "gamma", "b", "N", "v_rad"]
    suffixes = set()

    for key in kwargs:
        for prefix in expected_params:
            if key.startswith(f"{prefix}_"):
                suffix = key[len(prefix) + 1:]  # everything after the underscore
                suffixes.add(suffix)

    if not suffixes:
        raise ValueError("No species parameter groups like 'lambda_1st', 'N_2nd', etc. found.")

    all_lambda = []
    all_N = []
    all_v_rad = []
    all_b = []
    all_f = []
    all_gamma = []

    for suffix in sorted(suffixes, key=lambda x: int(re.findall(r'\d+', x)[0])):
        try:
            lambdas = kwargs[f'lambda_{suffix}']
            f = kwargs[f'f_{suffix}']
            gamma = kwargs[f'gamma_{suffix}']
            b = kwargs[f'b_{suffix}']
            N = kwargs[f'N_{suffix}']
            v_rad = kwargs[f'v_rad_{suffix}']
        except KeyError as e:
            raise ValueError(f"Missing parameter for species '{suffix}': {e}")

        lmb, N_, vr, b_, f_, g_ = process_species(lambdas, f, gamma, b, N, v_rad)

        all_lambda.append(lmb)
        all_N.append(N_)
        all_v_rad.append(vr)
        all_b.append(b_)
        all_f.append(f_)
        all_gamma.append(g_)

    # Combine everything
    lambda_use = np.concatenate(all_lambda)
    N_use = np.concatenate(all_N)
    v_rad_use = np.concatenate(all_v_rad)
    b_use = np.concatenate(all_b)
    f_use = np.concatenate(all_f)
    gamma_use = np.concatenate(all_gamma)

    return mother_function(
        wavegrid=np.asarray(wavegrid),
        lambda0=lambda_use,
        f=f_use,
        gamma=gamma_use,
        b=b_use,
        N=N_use,
        v_rad=v_rad_use,
        v_resolution=v_resolution,
        n_step=n_step
    )

