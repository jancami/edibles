import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import inspect
import collections
from math import floor
from lmfit import Model
from lmfit.models import update_param_vals

from edibles.utils.voigt_profile import voigt_absorption_line
from edibles.models import ContinuumModel

from pathlib import Path
from edibles import DATADIR
from edibles import PYTHONDIR

class ISLineFitter():
    def __init__(self, wave, flux, SNR=None, v_resolution=3.0):
        # wave and flux from edibles spectrum or elsewhere
        # load by separate method to be built
        self.wave = wave
        self.flux = flux
        self.SNR = SNR # Don't know yet how to use SNR in LMFIT
        self.v_res = v_resolution

        # attribute to archive model-fitting history
        self.model_all = []
        self.model_old = None # model with n components
        self.model_new = None # model with n+1 components
        
        # read in atomic line data frame
        folder = Path(PYTHONDIR+"/data")
        filename = folder / "auxiliary_data" / "line_catalogs" / "edibles_linelist_atoms.csv"
        self.species_df=pd.read_csv(filename)

    def getData2Fit(self, lam_0, windowsize):
        # clip and return spectral data around target lines
        return wave2fit, flux2fit, SNR2fit

    def baysianCriterion(self):
        # do the Baysian analysis to determine if self.model_new is better than self.model_old
        # if pass: keep adding new components
        #          self.model_old = self.model_new,
        #          return False
        # if not pass: return True to exit
        pass

    def fit(self, species="KI", n_anchors=5,Wave=None, WaveMin=None, WaveMax=None, OscillatorStrength=None, OscillatorStrengthMin=None, OscillatorStrengthMax=None,Gamma=None, GammaMin=None, GammaMax=None):
        # Do the fitting
        # 1. Get atomic data using Heather's method
        # 2. Clip spectral data around target lines
        # 3. Build model to be fit
        # 4. Fit, compare, repeat!

        ######################
        # get lam_0, fjj, gamma from Heather's code.
        # I still think get Nmag from data table, rather than estimating it,
        # would be a good idea...
        # For now use default, already embedded in ISLineModel
        
        spec_name,lam_0,fjj,gamma=self.select_species_data(species=species,Wave=Wave, WaveMin=WaveMin, WaveMax=WaveMax, OscillatorStrength=OscillatorStrength, OscillatorStrengthMin=OscillatorStrengthMin, OscillatorStrengthMax=OscillatorStrengthMax,Gamma=Gamma, GammaMin=GammaMin, GammaMax=GammaMax)
       
        ######################

        ######################
        # get wave2fit, flux2fit, and possibily SNR2fit
        # Require new method, self.getData2Fit(lam_0, windowsize)
        ######################
        wave2fit, flux2fit, SNR2fit = self.getData2Fit(lam_0, windowsize=5)

        while True:
            model2fit, pars_guess = self.buildModel(lam_0, fjj, gamma, Nmag, n_anchors)
            result = model2fit.fit(data=flux2fit, params=pars_guess, x=wave2fit)
            self.model_all.append(model2fit)
            self.model_new = model2fit
            # should we append parameters rather than models? Check how Baysian works...
            if self.baysianCriterion():
                break

        return result.params

    def buildModel(self, lam_0, fjj, gamma, Nmag, n_anchors, n_components = None):
        # build continuum and line model then combine them
        # we can reuse continuum model to boost efficiency?

        # Continuum first
        continuum_model = ContinuumModel(n_anchors = n_anchors)
        pars_guess = continuum_model.guess(self.flux, self.wave)
        model2fit = continuum_model

        # check n_components before building line-model
        # only include
        if n_components is None:
            n_components = len(self.model_all)
        if n_components > 0:
            line_model = ISLineModel(n_components,
                                     lam_0 = lam_0,
                                     fjj = fjj,
                                     gamma = gamma,
                                     Nmag = Nmag,
                                     v_res = self.v_res)

            ##############################
            # to do:
            # Guessing v_offs form the residual of self.model_old using corr method
            # Not available now, I'll just say V_off = 0s
            #
            # If we do not load Nmag from data table, we can make model.guess() more
            # complex to determine the N from residual of self.model_old?
            ##############################
            V_off = [0.0] * n_components
            pars_guess.update(line_model.guess(V_off=V_off))
            model2fit = model2fit * line_model

        return model2fit, pars_guess

    def plotModel(self):
        # method to make plots
        pass
        
        
    def select_species_data(self,species=None,Wave=None, WaveMin=None, WaveMax=None, OscillatorStrength=None, OscillatorStrengthMin=None, OscillatorStrengthMax=None,Gamma=None, GammaMin=None, GammaMax=None):
        '''This method will provide a filtered list of species information that matches
        the specified criteria on sightline/target parameters as well as
        on observational criteria (e.g. wavelength range).
    
         '''
        
        bool_species_matches = np.zeros(len(self.species_df.index),dtype=bool)
        
        if species is None:
            bool_species_matches = np.ones(len(self.species_df.index),dtype=bool)
        elif (isinstance(species, np.ndarray) | isinstance(species, list)):
            
            for thisobject in species:

                bool_species_matches = (self.species_df.Species.str.contains(thisobject)) | (bool_species_matches)
                
        else:
            
            bool_species_matches = self.species_df.Species == species

        bool_wave_matches = np.ones(len(self.species_df.index),dtype=bool)
        if Wave:
            bool_wave_matches = (self.species_df.WavelengthAir == Wave)
        if WaveMin:
            bool_wave_matches = (self.species_df.WavelengthAir > WaveMin) & (bool_wave_matches)
        if WaveMax:
            bool_wave_matches = (self.species_df.WavelengthAir < WaveMax) & (bool_wave_matches)
            
        
        bool_osc_matches = np.ones(len(self.species_df.index),dtype=bool)
        if OscillatorStrength:
            bool_osc_matches = (self.species_df.OscillatorStrength == OscillatorStrength)
        if OscillatorStrengthMin:
            bool_osc_matches = (self.species_df.OscillatorStrength > OscillatorStrengthMin) & (bool_osc_matches)
        if OscillatorStrengthMax:
            bool_osc_matches = (self.species_df.OscillatorStrength < OscillatorStrengthMax) & (bool_osc_matches)
            
        
        bool_gamma_matches = np.ones(len(self.species_df.index),dtype=bool)
        if Gamma:
            bool_gamma_matches = (self.species_df.Gamma == Gamma)
        if GammaMin:
            bool_gamma_matches = (self.species_df.Gamma > GammaMin) & (bool_gamma_matches)
        if GammaMax:
            bool_gamma_matches = (self.species_df.Gamma < GammaMax) & (bool_gamma_matches)
            
            
            
        ind = np.where(bool_species_matches & bool_wave_matches & bool_osc_matches & bool_gamma_matches)[0]
        self.species_list=self.species_df['Species'].iloc[ind].to_numpy()
        
        self.air_wavelength=self.species_df['WavelengthAir'].iloc[ind].to_numpy()
        self.oscillator_strength=self.species_df['OscillatorStrength'].iloc[ind].to_numpy()
        self.gamma=self.species_df['Gamma'].iloc[ind].to_numpy()
        return(self.species_list,self.air_wavelength,self.oscillator_strength,self.gamma)

        
    
    def determine_vrad_from_correlation(self,wave, flux, model):
        """
        Function to calculate the correlation between an observed spectrum and a model as a function of
        radial velocity and return the radial velocity with the highest correlation coefficient. 
        Args:
            wave (float64): array of wavelengths
            flux (float64): Flux (observed)
            model(float64): model
        Returns:
            vrad_best: radial velocity corresponding to highest correlation. 
        """
        # Create the grid of velocities at which to calculate the correlation. 
        # Using a step size of 0.1 km/s here, and a range of -50 to 50; this should 
        # suffice for most sightlines. 
        v_rad_grid = np.arange(-50.,50.,.1) # in km / s
        all_corr = v_rad_grid * 0.
        #print(v_rad_grid)
        for loop in range(len(v_rad_grid)):
            v_rad = v_rad_grid[loop]
            Doppler_factor = 1. + v_rad / cst.c.to("km/s").value
            #print(Doppler_factor)
            new_wave = wave * Doppler_factor
            
            # Interpolate shifted model to original wavelength grid
            interpolationfunction = interp1d(
            new_wave, model, kind="cubic", fill_value="extrapolate")
            interpolatedModel = interpolationfunction(wave)
            
            # Calculate correlation coefficient
            this_c, _ = pearsonr(flux, interpolatedModel)
            all_corr[loop] = this_c
            # Return the radial velocity at the maximum correlation.  
        v_rad_best = v_rad_grid[np.argmax(all_corr)]    

        return v_rad_best



class ISLineModel(Model):
    def __init__(self, n_components,
                 lam_0 = [3302.369, 3302.978],
                 fjj = [8.26e-03, 4.06e-03],
                 gamma = [6.280e7, 6.280e7],
                 Nmag = 14,
                 v_res = 3.0,
                 independent_vars=["x"],
                 prefix="",
                 nan_policy="raise",
                 n_step=25,
                 **kwargs):
        """
        :param n_components: int, number of velocity components
        :param lam_0: list, air wavelength of target line. lam_0, fjj and gamma should have same length
        :param fjj: list, oscillator strengths
        :param gamma: list, Gamma parameter related to broadening.
        :param Nmag: int, magnitude of column density
        :param v_res: float, resolution in km/s
        :param independent_vars: from lmfit and Klay's code
        :param prefix: from lmfit and Klay's code
        :param nan_policy: from lmfit and Klay's code
        :param n_step: int, no. of points in 1*FWHM during calculation. Under-sample losses information
        but over-sample losses efficiency.
        :param kwargs: ???
        """
        self.n_components, self.lam_0, self.fjj, self.gamma, self.Nmag, self.n_setp = \
            self.__inputCheck(n_components, lam_0, fjj, gamma, Nmag, n_step)
        self.v_res = v_res

        kwargs.update({"prefix": prefix,
                       "nan_policy": nan_policy,
                       "independent_vars": independent_vars})

        # creat b, N, and V according to no. of components
        self.b_names = ["b_Cloud%i" % (i) for i in range(n_components)]
        self.N_names = ["N_Cloud%i" % (i) for i in range(n_components)]
        self.V_names = ["V_off_Cloud%i" % (i) for i in range(n_components)]
        kwargs["param_names"] = self.b_names + self.N_names + self.V_names
        params = {}
        for name in kwargs["param_names"]:
            if name[0] == "b":
                params[name] = 1.0
            if name[0] == "N":
                params[name] = 1.0
            if name[0] == "V":
                params[name] = 0.0

        # Other than the default parameters from Cloud0, other parameters are passed in kwargs
        def calcISLineModel(x, b_Cloud0=1.0, N_Cloud0=1.0, V_off_Cloud0=0.0, **kwargs):
            lambda0 = self.lam_0 * self.n_components
            f = self.fjj * self.n_components
            gamma = self.gamma * self.n_components
            N_mag = self.Nmag
            v_resolution = self.v_res

            # parse parameters
            bs = [b_Cloud0] * len(self.lam_0)
            Ns = [N_Cloud0 * 10 ** N_mag] * len(self.lam_0)
            V_offs = [V_off_Cloud0] * len(self.lam_0)

            # just something to print so we know it's working...
            print("Velocity of Cloud_0:", np.unique(V_offs))

            for name in kwargs.keys():
                if name[0] == "b":
                    bs = bs + [kwargs[name]] * len(self.lam_0)
                if name[0] == "N":
                    Ns = Ns + [kwargs[name] * 10 ** N_mag] * len(self.lam_0)
                if name[0] == "V":
                    V_offs = V_offs + [kwargs[name]] * len(self.lam_0)

            # no problem for convolution since we only call voigt_absorption_line once.
            flux = voigt_absorption_line(
                x,
                lambda0=lambda0,
                b=bs,
                N=Ns,
                f=f,
                gamma=gamma,
                v_rad=V_offs,
                v_resolution=v_resolution,
                n_step=self.n_setp
            )

            return flux

        # play with the signature of the calculation function
        # I don't really understand what is happening...
        sig = inspect.signature(calcISLineModel)
        base_b = inspect.signature(calcISLineModel).parameters["b_Cloud0"]
        base_N = inspect.signature(calcISLineModel).parameters["N_Cloud0"]
        base_V_off = inspect.signature(calcISLineModel).parameters["V_off_Cloud0"]

        d = {'x': sig.parameters['x']}
        for i in range(n_components):
            b_key = "b_Cloud" + str(i)
            b_val = base_b.replace(name=b_key)
            d[b_key] = b_val

            N_key = "N_Cloud" + str(i)
            N_val = base_N.replace(name=N_key)
            d[N_key] = N_val

            V_off_key = "V_off_Cloud" + str(i)
            V_off_val = base_V_off.replace(name=V_off_key)
            d[V_off_key] = V_off_val

        d = collections.OrderedDict(d)
        calcISLineModel.__signature__ = sig.replace(parameters=tuple(d.values()))

        super().__init__(calcISLineModel, **kwargs)

    def guess(self, V_off=[0.0], **kwargs):
        # For now just type in V_off but we can include v_correlation in the future

        assert len(V_off) == self.n_components, "Number of components do not match."

        pars = self.make_params()
        for i, v in enumerate(V_off):
            pars["%sb_Cloud%i" % (self.prefix, i)].set(value=1.0, min=0, max=10)
            pars["%sN_Cloud%i" % (self.prefix, i)].set(value=1.0, min=0, max=1000)
            pars["%sV_off_Cloud%i" % (self.prefix, i)].set(value=v, min=v-50, max=v+50)
            # we can further constrain min and max on V_off if we have good estimate.

        return update_param_vals(pars, self.prefix, **kwargs)

    def __inputCheck(self, n_components, lam_0, fjj, gamma, Nmag, n_step):
        # n_components should be int
        if not isinstance(n_components, int):
            raise TypeError("n_components (%.1f) must be an integer!" % n_components)

        # lam_0, fjj, gamma should be LIST have same length
        if not isinstance(lam_0, list) or not isinstance(fjj, list) or not isinstance(gamma, list):
            raise TypeError("lam_0, fjj, gamma should be list")
        len_array = [len(lam_0), len(fjj), len(gamma)]
        if not np.max(len_array) == np.max(len_array):
            raise TypeError("lam_0, fjj, gamma should have the same length")

        # Nmag and n_step should be int but just in case it's a float
        Nmag = floor(Nmag)
        n_step = floor(n_step)
        return n_components, lam_0, fjj, gamma, Nmag, n_step



if __name__ == "__main__":
    print("Hello Word!")

    ##Random values for flux and wave to init. Remove before public push, or can be changed by others.
    import random
    wave=np.linspace(0,10,11)
    flux=np.asarray((random.sample(range(100), k=len(wave))))/100
    ####################################################
    fit_test=ISLineFitter(wave,flux)
    test_species_info=fit_test.select_species_data(OscillatorStrengthMin=0.1)
    #fit_test.fit(species=['Na'])


