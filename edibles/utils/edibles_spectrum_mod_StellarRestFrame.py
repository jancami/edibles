import numpy as np
from astropy.io import fits
import astropy.constants as cst
import matplotlib.pyplot as plt
from edibles.edibles import DATADIR
import pandas as pd
from edibles.edibles import PYTHONDIR
from edibles.edibles.utils.edibles_spectrum import EdiblesSpectrum
from edibles.edibles.utils.edibles_oracle import EdiblesOracle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from PyAstronomy import pyasl


class EdiblesSpectrum:
    """
    This class takes a spectrum file from EDIBLES,
    reads the header and data, and creates a DataFrame.

    The class will also contain a set of methods to operate on the data.

    :param filename: Name of the file, starting with the target
    :type filename: str
    :param header: The header of the FITS file from the target observation
    :type header: Object (astropy.io.fits.header.Header)
    :param target: The name of the target
    :type target: str
    :param date: Date of the target observation
    :type date: str
    :param v_bary: Barycentric velocity of the target star
    :type v_bary: float
    :param df: Pandas array containing geocentric and barycentric wavelength, and flux
    :type df: Pandas array (pandas.core.series.Series)
    :param wave: The wavelength grid for the spectrum, geocentric reference frame
    :type wave: Pandas array (pandas.core.series.Series)
    :param wave_units: The units of the wavelength array
    :type wave_units: str
    :param bary_wave: The wavelength grid for the spectrum, barycentric reference frame
    :type bary_wave: Pandas array (pandas.core.series.Series)
    :param flux: The flux data for the spectrum
    :type flux: Pandas array (pandas.core.series.Series)
    :param flux_units: The units of the flux data
    :type flux_units: str

    """

    def __init__(self, filename):
        """
        Filename is relative to the DR3 directory
        """
        self.filename = DATADIR + filename
        self.loadSpectrum()

    def loadSpectrum(self):
        # Assume the file is a DR3 product here.
        hdu = fits.open(self.filename)
        self.header = hdu[0].header
        self.target = self.header["OBJECT"]
        self.date = self.header["DATE-OBS"]

        flux = hdu[0].data
        crval1 = self.header["CRVAL1"]
        cdelt1 = self.header["CDELT1"]
        lenwave = len(flux)
        grid = np.arange(0, lenwave, 1)
        wave = (grid) * cdelt1 + crval1
        self.v_bary = self.header["HIERARCH ESO QC VRAD BARYCOR"]
        bary_wave = wave + (self.v_bary / cst.c.to("km/s").value) * wave

        d = {
            "wave": wave.tolist(),
            "bary_wave": bary_wave.tolist(),
            "flux": flux.tolist(),
        }
        self.df = pd.DataFrame(data=d)

        self.wave = self.df["wave"]
        self.wave_units = "AA"

        self.bary_wave = self.df["bary_wave"]

        self.flux = self.df["flux"]
        self.flux_units = "arbitrary"

    def getSpectrum(self, xmin=None, xmax=None):
        """
        Function to get the wavelength and flux arrays of a particular target.
        If xmin/xmax are not called, the data for the entire spectrum will be returned.

        Args:
            xmin (float): minimum wavelength (Optional)
            xmax (float): Maximum wavelength (Optional)
            bary (bool): Barycentric rest frame, default=False

        Returns:
            ndarray: wavelength grid
            ndarray: flux grid

        """

        if (xmin is not None) and (xmax is not None):
            assert xmin < xmax, "xmin must be less than xmax"

            df_subset = self.df[self.df["bary_wave"].between(xmin, xmax)]

            return df_subset

        return self.df

if __name__ == "__main__":
	pythia = EdiblesOracle()
	print(pythia)
	ListNaweak = pythia.GetObsListByWavelength(3302)
	ListNastrong = pythia.GetObsListByWavelength(5890)
	ListK = pythia.GetObsListByWavelength(7680)
	ListNaweak_O12 = [x for x in ListNaweak if not "O11" in x]
	ListNastrong_O4 = [x for x in ListNastrong if not "O5" in x]
	ListK_O13 = [x for x in ListK if not "O12" in x]
	#print(ListNaweak)
	#print(ListNastrong)
	#print(ListK)
	#print(ListNaweak_O12)
	#print(len(ListNaweak_O12))
	#print(ListNastrong_O4)
	#print(len(ListNastrong_O4))
	#print(ListK_O13)
	#print(len(ListK_O13))
	
	#ISM Frames Database:
	cloudvelocitiesdffilename = PYTHONDIR + '/edibles/data/ISM_Cloud_Velocities_inKms_Approximate_Database_withtargetsasindex.dat'
	cloudvelocitiesdf = pd.read_csv(cloudvelocitiesdffilename, delim_whitespace=True)
	
	#Stellar Rest Frame Targetlist
	ControlledTargetListNoStellarBands678 = ["/HD37061/RED_564/HD37061_w564_redl_20190101_O15.fits", "/HD37061/RED_564/HD37061_w564_redl_20160912_O15.fits", "/HD153919/RED_564/HD153919_w564_redl_20160707_O15.fits", "/HD185859/RED_564/HD185859_w564_redl_20160814_O15.fits", "/HD185859/RED_564/HD185859_w564_redl_20150920_O15.fits", "/HD185859/RED_564/HD185859_w564_redl_20160813_O15.fits", "/HD148937/RED_564/HD148937_w564_redl_20170424_O15.fits", "/HD148937/RED_564/HD148937_w564_redl_20150817_O15.fits", "/HD152424/RED_564/HD152424_w564_redl_20160411_O15.fits", "/HD149404/RED_564/HD149404_w564_redl_20180630_O15.fits", "/HD149404/RED_564/HD149404_w564_redl_20180623_O15.fits", "/HD149404/RED_564/HD149404_w564_redl_20170418_O15.fits", "/HD149404/RED_564/HD149404_w564_redl_20170507_O15.fits", "/HD154043/RED_564/HD154043_w564_redl_20170502_O15.fits", "/HD168076/RED_564/HD168076_w564_redl_20180912_O15.fits", "/HD168076/RED_564/HD168076_w564_redl_20180911_O15.fits", "/HD156201/RED_564/HD156201_w564_redl_20180830_O15.fits", "/HD75860/RED_564/HD75860_w564_redl_20170429_O15.fits", "/HD75860/RED_564/HD75860_w564_redl_20170430_O15.fits", "/HD186841/RED_564/HD186841_w564_redl_20160910_O15.fits", "/HD186841/RED_564/HD186841_w564_redl_20160909_O15.fits", "/HD112272/RED_564/HD112272_w564_redl_20170430_O15.fits", "/HD112272/RED_564/HD112272_w564_redl_20170501_O15.fits", "/HD170938/RED_564/HD170938_w564_redl_20180911_O15.fits", "/HD167971/RED_564/HD167971_w564_redl_20140921_O15.fits", "/HD169454/RED_564/HD169454_w564_redl_20160714_O15.fits", "/HD169454/RED_564/HD169454_w564_redl_20160724_O15.fits", "/HD147889/RED_564/HD147889_w564_redl_20140928_O15.fits"]
	print(len(ControlledTargetListNoStellarBands678))
	print(ControlledTargetListNoStellarBands678)
	
	#Stellar Rest Frames Database:
	stellarrestframedffilename = PYTHONDIR + '/edibles/data/Stellar_Rest_Frame_Velocities_OII5160_withtargetsasindex.txt'
	stellarrestframedf = pd.read_csv(stellarrestframedffilename, delim_whitespace=True)
	print(stellarrestframedf.stellar_rest_velocity)
	
	#THE 2 WEAK NA LINES AROUND 3302A:
	counter = -1
	dataframeOII = pd.DataFrame(columns=['Target', 'OII wavshift in A', 'OII vel in km/s'])
	pd.set_option('display.max_rows', None)
	for filename in ControlledTargetListNoStellarBands678:
		counter = counter + 1
		sp = EdiblesSpectrum(filename)
		#print(sp)
		print(sp.target)
		#print("Barycentric Velocity is", sp.v_bary)
		
		if sp.target not in cloudvelocitiesdf.index:
			cloudvelocityoftarget = 0
			print("De sp.target naam komt niet overeen met je naam in de Cloud Velocity Database (bij verschillende sp.target names per 1 Sightline), hier kan je niks mee, is nu in Barycentric frame")
		else:
			cloudvelocityoftarget = cloudvelocitiesdf.cloud_velocity[sp.target]
			print("Cloud Velocity of {} (in km/s) =".format(sp.target), cloudvelocityoftarget)
		
		OII = 5159.7
		OIIspectrum = sp.getSpectrum(xmin=5150, xmax=5162)
		minflux = np.argmin(OIIspectrum['flux'])
		
		wave=OIIspectrum.wave
		bary_wave = OIIspectrum.bary_wave
		ism_wave = bary_wave + (bary_wave * (cloudvelocityoftarget/(3*10**5)))
		#print("ism_wave", ism_wave)
		
		if sp.target not in stellarrestframedf.index:
			stellarrestvelocityoftarget = 0
			print("Deze Had geen Stellar Line, dus onverschoven")
		else:
			stellarrestvelocityoftarget = stellarrestframedf.stellar_rest_velocity[sp.target]
			print("Stellar Rest Velocity of {} (in km/s) =".format(sp.target), stellarrestvelocityoftarget)
		
		stellar_rest_wave = ism_wave + (ism_wave * (stellarrestvelocityoftarget/(3*10**5)))
		#print("stellar_rest_wave", stellar_rest_wave)
		
		xminvalue = ism_wave[minflux]
		wavshift = OII - xminvalue
		print("wavshift", wavshift)
		
		velocity = wavshift/OII * 3*10**5
		print("velocity", velocity)
		
		dataframeOII.loc[counter] = [sp.target, wavshift, velocity]
		#"""
		plt.axvline(5159.7, color='goldenrod', linewidth = 3.0, label="Strong Stellar Line O-II for Stellar Rest Frame")
		plt.axvline(5155.7, color='darkgoldenrod', linewidth = 3.0, label="Strong Stellar Line Fe-III for Stellar Rest Frame")
		#plt.plot(wave, OIIspectrum["flux"], color='b', label="Geocentric Subset {}".format(sp.target))
		#plt.plot(bary_wave, OIIspectrum["flux"], color='c', label="Barycentric Subset {}".format(sp.target))
		plt.plot(ism_wave, OIIspectrum["flux"], color='k', label="ISM Subset {}".format(sp.target))
		plt.plot(stellar_rest_wave, OIIspectrum["flux"], color='darkgreen', alpha=0.8, label="Stellar Rest Frame Subset {}".format(sp.target))
		plt.legend(loc)
		plt.show()
		#"""
	print("OII")
	print("_")
	print("NOTE: A +/Positive value in wavshift means that the spectrum HAS TO BE shifted to the left/blueshifted to correspond with the Na line. Or differently put, IS shifted to the right/redshifted")
	print("NOTE: Similarly of course, A -/Negative value in wavshift means that the spectrum HAS TO BE shifted to the right/redshifted to correspond with the Na line. Or differently put, IS shifted to the left/blueshifted")
	print("_")
	print(dataframeOII)
	#NOTE VOOR HIERBOVEN QUA WAVSHIFT: 
	#Een +/pos getal correspondeert met een spectrum dat NAAR LINKS verschoven moet worden, DEF 
	#Een -/neg getal correspondeert met een spectrum dat NAAR RECHTS verschoven moet worden, DEF
	
	

	
	#DATA OVER NA en K LINES:
	"""
	#HD37367:
	filename = "/HD37367/BLUE_346/HD37367_w346_blue_20141011_O12.fits"
	#filename = "/HD37367/RED_564/HD37367_w564_redu_20141011_O4.fits"
	#filename = "/HD37367/RED_860/HD37367_w860_redl_20141011_O13.fits"
	
	#HD81188: 
	#filename = "/HD81188/BLUE_346/HD81188_w346_blue_20170418_O12.fits"
	#filename = "/HD81188/RED_564/HD81188_w564_redu_20170418_O4.fits"
	#filename = "/HD81188/RED_860/HD81188_w860_redl_20170418_O13.fits"
	
	#HD164353:
	#filename = "/HD164353/BLUE_346/HD164353_w346_blue_20170615_O12.fits"
	#filename = "/HD164353/RED_564/HD164353_w564_redu_20170615_O4.fits"
	#filename = "/HD164353/RED_860/HD164353_w860_redl_20170414_O13.fits"
	
	sp = EdiblesSpectrum(filename)
	print(sp.target)
	print("Barycentric Velocity is", sp.v_bary)
	#spbeginslicedflux = sp.flux[200:]
	#plt.plot(sp.wave[200:], sp.flux[200:], color = 'k', label="Geocentric")
	
	Na1 = 3302.37
	Na2 = 3302.98
	Na3 = 5889.950
	Na4 = 5895.924
	K5 = 7664.8991
	K6 = 7698.9645
	
	plt.axvline(3302.37, color='c', label="NaI 3302.37 uncontaminated line #1")
	plt.axvline(3302.98, color='m', label="NaI 3302.98 uncontaminated line #2")
	#plt.axvline(5889.950, color='g', label="NaI 5889.950 strong line #1")
	#plt.axvline(5895.924, color='r', label="NaI 5895.924 strong line #2")
	#plt.axvline(7664.8991, color='y', label="KI 7664.8991 strong line #1")
	#plt.axvline(7698.9645, color='c', label="KI 7698.9645 strong line #2")
	plt.legend()
	
	subset = sp.getSpectrum(xmin=3301.8, xmax=3302.7)
	subset = sp.getSpectrum(xmin=3302.7, xmax=3303.6)
	#subset = sp.getSpectrum(xmin=5888.5, xmax=5891.0)
	#subset = sp.getSpectrum(xmin=5895.5, xmax=5897.0)
	#subset = sp.getSpectrum(xmin=7663.5, xmax=7665.7)
	#subset = sp.getSpectrum(xmin=7697.5, xmax=7700.0)
	
	minflux = np.argmin(subset['flux'])
	xminvalue = subset.bary_wave[minflux]
	print(xminvalue)
	print("Voor Interstellar Frame moet Spectrum {} naar Links".format(xminvalue-Na1))
	print("Cloud Velocity hierbij is".format())
	
	#plt.plot(subset["wave"], subset["flux"], color='k', label="Geocentric Subset {}".format(filename[1:9]))
	plt.plot(subset["bary_wave"], subset["flux"], color='b', label="Barycentric Subset {}".format(filename[1:9]))
	plt.legend()
	
	plt.show()
	"""

