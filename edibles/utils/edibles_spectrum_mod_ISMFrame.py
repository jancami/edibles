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
	
	cloudvelocitiesdffilename = PYTHONDIR + '/edibles/data/ISM_Cloud_Velocities_inKms_Approximate_Database_withtargetsasindex.dat'
	cloudvelocitiesdf = pd.read_csv(cloudvelocitiesdffilename, delim_whitespace=True)
	
	#THE 2 WEAK NA LINES AROUND 3302A:
	counter = -1
	dataframeNaweak = pd.DataFrame(columns=['Target', 'Na wl1  in A', 'Na wl2 in A', 'Na wl1 vel in km/s', 'Na wl2 vel in km/s', 'Average Vel in km/s', '~wavshift in A around C6H Bands 1,2'])
	pd.set_option('display.max_rows', None)
	for filename in ListK_O13:
		counter = counter + 1
		sp = EdiblesSpectrum(filename)
		#print(sp)
		#print(sp.target)
		#print("Barycentric Velocity is", sp.v_bary)
		
		if sp.target not in cloudvelocitiesdf.index:
			cloudvelocityoftarget = 0
			print("De sp.target naam komt niet overeen met je naam in de Cloud Velocity Database (bij verschillende sp.target names per 1 Sightline), hier kan je niks mee, is nu in Barycentric frame")
		else:
			cloudvelocityoftarget = cloudvelocitiesdf.cloud_velocity[sp.target]
			print("Cloud Velocity of {} (in km/s) =".format(sp.target), cloudvelocityoftarget)
		
		Na1 = 3302.37
		subset1 = sp.getSpectrum(xmin=3301.8, xmax=3302.7)
		minflux1 = np.argmin(subset1['flux'])
		xminvalue1 = subset1.bary_wave[minflux1]
		wavshift1 = Na1 - xminvalue1
		
		bary_wave1 = subset1.bary_wave
		#print("bary_wave1", bary_wave1)
		ism_wave1 = bary_wave1 + (bary_wave1 * (cloudvelocityoftarget/(3*10**5)))
		#print("ism_wave1", ism_wave1)
		
		Na2 = 3302.98
		subset2 = sp.getSpectrum(xmin=3302.8, xmax=3303.6)
		minflux2 = np.argmin(subset2['flux'])
		xminvalue2 = subset2.bary_wave[minflux2]
		wavshift2 = Na2 - xminvalue2
		
		bary_wave2 = subset2.bary_wave
		#print("bary_wave2", bary_wave2)
		ism_wave2 = bary_wave2 + (bary_wave2 * (cloudvelocityoftarget/(3*10**5)))
		#print("ism_wave2", ism_wave2)
		
		vel1 = wavshift1/Na1 * 3*10**5
		vel2 = wavshift2/Na2 * 3*10**5
		
		if abs(wavshift1-wavshift2) <= 0.10:
			averagevel = (vel1 + vel2)/2
			C6HBands12wavshift = 5265.0 * (averagevel/(3*10**5))
		else:
			averagevel = "too far apart"
			C6HBands12wavshift = "x"
		
		dataframeNaweak.loc[counter] = [sp.target, wavshift1, wavshift2, vel1, vel2, averagevel, C6HBands12wavshift]
		"""
		plt.axvline(3302.37, color='c', label="NaI 3302.37 uncontaminated line #1")
		plt.axvline(3302.98, color='m', label="NaI 3302.98 uncontaminated line #2")
		plt.plot(subset1["bary_wave"], subset1["flux"], color='b', label="Barycentric Subset Links {}".format(sp.target))
		plt.plot(ism_wave1, subset1["flux"], color='r', label="ISM Subset {}".format(sp.target))
		plt.legend()
		plt.show()
		plt.axvline(3302.37, color='c', label="NaI 3302.37 uncontaminated line #1")
		plt.axvline(3302.98, color='m', label="NaI 3302.98 uncontaminated line #2")
		plt.plot(subset2["bary_wave"], subset2["flux"], color='b', label="Barycentric Subset Rechts {}".format(sp.target))
		plt.plot(ism_wave2, subset2["flux"], color='r', label="ISM Subset {}".format(sp.target))
		plt.legend()
		plt.show()
		"""
	print("The 2 Weak Na Lines ~3302")
	print("_")
	print("NOTE: A +/Positive value in wavshift means that the spectrum HAS TO BE shifted to the left/blueshifted to correspond with the Na line. Or differently put, IS shifted to the right/redshifted")
	print("NOTE: Similarly of course, A -/Negative value in wavshift means that the spectrum HAS TO BE shifted to the right/redshifted to correspond with the Na line. Or differently put, IS shifted to the left/blueshifted")
	print("_")
	print(dataframeNaweak)
	#NOTE VOOR HIERBOVEN QUA WAVSHIFT: 
	#Een +/pos getal correspondeert met een spectrum dat NAAR LINKS verschoven moet worden, DEF 
	#Een -/neg getal correspondeert met een spectrum dat NAAR RECHTS verschoven moet worden, DEF
	
	
	
	#THE 2 STRONG NA LINES AROUND 5900A:
	counter = -1
	dataframeNastrong = pd.DataFrame(columns=['Target', 'Na sl1  in A', 'Na sl2 in A', 'Na sl1 vel in km/s', 'Na sl2 vel in km/s', 'Average Vel in km/s', '~wavshift in A around C6H Bands 1,2'])
	for filename in ListNastrong_O4:
		counter = counter + 1
		sp = EdiblesSpectrum(filename)
		#print(sp.target)
		#print("Barycentric Velocity is", sp.v_bary)
		
		if sp.target not in cloudvelocitiesdf.index:
			cloudvelocityoftarget = 0
			print("De sp.target naam komt niet overeen met je naam in de Cloud Velocity Database (bij verschillende sp.target names per 1 Sightline), hier kan je niks mee, is nu in Barycentric frame")
		else:
			cloudvelocityoftarget = cloudvelocitiesdf.cloud_velocity[sp.target]
			print("Cloud Velocity of {} (in km/s) =".format(sp.target), cloudvelocityoftarget)
		
		Na3 = 5889.950
		subset3 = sp.getSpectrum(xmin=5888.5, xmax=5892.0)
		minflux3 = np.argmin(subset3['flux'])
		xminvalue3 = subset3.bary_wave[minflux3]
		wavshift3 = Na3 - xminvalue3
		
		bary_wave3 = subset3.bary_wave
		print("bary_wave3", bary_wave3)
		ism_wave3 = bary_wave3 + (bary_wave3 * (cloudvelocityoftarget/(3*10**5)))
		print("ism_wave3", ism_wave3)
		
		Na4 = 5895.924
		subset4 = sp.getSpectrum(xmin=5894.5, xmax=5897.0)
		minflux4 = np.argmin(subset4['flux'])
		xminvalue4 = subset4.bary_wave[minflux4]
		wavshift4 = Na4 - xminvalue4
		
		bary_wave4 = subset4.bary_wave
		print("bary_wave4", bary_wave4)
		ism_wave4 = bary_wave4 + (bary_wave4 * (cloudvelocityoftarget/(3*10**5)))
		print("ism_wave4", ism_wave4)
		
		vel3 = wavshift3/Na3 * 3*10**5
		vel4 = wavshift4/Na4 * 3*10**5
		
		if abs(wavshift3-wavshift4) <= 0.10:
			averagevel = (vel3 + vel4)/2
			C6HBands12wavshift = 5265.0 * (averagevel/(3*10**5))
		else:
			averagevel = "too far apart"
			C6HBands12wavshift = "x"
		
		dataframeNastrong.loc[counter] = [sp.target, wavshift3, wavshift4, vel3, vel4, averagevel, C6HBands12wavshift]
		
		#"""
		plt.axvline(5889.950, color='c', label="NaI 5889.950 strong line #1")
		plt.axvline(5895.924, color='m', label="NaI 5895.924 strong line #2")
		plt.plot(subset3["bary_wave"], subset3["flux"], color='b', label="Barycentric Subset {}".format(sp.target))
		plt.plot(ism_wave3, subset3["flux"], color='r', label="ISM Subset {}".format(sp.target))
		plt.legend()
		plt.show()
		plt.axvline(5889.950, color='c', label="NaI 5889.950 strong line #1")
		plt.axvline(5895.924, color='m', label="NaI 5895.924 strong line #2")
		plt.plot(subset4["bary_wave"], subset4["flux"], color='b', label="Barycentric Subset {}".format(sp.target))
		plt.plot(ism_wave4, subset4["flux"], color='r', label="ISM Subset {}".format(sp.target))
		plt.legend()
		plt.show()
		#"""
	print("The 2 Strong Na Lines ~5900")
	print("_")
	print("NOTE: A +/Positive value in wavshift means that the spectrum HAS TO BE shifted to the left/blueshifted to correspond with the Na line. Or differently put, IS shifted to the right/redshifted")
	print("NOTE: Similarly of course, A -/Negative value in wavshift means that the spectrum HAS TO BE shifted to the right/redshifted to correspond with the Na line. Or differently put, IS shifted to the left/blueshifted")
	print("_")
	print(dataframeNastrong)
	#NOTE VOOR HIERBOVEN QUA WAVSHIFT: 
	#Een +/pos getal correspondeert met een spectrum dat NAAR LINKS verschoven moet worden, DEF 
	#Een -/neg getal correspondeert met een spectrum dat NAAR RECHTS verschoven moet worden, DEF
	
	
	
	#THE 2 K LINES AROUND 7680A:
	counter = -1
	dataframeK = pd.DataFrame(columns=['Target', 'K1 in A', 'K2 in A', 'K1 vel in km/s', 'K2 vel in km/s', 'Average Vel in km/s', '~wavshift in A around C6H Bands 1,2'])
	for filename in ListK_O13:
		counter = counter + 1
		sp = EdiblesSpectrum(filename)
		#print(sp.target)
		#print("Barycentric Velocity is", sp.v_bary)
		
		if sp.target not in cloudvelocitiesdf.index:
			cloudvelocityoftarget = 0
			print("De sp.target naam komt niet overeen met je naam in de Cloud Velocity Database (bij verschillende sp.target names per 1 Sightline), hier kan je niks mee, is nu in Barycentric frame")
		else:
			cloudvelocityoftarget = cloudvelocitiesdf.cloud_velocity[sp.target]
			print("Cloud Velocity of {} (in km/s) =".format(sp.target), cloudvelocityoftarget)
		
		K5 = 7664.8991
		subset5 = sp.getSpectrum(xmin=7663.5, xmax=7665.7)
		minflux5 = np.argmin(subset5['flux'])
		xminvalue5 = subset5.bary_wave[minflux5]
		wavshift5 = K5 - xminvalue5
		wavshift5 = 'x'
		vel5 = 'x'
		
		bary_wave5 = subset5.bary_wave
		print("bary_wave5", bary_wave5)
		ism_wave5 = bary_wave5 + (bary_wave5 * (cloudvelocityoftarget/(3*10**5)))
		print("ism_wave5", ism_wave5)
		
		K6 = 7698.9645
		subset6 = sp.getSpectrum(xmin=7697.5, xmax=7700.0)
		minflux6 = np.argmin(subset6['flux'])
		xminvalue6 = subset6.bary_wave[minflux6]
		wavshift6 = K6 - xminvalue6
		
		bary_wave6 = subset6.bary_wave
		print("bary_wave6", bary_wave6)
		ism_wave6 = bary_wave6 + (bary_wave6 * (cloudvelocityoftarget/(3*10**5)))
		print("ism_wave6", ism_wave6)
		
		vel6 = wavshift6/K6 * 3*10**5
		print(sp.target)
		print(vel6)
		
		averagevel = 'x'
		C6HBands12wavshift = 5265.0 * (vel6/(3*10**5))
		
		dataframeK.loc[counter] = [sp.target, wavshift5, wavshift6, vel5, vel6, averagevel, C6HBands12wavshift]
		"""
		#plt.axvline(7664.8991, color='y', label="KI 7664.8991 strong line #1")
		#plt.axvline(7698.9645, color='c', label="KI 7698.9645 strong line #2")
		#plt.plot(subset5["bary_wave"], subset5["flux"], color='b', label="Barycentric Subset L {}".format(filename[1:9]))
		#plt.legend()
		#plt.show()
		#plt.axvline(7664.8991, color='y', label="KI 7664.8991 strong line #1")
		plt.axvline(7698.9645, color='c', label="KI 7698.9645 strong line #2")
		plt.plot(subset6["bary_wave"], subset6["flux"], color='b', label="Barycentric Subset R {}".format(filename[1:9]))
		plt.legend()
		plt.show()
		"""
	print("The 2 K Lines ~7680")
	print("_")
	print("NOTE: A +/Positive value in wavshift means that the spectrum HAS TO BE shifted to the left/blueshifted to correspond with the Na line. Or differently put, IS shifted to the right/redshifted")
	print("NOTE: Similarly of course, A -/Negative value in wavshift means that the spectrum HAS TO BE shifted to the right/redshifted to correspond with the Na line. Or differently put, IS shifted to the left/blueshifted")
	print("_")
	print(dataframeK)
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

