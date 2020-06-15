from edibles.edibles import DATADIR
from edibles.edibles import PYTHONDIR
from edibles.edibles.utils.edibles_spectrum import EdiblesSpectrum
from edibles.edibles.utils.edibles_oracle import EdiblesOracle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from PyAstronomy import pyasl



#OBTAINING THE LABSPECTRUM/SIMULATED SPECTRUM OF A MOLECULE:
# Example: Look at Farid's spectrum of pentacene. 
###labfilename = PYTHONDIR + '/edibles/data/Labdata/CRDS/PENTACENE.DAT'

#Nieuwe Main C6H 19bands Simulated Data 20k, 5x more smooth/unresolved:
labfilename = PYTHONDIR + '/edibles/data/Labdata/C6H/C6H_Data_All19Bands_20K_Resolvingpower5timeslesstoGau03.dat'

#Main C6H 19Bands Simulated Data 20K:
#labfilename = PYTHONDIR + '/edibles/data/Labdata/C6H/C6H_Data_All19Bands_20K.dat'

#Temperature Changes for C6H Spectrum:
#labfilename = PYTHONDIR + '/edibles/data/Labdata/C6H/C6H_Data_All19Bands_2K.dat'
#labfilename = PYTHONDIR + '/edibles/data/Labdata/C6H/C6H_Data_All19Bands_5K.dat'
#labfilename = PYTHONDIR + '/edibles/data/Labdata/C6H/C6H_Data_All19Bands_10K.dat'
#labfilename = PYTHONDIR + '/edibles/data/Labdata/C6H/C6H_Data_All19Bands_15K.dat'
#labfilename = PYTHONDIR + '/edibles/data/Labdata/C6H/C6H_Data_All19Bands_20K.dat'
#labfilename = PYTHONDIR + '/edibles/data/Labdata/C6H/C6H_Data_All19Bands_25K.dat'
#labfilename = PYTHONDIR + '/edibles/data/Labdata/C6H/C6H_Data_All19Bands_50K.dat'
#labfilename = PYTHONDIR + '/edibles/data/Labdata/C6H/C6H_Data_All19Bands_200K.dat'
#labfilename = PYTHONDIR + '/edibles/data/Labdata/C6H/C6H_Data_All19Bands_2000K.dat'

#C8H Data 2 bands:
#labfilename = PYTHONDIR + '/edibles/data/Labdata/C8H/C8H_Bands1_2_15K.dat'

#C10H Data 1 band:
#labfilename = PYTHONDIR + '/edibles/data/Labdata/C10H/C10H_Band1_15K.dat'

#TRANSFORMING THE SIMULATED MOLECULE SPECTRUM/DATA TO COMPARE WITH THE EDIBLES DATA:
# The file contains wavenumbers (in vacuum) and intensities. Read those in as pandas dataframe. 
labspec = pd.read_csv(labfilename, delim_whitespace=True)
# Add a column to contain the (air) wavelength in AA. 
labspec['wavelength'] = pyasl.vactoair2(1e8/labspec['wno'],mode='ciddor')
normint = labspec.int - np.median(labspec.int)
#Doe de y-normalisatie naar beneden, want het zijn absorptielijnen, dus *-1 hieronder
labspec['norm'] = normint / np.max(normint) * -1

#PROOF IN TERMINAL THAT THE CONVERSION FROM WNO IN VAC > WAV IN AIR GOED GAAT:
#print("Labspec in Wavenumbers", labspec['wno'])
#print("Labspec in Wavelengths in Vacuum", 1e8/labspec['wno'])
#print("Labspec in Wavelength in Air", labspec['wavelength'])




#EDIBLES ORACLE TARGET LIST: Om een lijst te obtainen van alle Sightlines in een bepaalde
#Wavelength range. Alleen nuttig voor een initiele search through!!!.
#C6H LIST FOR BANDS 1,2:
#"""
# Pentacene has strongest bands at about 5338 and 5361 AA. Let's search what EDIBLES spectra
# We have in that wavelength range and compare. 
#plotrange = [5260,5270]
plotrange = [5230,5300]
pythia = EdiblesOracle()
List = pythia.GetObsListByWavelength(5265.5)
#Reversing the list so that the Absorption Peaks around above wavelength come first
#Sclice above List to only contain O17's, this gives a much clearer view. Could in the future also
#be manipulated to only contain O18's for example, or remove it to contain both.
ListO17 = [x for x in List if not "O18" in x]
#print(ReverseListO17)


#"""
#C6H LIST FOR BANDS 3,4,5:
"""
#Om de 3e, 4e en 5e band tegelijkertijd te kunnen vergelijken:
plotrange = [5190, 5250]
pythia = EdiblesOracle()
List2 = pythia.GetObsListByWavelength(5220)
ReverseList2 = List2[::-1]
ReverseListO16 = [x for x in ReverseList2 if not "O17" in x]
print(ReverseListO16)
counter = 0
"""
#C6H LIST FOR BAND 19:
"""
#Om alle C6H Band 19's te kunnen vergelijken, dit lijkt based on Xavier's intensity strengths
#Paper de 2e sterkste band te zijn na Bands 1,2!!!. Dus de beste vergelijking mocht er
#Iets zijn met bands 1,2 en ook een goede check daarvoor!!!.
plotrange = [4740, 4750]
pythia = EdiblesOracle()
ListBand19 = pythia.GetObsListByWavelength(4745)
ReverseListBand19 = ListBand19[::-1]
ReverseListBand19_O4 = [x for x in ReverseListBand19 if not "O5" in x]
ReverseListBand19_O4x = [x for x in ReverseListBand19_O4 if not "O26" in x]
ReverseListBand19_O4xx = [x for x in ReverseListBand19_O4x if not "O27" in x]
#print(ReverseListBand19_O4xx)
counter = 0
"""
#C8H List:
"""
#C8H Main Bands List: 
plotrange = [6250, 6270]
pythia = EdiblesOracle()
ListC8H = pythia.GetObsListByWavelength(6262)
ReverseListC8H = ListC8H[::-1]
ReverseListC8HO10 = [x for x in ReverseListC8H if not "O11" in x]
print(ReverseListC8HO10)
counter = 0
"""

#C10H List:
"""
#C10H Main Bands List: 
plotrange = [7120, 7160]
pythia = EdiblesOracle()
ListC8H = pythia.GetObsListByWavelength(7141)
ReverseListC10H = ListC8H[::-1]
ReverseListC10HO7 = [x for x in ReverseListC10H if not "O6" in x]
print(ReverseListC10HO7)
counter = 0
"""



#CLOUD VELOCITIES DATAFRAME INTEGRATED:
#Old:
"""
HDtargetlist=['HD37367', 'HD81188']
cloudvelocitylist=[-15, -8]
testcloudvelocitydataframe = pd.DataFrame(cloudvelocitylist, columns=['Cloud_Velocity'], index=HDtargetlist)
print(testcloudvelocitydataframe)
cloudvelocity = testcloudvelocitydataframe.Cloud_Velocity['HD37367']
"""
#Not Handy:
"""
cloudvelocitiesfilename = PYTHONDIR + '/edibles/data/ISM_Cloud_Velocities_inKms_Approximate_Database.dat'
cloudvelocities = pd.read_csv(cloudvelocitiesfilename, delim_whitespace=True)
print(cloudvelocities)
print(cloudvelocities.target)
print(cloudvelocities.cloud_velocity)
print(cloudvelocities.comment)
"""

plotrange = [5255,5270]
pythia = EdiblesOracle()
List = pythia.GetObsListByWavelength(5265.5)
ListO17 = [x for x in List if not "O18" in x]
print(List)
print(ListO17)
print(len(ListO17))
SelectedListO17 = [x for x in ListO17 if "183143" in x or "147889" in x or "169454" in x or "167971" in x or "170938" in x or "186841" in x or "164740" in x or "186745" in x or "75860" in x or "156201" in x or "112272" in x or "61827" in x or "165319" in x or "168076" in x or "154043" in x or "73882" in x or "149404" in x or "152424" in x or "148937" in x or "161056" in x or "80558" in x or "185859" in x or "167838" in x or "153919" in x or "43384" in x or "37061" in x]
#SelectedListO17 = [x for x in ListO17 if "183143" in x or "147889" in x or "169454" in x or "167971" in x or "170938" in x or "186841" in x or "164740" in x or "186745" in x or "75860" in x or "156201" in x or "112272" in x or "61827" in x or "165319" in x or "168076" in x or "154043" in x or "73882" in x or "149404" in x or "152424" in x or "148937" in x or "157978" in x or "161056" in x or "80558" in x or "185859" in x or "167838" in x or "153919" in x or "43384" in x or "37061" in x]
#"Too Variable" Targetlist: HD101065, HD147084, HD50820, HD157978
print(SelectedListO17)



#USE THIS!!!:
cloudvelocitiesdffilename = PYTHONDIR + '/edibles/data/ISM_Cloud_Velocities_inKms_Approximate_Database_withtargetsasindex.dat'
cloudvelocitiesdf = pd.read_csv(cloudvelocitiesdffilename, delim_whitespace=True)
#pd.set_option('display.max_rows', None)
#print(cloudvelocitiesdf)
#print(cloudvelocitiesdf.index)
#print(cloudvelocitiesdf.index[0])
#print(cloudvelocitiesdf.cloud_velocity)
#print(cloudvelocitiesdf.comment)


#INSPECTING FOR-LOOP: That goes through all the files in the list by wavelength
#Unfortunately the EdiblesOracle tool can only select on Wavelength, so we can't group
#Different wavelength ranges (Different C6H Bands) after eachother.
#A new python piece in this code is needed for that!!!.
#LET OP: Verander de list hieronder per stukje wat je wilt bekijken!!!.
#"""
counter = 0
trial = 0
for filename in SelectedListO17: 
	counter = counter + 1
	print(counter)
	sp = EdiblesSpectrum(filename)
	print(sp.target)
	
	if sp.target not in cloudvelocitiesdf.index:
		cloudvelocityoftarget = 0
		print("De sp.target naam komt niet overeen met je naam in de Cloud Velocity Database (bij verschillende sp.target names per 1 Sightline), hier kan je niks mee, is nu in Barycentric frame")
	else:
		cloudvelocityoftarget = cloudvelocitiesdf.cloud_velocity[sp.target]
		print("Cloud Velocity of {} (in km/s) =".format(sp.target), cloudvelocityoftarget)
	
	wave = sp.wave
	flux = np.clip(sp.flux, 0, None) 
	bool_keep = (wave > plotrange[0]) & (wave < plotrange[1])
	plotwave = wave[bool_keep]
	plotflux = flux[bool_keep]
	normflux = plotflux / np.median(plotflux)
	#plt.figure()
	#C6H Bands1,2:
	plt.xlim(5255,5270)
	#C6H Bands 3,4,5:
	#plt.xlim(5190,5250)
	#C6H Band19:
	#plt.xlim(4740, 4750)
	#C8H:
	#plt.xlim(6250, 6270)
	#C10H:
	#plt.xlim(7120, 7160)
	ylim = [(np.min(normflux) + trial), (np.max(normflux) + trial)]
	#plt.ylim(ylim)
	
	trialplotnormflux = normflux + trial
	#trial = trial + 0.02
	
	#Edibles Data Geocentric Subset Spectrum Plot:
	#plt.plot(plotwave, trialplotnormflux, alpha=1, label="Geocentric Frame {}".format(sp.target))

	#Edibles Data Barycentric Subset Spectrum Plot:
	#C6H Bands 1,2:
	subset = sp.getSpectrum(xmin=5255, xmax=5270)
	#C6H Bands 3,4,5:
	#subset = sp.getSpectrum(xmin=5190, xmax=5250)
	#C6H Band 19:
	#subset = sp.getSpectrum(xmin=plotrange[0], xmax=plotrange[1])
	#C8H:
	#subset = sp.getSpectrum(xmin=6250, xmax=6270)
	#C10H:
	#subset = sp.getSpectrum(xmin=7120, xmax=7160)
	bary_wave = subset.bary_wave
	subsetflux = subset.flux 
	normsubsetflux = subsetflux / np.median(subsetflux)
	
	trialplotnormsubsetflux = normsubsetflux + trial
	trial = trial + 0.01
	
	#bool_keep2 = (bary_wave > plotrange[0]) & (bary_wave < plotrange[1])
	#plotbary_wave = bary_wave[bool_keep2]
	###plt.plot(bary_wave, trialplotnormsubsetflux, alpha=1, label="Barycentric Frame {}".format(sp.target))
	
	#Edibles Data ISM Cloud Rest Frame Spectrum Plot:
	ism_wave_formula = subset.bary_wave + (subset.bary_wave * (cloudvelocityoftarget/(3*10**5)))
	#print(ism_wave_formula)
	#bool_keep2 = (bary_wave > plotrange[0]) & (bary_wave < plotrange[1])
	#plotbary_wave = bary_wave[bool_keep2]
	plt.plot(ism_wave_formula, trialplotnormsubsetflux, alpha=1, label="ISM Frame {}".format(sp.target))
	
	plt.xlabel("Wavelength in Angstrom")
	plt.ylabel("Arbitrary Normalized Flux")
	plt.legend(framealpha=0.6, loc='upper center', ncol=8, fontsize= 'x-small')
	plt.grid()
	#plt.show()
	
#Labspectrum Plot:
# Rescale lab spectrum to plot range
dynrange = ylim[1]-ylim[0]
plt.plot(labspec.wavelength, labspec.norm * dynrange*2 + 1, color='k', alpha=1, label="C6H Bands 1,2")
	
#Plotting Everything Together:
plt.xlabel("Wavelength in Angstrom")
plt.ylabel("Arbitrary Normalized Flux")
#plt.ylim(0.96, 1.03)
plt.legend(framealpha=0.6, loc='lower center', ncol=8, fontsize= 'x-small')
plt.grid()
plt.show()
#"""



#CONTROLE FUNCTIES C6H BAND RANGES: Om het spectrum te kunnen plotten in een bepaalde O-Range 
#Gericht op C6H Vooralsnog, voor een nieuw molecuul zijn wss enigszins andere 
#Ranges nodig.
def spectrumO17(filename): 
	plotrange = [5250,5280]
	sp = EdiblesSpectrum(filename)
	print(sp.target+", Range O17, C6HBands 1,2 (Main Bands!!!)")
	
	cloudvelocityoftarget = cloudvelocitiesdf.cloud_velocity[sp.target]
	print("Cloud Velocity of {} (in km/s) =".format(sp.target), cloudvelocityoftarget)
	
	wave = sp.wave
	flux = np.clip(sp.flux, 0, None) 
	bool_keep = (wave > plotrange[0]) & (wave < plotrange[1])
	plotwave = wave[bool_keep]
	plotflux = flux[bool_keep]
	normflux = plotflux / np.median(plotflux)
	plt.figure()
	plt.xlim(5250,5280)
	ylim = [np.min(normflux), np.max(normflux)]
	plt.ylim(ylim)
	
	#Edibles Data Geocentric Subset Spectrum Plot:
	###plt.plot(plotwave, normflux, color='k', alpha=0.3, label="EDIBLES Data {}, GEOCENTRIC RESTFRAME".format(sp.target))

	#Edibles Data Barycentric Subset Spectrum Plot:
	subset = sp.getSpectrum(xmin=plotrange[0], xmax=plotrange[1])
	bary_wave = subset.bary_wave
	subsetflux = subset.flux 
	normsubsetflux = subsetflux / np.median(subsetflux)
	#bool_keep2 = (bary_wave > plotrange[0]) & (bary_wave < plotrange[1])
	#plotbary_wave = bary_wave[bool_keep2]
	plt.plot(bary_wave, normsubsetflux, color='b', alpha=0.5, label="EDIBLES Data {}, BARYCENTRIC RESTFRAME".format(sp.target))
	
	#Edibles Data ISM Cloud Rest Frame Spectrum Plot:
	ism_wave_formula = subset.bary_wave + (subset.bary_wave * (cloudvelocityoftarget/(3*10**5)))
	#print(ism_wave_formula)
	#bool_keep2 = (bary_wave > plotrange[0]) & (bary_wave < plotrange[1])
	#plotbary_wave = bary_wave[bool_keep2]
	plt.plot(ism_wave_formula, normsubsetflux, color='r', alpha=1, label="EDIBLES Data {}, ISM CLOUD REST FRAME".format(sp.target))
	
	#Labspectrum Plot:
	# Rescale lab spectrum to plot range
	dynrange = ylim[1]-ylim[0]
	plt.plot(labspec.wavelength, labspec.norm * dynrange/2 + 1, color='g', alpha=0.6, label="C6H Lab Spectrum Bands 1,2,3")
	
	#Plotting Everything Together:
	plt.xlabel("Wavelength in Angstrom")
	plt.ylabel("Arbitrary Normalized Flux")
	plt.legend()
	plt.grid()
	plt.show()

def spectrumO16(filename): 
	plotrange = [5190,5250]
	sp = EdiblesSpectrum(filename)
	print(sp.target+", Range O16, C6HBands 3,4,5")
	
	cloudvelocityoftarget = cloudvelocitiesdf.cloud_velocity[sp.target]
	#print("Cloud Velocity of {} (in km/s) =".format(sp.target), cloudvelocityoftarget)
	
	wave = sp.wave
	flux = np.clip(sp.flux, 0, None) 
	bool_keep = (wave > plotrange[0]) & (wave < plotrange[1])
	plotwave = wave[bool_keep]
	plotflux = flux[bool_keep]
	normflux = plotflux / np.median(plotflux)
	plt.figure()
	plt.xlim(5190,5250)
	ylim = [np.min(normflux), np.max(normflux)]
	plt.ylim(ylim)
	
	#Edibles Data Geocentric Subset Spectrum Plot:
	###plt.plot(plotwave, normflux, color='k', alpha=0.3, label="EDIBLES Data {}, GEOCENTRIC RESTFRAME".format(sp.target))

	#Edibles Data Barycentric Subset Spectrum Plot:
	subset = sp.getSpectrum(xmin=plotrange[0], xmax=plotrange[1])
	bary_wave = subset.bary_wave
	subsetflux = subset.flux 
	normsubsetflux = subsetflux / np.median(subsetflux)
	#bool_keep2 = (bary_wave > plotrange[0]) & (bary_wave < plotrange[1])
	#plotbary_wave = bary_wave[bool_keep2]
	plt.plot(bary_wave, normsubsetflux, color='b', alpha=0.5, label="EDIBLES Data {}, BARYCENTRIC RESTFRAME".format(sp.target))
	
	#Edibles Data ISM Cloud Rest Frame Spectrum Plot:
	ism_wave_formula = subset.bary_wave + (subset.bary_wave * (cloudvelocityoftarget/(3*10**5)))
	#print(ism_wave_formula)
	#bool_keep2 = (bary_wave > plotrange[0]) & (bary_wave < plotrange[1])
	#plotbary_wave = bary_wave[bool_keep2]
	plt.plot(ism_wave_formula, normsubsetflux, color='r', alpha=1, label="EDIBLES Data {}, ISM CLOUD REST FRAME".format(sp.target))
	
	#Labspectrum Plot:
	# Rescale lab spectrum to plot range
	dynrange = ylim[1]-ylim[0]
	plt.plot(labspec.wavelength, labspec.norm * dynrange/4 + 1, color='g', alpha=0.6, label="C6H Lab Spectrum Bands 3,4,5")
	
	#Stellar Lines: 
	plt.axvline(5243, color='m', label="Fe-III Stellar Line")
	plt.axvline(5208.5, color='y', label="O-II Stellar Line")
	plt.axvline(5190.5, color='b', label="O-II Stellar Line")
	
	#Plotting Everything Together:
	plt.xlabel("Wavelength in Angstrom")
	plt.ylabel("Arbitrary Normalized Flux")
	plt.legend()
	plt.grid()
	plt.show()

def spectrumO15(filename): 
	plotrange = [5155,5195]
	sp = EdiblesSpectrum(filename)
	print(sp.target+", Range O15, C6HBands 6,7,8")
	
	cloudvelocityoftarget = cloudvelocitiesdf.cloud_velocity[sp.target]
	#print("Cloud Velocity of {} (in km/s) =".format(sp.target), cloudvelocityoftarget)
	
	wave = sp.wave
	flux = np.clip(sp.flux, 0, None) 
	bool_keep = (wave > plotrange[0]) & (wave < plotrange[1])
	plotwave = wave[bool_keep]
	plotflux = flux[bool_keep]
	normflux = plotflux / np.median(plotflux)
	plt.figure()
	plt.xlim(5155,5195)
	ylim = [np.min(normflux), np.max(normflux)]
	plt.ylim(ylim)
	
	#Edibles Data Geocentric Subset Spectrum Plot:
	###plt.plot(plotwave, normflux, color='k', alpha=0.3, label="EDIBLES Data {}, GEOCENTRIC RESTFRAME".format(sp.target))
	
	#Edibles Data Barycentric Subset Spectrum Plot:
	subset = sp.getSpectrum(xmin=plotrange[0], xmax=plotrange[1])
	bary_wave = subset.bary_wave
	subsetflux = subset.flux 
	normsubsetflux = subsetflux / np.median(subsetflux)
	#bool_keep2 = (bary_wave > plotrange[0]) & (bary_wave < plotrange[1])
	#plotbary_wave = bary_wave[bool_keep2]
	plt.plot(bary_wave, normsubsetflux, color='b', alpha=0.5, label="EDIBLES Data {}, BARYCENTRIC RESTFRAME".format(sp.target))
	
	#Edibles Data ISM Cloud Rest Frame Spectrum Plot:
	ism_wave_formula = subset.bary_wave + (subset.bary_wave * (cloudvelocityoftarget/(3*10**5)))
	#print(ism_wave_formula)
	#bool_keep2 = (bary_wave > plotrange[0]) & (bary_wave < plotrange[1])
	#plotbary_wave = bary_wave[bool_keep2]
	plt.plot(ism_wave_formula, normsubsetflux, color='r', alpha=1, label="EDIBLES Data {}, ISM CLOUD REST FRAME".format(sp.target))
	
	#Labspectrum Plot:
	# Rescale lab spectrum to plot range
	dynrange = ylim[1]-ylim[0]
	plt.plot(labspec.wavelength, labspec.norm * dynrange/4 + 1, color='g', alpha=0.6, label="C6H Lab Spectrum Bands 6,7,8")
	
	#Stellar Lines:
	plt.axvline(5176.0, color='y', label="O-II Stellar Line")
	
	#Plotting Everything Together:
	plt.xlabel("Wavelength in Angstrom")
	plt.ylabel("Arbitrary Normalized Flux")
	plt.legend()
	plt.grid()
	plt.show()

def spectrumO14(filename): 
	plotrange = [5095,5155]
	sp = EdiblesSpectrum(filename)
	print(sp.target+", Range O14, C6H Bands 9,10,11,12,13,14")
	
	cloudvelocityoftarget = cloudvelocitiesdf.cloud_velocity[sp.target]
	#print("Cloud Velocity of {} (in km/s) =".format(sp.target), cloudvelocityoftarget)
	
	wave = sp.wave
	flux = np.clip(sp.flux, 0, None) 
	bool_keep = (wave > plotrange[0]) & (wave < plotrange[1])
	plotwave = wave[bool_keep]
	plotflux = flux[bool_keep]
	normflux = plotflux / np.median(plotflux)
	plt.figure()
	plt.xlim(5095,5155)
	ylim = [np.min(normflux), np.max(normflux)]
	plt.ylim(ylim)
	
	#Edibles Data Geocentric Subset Spectrum Plot:
	###plt.plot(plotwave, normflux, color='k', alpha=0.3, label="EDIBLES Data {}, GEOCENTRIC RESTFRAME".format(sp.target))
	
	#Edibles Data Barycentric Subset Spectrum Plot:
	subset = sp.getSpectrum(xmin=plotrange[0], xmax=plotrange[1])
	bary_wave = subset.bary_wave
	subsetflux = subset.flux 
	normsubsetflux = subsetflux / np.median(subsetflux)
	#bool_keep2 = (bary_wave > plotrange[0]) & (bary_wave < plotrange[1])
	#plotbary_wave = bary_wave[bool_keep2]
	plt.plot(bary_wave, normsubsetflux, color='b', alpha=0.5, label="EDIBLES Data {}, BARYCENTRIC RESTFRAME".format(sp.target))
	
	#Edibles Data ISM Cloud Rest Frame Spectrum Plot:
	ism_wave_formula = subset.bary_wave + (subset.bary_wave * (cloudvelocityoftarget/(3*10**5)))
	#print(ism_wave_formula)
	#bool_keep2 = (bary_wave > plotrange[0]) & (bary_wave < plotrange[1])
	#plotbary_wave = bary_wave[bool_keep2]
	plt.plot(ism_wave_formula, normsubsetflux, color='r', alpha=1, label="EDIBLES Data {}, ISM CLOUD REST FRAME".format(sp.target))
	
	#Labspectrum Plot:
	# Rescale lab spectrum to plot range
	dynrange = ylim[1]-ylim[0]
	plt.plot(labspec.wavelength, labspec.norm * dynrange/4 + 1, color='g', alpha=0.6, label="C6H Lab Spectrum Bands 9,10,11,12,13,14")
	
	#Stellar Lines:
	plt.axvline(5145.0, color='m', label="C-II Stellar Line")
	plt.axvline(5144.0, color='y', label="C-II Stellar Line")
	plt.axvline(5139.0, color='b', label="Weak C-II Stellar Line")
	plt.axvline(5133.0, color='m', label="C-II Stellar Line")
	plt.axvline(5127.0, color='y', label="Fe-III Stellar Line")
	
	#Plotting Everything Together:
	plt.xlabel("Wavelength in Angstrom")
	plt.ylabel("Arbitrary Normalized Flux")
	plt.legend()
	plt.grid()
	plt.show()

def spectrumO13(filename): 
	plotrange = [5050,5120]
	sp = EdiblesSpectrum(filename)
	print(sp.target+", Range O13, C6HBands 15, 16")
	
	cloudvelocityoftarget = cloudvelocitiesdf.cloud_velocity[sp.target]
	#print("Cloud Velocity of {} (in km/s) =".format(sp.target), cloudvelocityoftarget)
	
	wave = sp.wave
	flux = np.clip(sp.flux, 0, None) 
	bool_keep = (wave > plotrange[0]) & (wave < plotrange[1])
	plotwave = wave[bool_keep]
	plotflux = flux[bool_keep]
	normflux = plotflux / np.median(plotflux)
	plt.figure()
	plt.xlim(5050,5120)
	ylim = [np.min(normflux), np.max(normflux)]
	plt.ylim(ylim)
	
	#Edibles Data Geocentric Subset Spectrum Plot:
	###plt.plot(plotwave, normflux, color='k', alpha=0.3, label="EDIBLES Data {}, GEOCENTRIC RESTFRAME".format(sp.target))
	
	#Edibles Data Barycentric Subset Spectrum Plot:
	subset = sp.getSpectrum(xmin=plotrange[0], xmax=plotrange[1])
	bary_wave = subset.bary_wave
	subsetflux = subset.flux 
	normsubsetflux = subsetflux / np.median(subsetflux)
	#bool_keep2 = (bary_wave > plotrange[0]) & (bary_wave < plotrange[1])
	#plotbary_wave = bary_wave[bool_keep2]
	plt.plot(bary_wave, normsubsetflux, color='b', alpha=0.5, label="EDIBLES Data {}, BARYCENTRIC RESTFRAME".format(sp.target))
	
	#Edibles Data ISM Cloud Rest Frame Spectrum Plot:
	ism_wave_formula = subset.bary_wave + (subset.bary_wave * (cloudvelocityoftarget/(3*10**5)))
	#print(ism_wave_formula)
	#bool_keep2 = (bary_wave > plotrange[0]) & (bary_wave < plotrange[1])
	#plotbary_wave = bary_wave[bool_keep2]
	plt.plot(ism_wave_formula, normsubsetflux, color='r', alpha=1, label="EDIBLES Data {}, ISM CLOUD REST FRAME".format(sp.target))
	
	#Labspectrum Plot:
	# Rescale lab spectrum to plot range
	dynrange = ylim[1]-ylim[0]
	plt.plot(labspec.wavelength, labspec.norm * dynrange/4 + 1, color='g', alpha=0.6, label="C6H Lab Spectrum Bands 15,16")
	
	#Plotting Everything Together:
	plt.xlabel("Wavelength in Angstrom")
	plt.ylabel("Arbitrary Normalized Flux")
	plt.legend()
	plt.grid()
	plt.show()

def spectrumO12(filename): 
	plotrange = [5010,5050]
	sp = EdiblesSpectrum(filename)
	print(sp.target+", Range O12, C6HBand 17")
	
	cloudvelocityoftarget = cloudvelocitiesdf.cloud_velocity[sp.target]
	#print("Cloud Velocity of {} (in km/s) =".format(sp.target), cloudvelocityoftarget)
	
	wave = sp.wave
	flux = np.clip(sp.flux, 0, None) 
	bool_keep = (wave > plotrange[0]) & (wave < plotrange[1])
	plotwave = wave[bool_keep]
	plotflux = flux[bool_keep]
	normflux = plotflux / np.median(plotflux)
	plt.figure()
	plt.xlim(5010,5050)
	ylim = [np.min(normflux), np.max(normflux)]
	plt.ylim(ylim)
	
	#Edibles Data Geocentric Subset Spectrum Plot:
	###plt.plot(plotwave, normflux, color='k', alpha=0.3, label="EDIBLES Data {}, GEOCENTRIC RESTFRAME".format(sp.target))

	#Edibles Data Barycentric Subset Spectrum Plot:
	subset = sp.getSpectrum(xmin=plotrange[0], xmax=plotrange[1])
	bary_wave = subset.bary_wave
	subsetflux = subset.flux 
	normsubsetflux = subsetflux / np.median(subsetflux)
	#bool_keep2 = (bary_wave > plotrange[0]) & (bary_wave < plotrange[1])
	#plotbary_wave = bary_wave[bool_keep2]
	plt.plot(bary_wave, normsubsetflux, color='b', alpha=0.5, label="EDIBLES Data {}, BARYCENTRIC RESTFRAME".format(sp.target))
	
	#Edibles Data ISM Cloud Rest Frame Spectrum Plot:
	ism_wave_formula = subset.bary_wave + (subset.bary_wave * (cloudvelocityoftarget/(3*10**5)))
	#print(ism_wave_formula)
	#bool_keep2 = (bary_wave > plotrange[0]) & (bary_wave < plotrange[1])
	#plotbary_wave = bary_wave[bool_keep2]
	plt.plot(ism_wave_formula, normsubsetflux, color='r', alpha=1, label="EDIBLES Data {}, ISM CLOUD REST FRAME".format(sp.target))
	
	#Labspectrum Plot:
	# Rescale lab spectrum to plot range
	dynrange = ylim[1]-ylim[0]
	plt.plot(labspec.wavelength, labspec.norm * dynrange/4 + 1, color='g', alpha=0.6, label="C6H Lab Spectrum Band 17")
	
	#Stellar Lines:
	plt.axvline(5016.0, color='m', label="He-I Stellar Line")
	
	#Plotting Everything Together:
	plt.xlabel("Wavelength in Angstrom")
	plt.ylabel("Arbitrary Normalized Flux")
	plt.legend()
	plt.grid()
	plt.show()

def spectrumO7(filename): 
	plotrange = [4810,4870]
	sp = EdiblesSpectrum(filename)
	print(sp.target+", Range O7, C6HBand 18")
	
	cloudvelocityoftarget = cloudvelocitiesdf.cloud_velocity[sp.target]
	#print("Cloud Velocity of {} (in km/s) =".format(sp.target), cloudvelocityoftarget)
	
	wave = sp.wave
	flux = np.clip(sp.flux, 0, None) 
	bool_keep = (wave > plotrange[0]) & (wave < plotrange[1])
	plotwave = wave[bool_keep]
	plotflux = flux[bool_keep]
	normflux = plotflux / np.median(plotflux)
	plt.figure()
	plt.xlim(4810,4870)
	ylim = [np.min(normflux), np.max(normflux)]
	plt.ylim(ylim)
	
	#Edibles Data Geocentric Subset Spectrum Plot:
	###plt.plot(plotwave, normflux, color='k', alpha=0.3, label="EDIBLES Data {}, GEOCENTRIC RESTFRAME".format(sp.target))

	#Edibles Data Barycentric Subset Spectrum Plot:
	subset = sp.getSpectrum(xmin=plotrange[0], xmax=plotrange[1])
	bary_wave = subset.bary_wave
	subsetflux = subset.flux 
	normsubsetflux = subsetflux / np.median(subsetflux)
	#bool_keep2 = (bary_wave > plotrange[0]) & (bary_wave < plotrange[1])
	#plotbary_wave = bary_wave[bool_keep2]
	plt.plot(bary_wave, normsubsetflux, color='b', alpha=0.5, label="EDIBLES Data {}, BARYCENTRIC RESTFRAME".format(sp.target))
	
	#Edibles Data ISM Cloud Rest Frame Spectrum Plot:
	ism_wave_formula = subset.bary_wave + (subset.bary_wave * (cloudvelocityoftarget/(3*10**5)))
	#print(ism_wave_formula)
	#bool_keep2 = (bary_wave > plotrange[0]) & (bary_wave < plotrange[1])
	#plotbary_wave = bary_wave[bool_keep2]
	plt.plot(ism_wave_formula, normsubsetflux, color='r', alpha=1, label="EDIBLES Data {}, ISM CLOUD REST FRAME".format(sp.target))
	
	#Labspectrum Plot:
	# Rescale lab spectrum to plot range
	dynrange = ylim[1]-ylim[0]
	plt.plot(labspec.wavelength, labspec.norm * dynrange/4 + 1, color='g', alpha=0.6, label="C6H Lab Spectrum Band 18")
	
	#Stellar Lines:
	plt.axvline(4861.0, color='m', label="O-II Stellar Line")
	plt.axvline(4842.5, color='y', label="Weak O-II Stellar Line")
	
	#Plotting Everything Together:
	plt.xlabel("Wavelength in Angstrom")
	plt.ylabel("Arbitrary Normalized Flux")
	plt.legend()
	plt.grid()
	plt.show()

def spectrumO4(filename): 
	plotrange = [4720,4760]
	sp = EdiblesSpectrum(filename)
	print(sp.target+", Range O4, C6HBand 19")
	
	cloudvelocityoftarget = cloudvelocitiesdf.cloud_velocity[sp.target]
	#print("Cloud Velocity of {} (in km/s) =".format(sp.target), cloudvelocityoftarget)
	
	wave = sp.wave
	flux = np.clip(sp.flux, 0, None) 
	bool_keep = (wave > plotrange[0]) & (wave < plotrange[1])
	plotwave = wave[bool_keep]
	plotflux = flux[bool_keep]
	normflux = plotflux / np.median(plotflux)
	plt.figure()
	plt.xlim(4720,4760)
	ylim = [np.min(normflux), np.max(normflux)]
	plt.ylim(ylim)
	
	#Edibles Data Geocentric Subset Spectrum Plot:
	###plt.plot(plotwave, normflux, color='k', alpha=0.3, label="EDIBLES Data {}, GEOCENTRIC RESTFRAME".format(sp.target))
	
	#Edibles Data Barycentric Subset Spectrum Plot:
	subset = sp.getSpectrum(xmin=plotrange[0], xmax=plotrange[1])
	bary_wave = subset.bary_wave
	subsetflux = subset.flux 
	normsubsetflux = subsetflux / np.median(subsetflux)
	#bool_keep2 = (bary_wave > plotrange[0]) & (bary_wave < plotrange[1])
	#plotbary_wave = bary_wave[bool_keep2]
	plt.plot(bary_wave, normsubsetflux, color='b', alpha=0.5, label="EDIBLES Data {}, BARYCENTRIC RESTFRAME".format(sp.target))
	
	#Edibles Data ISM Cloud Rest Frame Spectrum Plot:
	ism_wave_formula = subset.bary_wave + (subset.bary_wave * (cloudvelocityoftarget/(3*10**5)))
	#print(ism_wave_formula)
	#bool_keep2 = (bary_wave > plotrange[0]) & (bary_wave < plotrange[1])
	#plotbary_wave = bary_wave[bool_keep2]
	plt.plot(ism_wave_formula, normsubsetflux, color='r', alpha=1, label="EDIBLES Data {}, ISM CLOUD REST FRAME".format(sp.target))
	
	#Labspectrum Plot:
	# Rescale lab spectrum to plot range
	dynrange = ylim[1]-ylim[0]
	plt.plot(labspec.wavelength, labspec.norm * dynrange/4 + 1, color='g', alpha=0.6, label="C6H Lab Spectrum Band 19")
	
	#Stellar Line:
	plt.axvline(4751.0, color='m', label="O-II Stellar Line")
	
	#Plotting Everything Together:
	plt.xlabel("Wavelength in Angstrom")
	plt.ylabel("Arbitrary Normalized Flux")
	plt.legend()
	plt.grid()
	plt.show()




#TARGET LISTS: Deze worden bewaard in Notion, hier zorgt het alleen maar voor clutter
#1e Target List (voor de presentatie): Deze hebben we deze Lijst puur gebaseerd op Bands 1,2
#2e Target List: Deze is vnml gebaseerd op bands 3,4,5 en nog ge cross-referenced op Bands 1,2
#3e BEST TARGET LIST: Based on bovenstaande 2 lijsten en het GEHELE SPECTRUM. Ik durf te
#wedden dat we geen betere target list gaan krijgen voor C6H dan deze, dus hier wil ik 
#het eigenlijk bij laten qua C6H target list en eventueel verder gaan analyseren dieper
#induikend op deze BEST TARGET LIST.

#BANDWIDTH RANGES C6H:
#O17: 5250-5280: C6H Bands 1,2
#O16: 5190-5250: C6H Bands 3,4,5
#O15: 5155-5195: C6H Bands 6,7,8
#O14: 5095-5155: C6H Bands 9,10,11,12,13,14
#O13: 5060-5120: C6H Bands 15,16
#O12: 5010-5050: C6H Band 17
#O7: 4810-4870: C6H Band 18
#O4: 4720-4760: C6H Band 19



#FULL SPECTRUM FUNCTION C6H:
def fullspectrum(shortenedfilename):
	#print(shortenedfilename)
	O17name = "{}_O17.fits".format(shortenedfilename)
	#print(O17name)
	O16name = "{}_O16.fits".format(shortenedfilename)
	O15name = "{}_O15.fits".format(shortenedfilename)
	O14name = "{}_O14.fits".format(shortenedfilename)
	O13name = "{}_O13.fits".format(shortenedfilename)
	O12name = "{}_O12.fits".format(shortenedfilename)
	O7name = "{}_O7.fits".format(shortenedfilename)
	O4name = "{}_O4.fits".format(shortenedfilename)
	spectrumO17(O17name)
	spectrumO16(O16name)
	spectrumO15(O15name)
	spectrumO14(O14name)
	spectrumO13(O13name)
	spectrumO12(O12name)
	spectrumO7(O7name)
	spectrumO4(O4name)



#FULLSPECTRUM C6H FUNCTIE INGEVULD VOOR ALLE TARGETS: Let er dus op dat je het zo invuld als
#HD37367 hieronder, stoppende bij de datum, om het correct te krijgen!!!:
#NIEUWE BEST TARGET ORDER, BASED ON FULL SPECTRUM:
#"""
#MAIN BEST LIST, UPDATED NA BAND19:
fullspectrum("/HD63804/RED_564/HD63804_w564_redl_20190203")
#multi-plotted
fullspectrum("/HD37367/RED_564/HD37367_w564_redl_20141011")
fullspectrum("/HD81188/RED_564/HD81188_w564_redl_20170418")
fullspectrum("/HD164353/RED_564/HD164353_w564_redl_20170615")
fullspectrum("/HD164073/RED_564/HD164073_w564_redl_20180830")
fullspectrum("/HD63804/RED_564/HD63804_w564_redl_20190203")
fullspectrum("/HD80558/RED_564/HD80558_w564_redl_20150603")
fullspectrum("/HD170740/RED_564/HD170740_w564_redl_20140916")
#non multi-plotted
fullspectrum("/HD79186/RED_564/HD79186_w564_redl_20141013")
fullspectrum("/HD75860/RED_564/HD75860_w564_redl_20170429")
fullspectrum("/HD167838/RED_564/HD167838_w564_redl_20160808")

#BACKUP BEST LIST:
fullspectrum("/HD63804/RED_564/HD63804_w564_redl_20190203")
fullspectrum("/HD164073/RED_564/HD164073_w564_redl_20180830") 
fullspectrum("/HD148184/RED_564/HD148184_w564_redl_20180623")
fullspectrum("/HD24398/RED_564/HD24398_w564_redl_20140923")
fullspectrum("/HD154043/RED_564/HD154043_w564_redl_20170502")
fullspectrum("/HD144470/RED_564/HD144470_w564_redl_20150720")
fullspectrum("/HD158926/RED_564/lambdaSco_w564_redl_20150720")
fullspectrum("/HD114886/RED_564/HD114886_w564_redl_20170503")
fullspectrum("/HD61827/RED_564/HD61827_w564_redl_20180324")
fullspectrum("/HD101065/RED_564/HD101065_w564_redl_20170420")
fullspectrum("/HD157978/RED_564/HD157978_w564_redl_20160508")
#"""

#TOEKOMST: We kunnen weer alle targets afgaan op deze manier, maar wss gaan we daar
#niks boeiends meer uit vinden, dus THIS IS IT FOR C6H!!!.
#We hebben hierboven een intiele lijst gecreerd based on bands 3,4,5 mainly 
#Gecombineerd met Bands 1,2 en daaruit een main en backup BEST LIST uit 
#gevormd Based on het Volledige Spectrum!!!!!.
