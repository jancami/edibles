from edibles.edibles import DATADIR
from edibles.edibles import PYTHONDIR
from edibles.edibles.utils.edibles_spectrum import EdiblesSpectrum
from edibles.edibles.utils.edibles_oracle import EdiblesOracle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from PyAstronomy import pyasl

#Main C6H 19Bands Simulated Data 20K:
labfilename = PYTHONDIR + '/edibles/data/Labdata/C6H/C6H_Data_All19Bands_20K_Resolvingpower5timeslesstoGau03.dat'

#Temperature Changes for C6H Spectrum:
labfilename1 = PYTHONDIR + '/edibles/data/Labdata/C6H/C6H_Data_All19Bands_2K.dat'
labfilename2 = PYTHONDIR + '/edibles/data/Labdata/C6H/C6H_Data_All19Bands_5K.dat'
labfilename3 = PYTHONDIR + '/edibles/data/Labdata/C6H/C6H_Data_All19Bands_10K.dat'
labfilename4 = PYTHONDIR + '/edibles/data/Labdata/C6H/C6H_Data_All19Bands_15K.dat'
labfilename5 = PYTHONDIR + '/edibles/data/Labdata/C6H/C6H_Data_All19Bands_20K.dat'
labfilename6 = PYTHONDIR + '/edibles/data/Labdata/C6H/C6H_Data_All19Bands_25K.dat'
labfilename7 = PYTHONDIR + '/edibles/data/Labdata/C6H/C6H_Data_All19Bands_50K.dat'
labfilename8 = PYTHONDIR + '/edibles/data/Labdata/C6H/C6H_Data_All19Bands_200K.dat'
labfilename9 = PYTHONDIR + '/edibles/data/Labdata/C6H/C6H_Data_All19Bands_2000K.dat'

#C8H:
#labfilename = PYTHONDIR + '/edibles/data/Labdata/C8H/C8H_Bands1_2_15K.dat'
#C10H:
#labfilename = PYTHONDIR + '/edibles/data/Labdata/C10H/C10H_Band1_15K.dat'

#LABSPECTRUM/SIMULATED DATA (in wavenumbers and flux, so needs to be converted to air wavelength in angstrom)
#Main:
labspec = pd.read_csv(labfilename, delim_whitespace=True)
labspec['wavelength'] = pyasl.vactoair2(1e8/labspec['wno'],mode='ciddor')
normint = labspec.int - np.median(labspec.int)
labspec['norm'] = normint / np.max(normint) * -1
#1: 
labspec1 = pd.read_csv(labfilename1, delim_whitespace=True)
labspec1['wavelength'] = pyasl.vactoair2(1e8/labspec1['wno'],mode='ciddor')
normint1 = labspec1.int - np.median(labspec1.int)
labspec1['norm'] = normint1 / np.max(normint1) * -1
#2: 
labspec2 = pd.read_csv(labfilename2, delim_whitespace=True)
labspec2['wavelength'] = pyasl.vactoair2(1e8/labspec2['wno'],mode='ciddor')
normint2 = labspec2.int - np.median(labspec2.int)
labspec2['norm'] = normint2 / np.max(normint2) * -1
#3: 
labspec3 = pd.read_csv(labfilename3, delim_whitespace=True)
labspec3['wavelength'] = pyasl.vactoair2(1e8/labspec3['wno'],mode='ciddor')
normint3 = labspec3.int - np.median(labspec3.int)
labspec3['norm'] = normint3 / np.max(normint3) * -1
#4: 
labspec4 = pd.read_csv(labfilename4, delim_whitespace=True)
labspec4['wavelength'] = pyasl.vactoair2(1e8/labspec4['wno'],mode='ciddor')
normint4 = labspec4.int - np.median(labspec4.int)
labspec4['norm'] = normint4 / np.max(normint4) * -1
#5: 
labspec5 = pd.read_csv(labfilename5, delim_whitespace=True)
labspec5['wavelength'] = pyasl.vactoair2(1e8/labspec5['wno'],mode='ciddor')
normint5 = labspec5.int - np.median(labspec5.int)
labspec5['norm'] = normint5 / np.max(normint5) * -1
#6: 
labspec6 = pd.read_csv(labfilename6, delim_whitespace=True)
labspec6['wavelength'] = pyasl.vactoair2(1e8/labspec6['wno'],mode='ciddor')
normint6 = labspec6.int - np.median(labspec6.int)
labspec6['norm'] = normint6 / np.max(normint6) * -1
#7: 
labspec7 = pd.read_csv(labfilename7, delim_whitespace=True)
labspec7['wavelength'] = pyasl.vactoair2(1e8/labspec7['wno'],mode='ciddor')
normint7 = labspec7.int - np.median(labspec7.int)
labspec7['norm'] = normint7 / np.max(normint7) * -1
#8: 
labspec8 = pd.read_csv(labfilename8, delim_whitespace=True)
labspec8['wavelength'] = pyasl.vactoair2(1e8/labspec8['wno'],mode='ciddor')
normint8 = labspec8.int - np.median(labspec8.int)
labspec8['norm'] = normint8 / np.max(normint8) * -1
#9: 
labspec9 = pd.read_csv(labfilename9, delim_whitespace=True)
labspec9['wavelength'] = pyasl.vactoair2(1e8/labspec9['wno'],mode='ciddor')
normint9 = labspec9.int - np.median(labspec9.int)
labspec9['norm'] = normint9 / np.max(normint9) * -1


#STACKING:
#
#
#
#
#
#
#
#
#C8H LISTS:
"""
#C8H Main Bands List: 
plotrange = [6258, 6265]
#plotrange = [6261, 6264]
pythia = EdiblesOracle()
ListC8H = pythia.GetObsListByWavelength(6262)
ListC8HO10 = [x for x in ListC8H if not "O11" in x]
print(ListC8HO10)

SelectedListC8H = [x for x in ListC8HO10 if "101065" not in x and "147084" not in x and "50820" not in x and "157978" not in x and "133518" not in x]
High_Extinction_SelectedListC8H = [x for x in SelectedListC8H if "183143" in x or "169454" in x or "167971" in x or "170938" in x or "186841" in x or "164740" in x or "186745" in x or "75860" in x or "156201" in x or "112272" in x or "168076" in x or "154043" in x or "73882" in x or "149404" in x or "152424" in x or "148937" in x or "161056" in x or "80558" in x or "185859" in x or "167838" in x or "153919" in x or "43384" in x or "37061" in x]
"""

#C10H LISTS, WERKT NIET!!!:
"""
#C10H Main Bands List: 
plotrange = [7135, 7145]
pythia = EdiblesOracle()
ListC10H = pythia.GetObsListByWavelength(7130)
ListC10HO7 = [x for x in ListC10H if not "O6" in x]
print(ListC10HO7)

SelectedListC10H = [x for x in ListC10HO7 if "101065" not in x and "147084" not in x and "50820" not in x and "157978" not in x and "133518" not in x]
High_Extinction_SelectedListC10H = [x for x in SelectedListC10H if "183143" in x or "169454" in x or "167971" in x or "170938" in x or "186841" in x or "164740" in x or "186745" in x or "75860" in x or "156201" in x or "112272" in x or "168076" in x or "154043" in x or "73882" in x or "149404" in x or "152424" in x or "148937" in x or "161056" in x or "80558" in x or "185859" in x or "167838" in x or "153919" in x or "43384" in x or "37061" in x]
"""

#C6H LISTS FOR BANDS 3,4,5 (STILL IRRELEVANT DUE TO STELLAR LINES):
"""
plotrange = [5199, 5220]
pythia = EdiblesOracle()
List345 = pythia.GetObsListByWavelength(5220)
List345_O16 = [x for x in List345 if not "O17" in x]
print(List345_O16)

SelectedListBands345 = [x for x in List345_O16 if "101065" not in x and "147084" not in x and "50820" not in x and "157978" not in x and "133518" not in x]

#(not relevant)SELECTED, EXTINCTION ORDERED (Increasing E(B-V)) TARGETS WITH STELLAR FEATURES
ControlledTargetListStellarBands345 = ["/HD43384/RED_564/HD43384_w564_redl_20150930_O16.fits", "/HD167838/RED_564/HD167838_w564_redl_20180902_O16.fits", "/HD167838/RED_564/HD167838_w564_redl_20160808_O16.fits", "/HD80558/RED_564/HD80558_w564_redl_20160221_O16.fits", "/HD80558/RED_564/HD80558_w564_redl_20150603_O16.fits", "/HD61827/RED_564/HD61827_w564_redl_20180324_O16.fits", "/HD186745/RED_564/HD186745_w564_redl_20150920_O16.fits", "/HD186745/RED_564/HD186745_w564_redl_20160909_O16.fits", "/HD183143/RED_564/HD183143_w564_redl_20180903_O16.fits", "/HD183143/RED_564/HD183143_w564_redl_20180907_O16.fits"]
print(len(ControlledTargetListStellarBands345))

#(not relevant)SELECTED, EXTINCTION ORDERED (Increasing E(B-V)) TARGETS WITHOUT STELLAR FEATURES
ControlledTargetListNoStellarBands345 = ["/HD37061/RED_564/HD37061_w564_redl_20190101_O16.fits", "/HD37061/RED_564/HD37061_w564_redl_20160912_O16.fits", "/HD153919/RED_564/HD153919_w564_redl_20160707_O16.fits", "/HD185859/RED_564/HD185859_w564_redl_20160814_O16.fits", "/HD185859/RED_564/HD185859_w564_redl_20150920_O16.fits", "/HD185859/RED_564/HD185859_w564_redl_20160813_O16.fits", "/HD148937/RED_564/HD148937_w564_redl_20170424_O16.fits", "/HD148937/RED_564/HD148937_w564_redl_20150817_O16.fits", "/HD152424/RED_564/HD152424_w564_redl_20160411_O16.fits", "/HD149404/RED_564/HD149404_w564_redl_20180630_O16.fits", "/HD149404/RED_564/HD149404_w564_redl_20180623_O16.fits", "/HD149404/RED_564/HD149404_w564_redl_20170418_O16.fits", "/HD149404/RED_564/HD149404_w564_redl_20170507_O16.fits", "/HD154043/RED_564/HD154043_w564_redl_20170502_O16.fits", "/HD168076/RED_564/HD168076_w564_redl_20180912_O16.fits", "/HD168076/RED_564/HD168076_w564_redl_20180911_O16.fits", "/HD156201/RED_564/HD156201_w564_redl_20180830_O16.fits", "/HD75860/RED_564/HD75860_w564_redl_20170429_O16.fits", "/HD75860/RED_564/HD75860_w564_redl_20170430_O16.fits", "/HD186841/RED_564/HD186841_w564_redl_20160910_O16.fits", "/HD186841/RED_564/HD186841_w564_redl_20160909_O16.fits", "/HD112272/RED_564/HD112272_w564_redl_20170430_O16.fits", "/HD112272/RED_564/HD112272_w564_redl_20170501_O16.fits", "/HD170938/RED_564/HD170938_w564_redl_20180911_O16.fits", "/HD167971/RED_564/HD167971_w564_redl_20140921_O16.fits", "/HD169454/RED_564/HD169454_w564_redl_20160714_O16.fits", "/HD169454/RED_564/HD169454_w564_redl_20160724_O16.fits", "/HD147889/RED_564/HD147889_w564_redl_20140928_O16.fits"]
#ControlledTargetListNoStellar_O18 = ["/HD37061/RED_564/HD37061_w564_redl_20190101_O18.fits", "/HD37061/RED_564/HD37061_w564_redl_20160912_O18.fits", "/HD153919/RED_564/HD153919_w564_redl_20160707_O18.fits", "/HD185859/RED_564/HD185859_w564_redl_20160814_O18.fits", "/HD185859/RED_564/HD185859_w564_redl_20150920_O18.fits", "/HD185859/RED_564/HD185859_w564_redl_20160813_O18.fits", "/HD148937/RED_564/HD148937_w564_redl_20170424_O18.fits", "/HD148937/RED_564/HD148937_w564_redl_20150817_O18.fits", "/HD152424/RED_564/HD152424_w564_redl_20160411_O18.fits", "/HD149404/RED_564/HD149404_w564_redl_20180630_O18.fits", "/HD149404/RED_564/HD149404_w564_redl_20180623_O18.fits", "/HD149404/RED_564/HD149404_w564_redl_20170418_O18.fits", "/HD149404/RED_564/HD149404_w564_redl_20170507_O18.fits", "/HD154043/RED_564/HD154043_w564_redl_20170502_O18.fits", "/HD168076/RED_564/HD168076_w564_redl_20180912_O18.fits", "/HD168076/RED_564/HD168076_w564_redl_20180911_O18.fits", "/HD156201/RED_564/HD156201_w564_redl_20180830_O18.fits", "/HD75860/RED_564/HD75860_w564_redl_20170429_O18.fits", "/HD75860/RED_564/HD75860_w564_redl_20170430_O18.fits", "/HD186841/RED_564/HD186841_w564_redl_20160910_O18.fits", "/HD186841/RED_564/HD186841_w564_redl_20160909_O18.fits", "/HD112272/RED_564/HD112272_w564_redl_20170430_O18.fits", "/HD112272/RED_564/HD112272_w564_redl_20170501_O18.fits", "/HD170938/RED_564/HD170938_w564_redl_20180911_O18.fits", "/HD167971/RED_564/HD167971_w564_redl_20140921_O18.fits", "/HD169454/RED_564/HD169454_w564_redl_20160714_O18.fits", "/HD169454/RED_564/HD169454_w564_redl_20160724_O18.fits", "/HD147889/RED_564/HD147889_w564_redl_20140928_O18.fits"]
print(len(ControlledTargetListNoStellarBands345))
print(ControlledTargetListNoStellarBands345)

High_Extinction_SelectedListBands345 = [x for x in SelectedListBands345 if "183143" in x or "169454" in x or "167971" in x or "170938" in x or "186841" in x or "164740" in x or "186745" in x or "75860" in x or "156201" in x or "112272" in x or "168076" in x or "154043" in x or "73882" in x or "149404" in x or "152424" in x or "148937" in x or "161056" in x or "80558" in x or "185859" in x or "167838" in x or "153919" in x or "43384" in x or "37061" in x]
Low_Extinction_SelectedListBands345 = [x for x in SelectedListBands345 if "23016" in x or "157246" in x or "38771" in x or "148605" in x or "93030" in x or "66811" in x or "37128" in x or "158926" in x or "81188" in x]
Low_Extinction_SelectedListBands345_33AddedTargets = [x for x in SelectedListBands345 if "23016" in x or "157246" in x or "38771" in x or "148605" in x or "93030" in x or "66811" in x or "37128" in x or "158926" in x or "81188" in x or "40111" in x or "180554" in x or "57061" in x or "55879" in x or "36861" in x or "133518" in x or "34748" in x or "36822" in x or "75309" in x or "143275" in x or "49787" in x or "53975" in x or "94493" in x or "113904" in x or "135591" in x or "144470" in x or "37041" in x or "22951" in x or "104705" in x or "184915" in x or "167264" in x or "103779" in x or "91824" in x or "109399" in x or "93843" in x or "171957" in x or "145502" in x or "145502" in x or "54439" in x or "155806" in x]
"""

#C6H LISTS FOR BANDS6,7,8:
"""
plotrange = [5140, 5180]
#plotrange = [5165, 5178]
#plotrange = [5150, 5165]
pythia = EdiblesOracle()
List678 = pythia.GetObsListByWavelength(5171)
List678_O15 = [x for x in List678 if not "O16" in x]
List678_O15x = [x for x in List678_O15 if not "O14" in x]
print(List678_O15x)

SelectedListBands678 = [x for x in List678_O15x if "101065" not in x and "147084" not in x and "50820" not in x and "157978" not in x and "133518" not in x]

#Stellar Line Check:
#SELECTED, EXTINCTION ORDERED (Increasing E(B-V)) TARGETS WITHOUT STELLAR FEATURES
ControlledTargetListNoStellarBands678 = ["/HD37061/RED_564/HD37061_w564_redl_20190101_O15.fits", "/HD37061/RED_564/HD37061_w564_redl_20160912_O15.fits", "/HD153919/RED_564/HD153919_w564_redl_20160707_O15.fits", "/HD185859/RED_564/HD185859_w564_redl_20160814_O15.fits", "/HD185859/RED_564/HD185859_w564_redl_20150920_O15.fits", "/HD185859/RED_564/HD185859_w564_redl_20160813_O15.fits", "/HD148937/RED_564/HD148937_w564_redl_20170424_O15.fits", "/HD148937/RED_564/HD148937_w564_redl_20150817_O15.fits", "/HD152424/RED_564/HD152424_w564_redl_20160411_O15.fits", "/HD149404/RED_564/HD149404_w564_redl_20180630_O15.fits", "/HD149404/RED_564/HD149404_w564_redl_20180623_O15.fits", "/HD149404/RED_564/HD149404_w564_redl_20170418_O15.fits", "/HD149404/RED_564/HD149404_w564_redl_20170507_O15.fits", "/HD154043/RED_564/HD154043_w564_redl_20170502_O15.fits", "/HD168076/RED_564/HD168076_w564_redl_20180912_O15.fits", "/HD168076/RED_564/HD168076_w564_redl_20180911_O15.fits", "/HD156201/RED_564/HD156201_w564_redl_20180830_O15.fits", "/HD75860/RED_564/HD75860_w564_redl_20170429_O15.fits", "/HD75860/RED_564/HD75860_w564_redl_20170430_O15.fits", "/HD186841/RED_564/HD186841_w564_redl_20160910_O15.fits", "/HD186841/RED_564/HD186841_w564_redl_20160909_O15.fits", "/HD112272/RED_564/HD112272_w564_redl_20170430_O15.fits", "/HD112272/RED_564/HD112272_w564_redl_20170501_O15.fits", "/HD170938/RED_564/HD170938_w564_redl_20180911_O15.fits", "/HD167971/RED_564/HD167971_w564_redl_20140921_O15.fits", "/HD169454/RED_564/HD169454_w564_redl_20160714_O15.fits", "/HD169454/RED_564/HD169454_w564_redl_20160724_O15.fits", "/HD147889/RED_564/HD147889_w564_redl_20140928_O15.fits"]
print(len(ControlledTargetListNoStellarBands678))
print(ControlledTargetListNoStellarBands678)

High_Extinction_SelectedListBands678 = [x for x in SelectedListBands678 if "183143" in x or "169454" in x or "167971" in x or "170938" in x or "186841" in x or "164740" in x or "186745" in x or "75860" in x or "156201" in x or "112272" in x or "168076" in x or "154043" in x or "73882" in x or "149404" in x or "152424" in x or "148937" in x or "161056" in x or "80558" in x or "185859" in x or "167838" in x or "153919" in x or "43384" in x or "37061" in x]
Low_Extinction_SelectedListBands678 = [x for x in SelectedListBands678 if "23016" in x or "157246" in x or "38771" in x or "148605" in x or "93030" in x or "66811" in x or "37128" in x or "158926" in x or "81188" in x]
Low_Extinction_SelectedListBands678_33AddedTargets = [x for x in SelectedListBands678 if "23016" in x or "157246" in x or "38771" in x or "148605" in x or "93030" in x or "66811" in x or "37128" in x or "158926" in x or "81188" in x or "40111" in x or "180554" in x or "57061" in x or "55879" in x or "36861" in x or "133518" in x or "34748" in x or "36822" in x or "75309" in x or "143275" in x or "49787" in x or "53975" in x or "94493" in x or "113904" in x or "135591" in x or "144470" in x or "37041" in x or "22951" in x or "104705" in x or "184915" in x or "167264" in x or "103779" in x or "91824" in x or "109399" in x or "93843" in x or "171957" in x or "145502" in x or "145502" in x or "54439" in x or "155806" in x]
"""

#C6H LISTS FOR BANDS9,10,11,12,12,13,14:
"""
plotrange = [5095,5155]
pythia = EdiblesOracle()
List910 = pythia.GetObsListByWavelength(5100)
List910_O14 = [x for x in List910 if not "O15" in x]
List910_O14x = [x for x in List910_O14 if not "O13" in x]
print(List910_O14x)

ControlledTargetListNoStellarBands914 = ["/HD37061/RED_564/HD37061_w564_redl_20190101_O14.fits", "/HD37061/RED_564/HD37061_w564_redl_20160912_O14.fits", "/HD153919/RED_564/HD153919_w564_redl_20160707_O14.fits", "/HD185859/RED_564/HD185859_w564_redl_20160814_O14.fits", "/HD185859/RED_564/HD185859_w564_redl_20150920_O14.fits", "/HD185859/RED_564/HD185859_w564_redl_20160813_O14.fits", "/HD148937/RED_564/HD148937_w564_redl_20170424_O14.fits", "/HD148937/RED_564/HD148937_w564_redl_20150817_O14.fits", "/HD152424/RED_564/HD152424_w564_redl_20160411_O14.fits", "/HD149404/RED_564/HD149404_w564_redl_20180630_O14.fits", "/HD149404/RED_564/HD149404_w564_redl_20180623_O14.fits", "/HD149404/RED_564/HD149404_w564_redl_20170418_O14.fits", "/HD149404/RED_564/HD149404_w564_redl_20170507_O14.fits", "/HD154043/RED_564/HD154043_w564_redl_20170502_O14.fits", "/HD168076/RED_564/HD168076_w564_redl_20180912_O14.fits", "/HD168076/RED_564/HD168076_w564_redl_20180911_O14.fits", "/HD156201/RED_564/HD156201_w564_redl_20180830_O14.fits", "/HD75860/RED_564/HD75860_w564_redl_20170429_O14.fits", "/HD75860/RED_564/HD75860_w564_redl_20170430_O14.fits", "/HD186841/RED_564/HD186841_w564_redl_20160910_O14.fits", "/HD186841/RED_564/HD186841_w564_redl_20160909_O14.fits", "/HD112272/RED_564/HD112272_w564_redl_20170430_O14.fits", "/HD112272/RED_564/HD112272_w564_redl_20170501_O14.fits", "/HD170938/RED_564/HD170938_w564_redl_20180911_O14.fits", "/HD167971/RED_564/HD167971_w564_redl_20140921_O14.fits", "/HD169454/RED_564/HD169454_w564_redl_20160714_O14.fits", "/HD169454/RED_564/HD169454_w564_redl_20160724_O14.fits", "/HD147889/RED_564/HD147889_w564_redl_20140928_O14.fits"]

SelectedListBands910 = [x for x in List910_O14x if "101065" not in x and "147084" not in x and "50820" not in x and "157978" not in x and "133518" not in x]

High_Extinction_SelectedListBands910 = [x for x in SelectedListBands910 if "183143" in x or "169454" in x or "167971" in x or "170938" in x or "186841" in x or "164740" in x or "186745" in x or "75860" in x or "156201" in x or "112272" in x or "168076" in x or "154043" in x or "73882" in x or "149404" in x or "152424" in x or "148937" in x or "161056" in x or "80558" in x or "185859" in x or "167838" in x or "153919" in x or "43384" in x or "37061" in x]
Low_Extinction_SelectedListBands910 = [x for x in SelectedListBands910 if "23016" in x or "157246" in x or "38771" in x or "148605" in x or "93030" in x or "66811" in x or "37128" in x or "158926" in x or "81188" in x]
Low_Extinction_SelectedListBands910_33AddedTargets = [x for x in SelectedListBands910 if "23016" in x or "157246" in x or "38771" in x or "148605" in x or "93030" in x or "66811" in x or "37128" in x or "158926" in x or "81188" in x or "40111" in x or "180554" in x or "57061" in x or "55879" in x or "36861" in x or "133518" in x or "34748" in x or "36822" in x or "75309" in x or "143275" in x or "49787" in x or "53975" in x or "94493" in x or "113904" in x or "135591" in x or "144470" in x or "37041" in x or "22951" in x or "104705" in x or "184915" in x or "167264" in x or "103779" in x or "91824" in x or "109399" in x or "93843" in x or "171957" in x or "145502" in x or "145502" in x or "54439" in x or "155806" in x]
"""

#C6H Bands15,16 (O13 range):
"""
plotrange = [5080, 5095]
pythia = EdiblesOracle()
List1516 = pythia.GetObsListByWavelength(5060)
List1516_O13 = [x for x in List1516 if not "O14" in x]
List1516_O13x = [x for x in List1516_O13 if not "O12" in x]
print(List1516_O13x)

ControlledTargetListNoStellarBands1516 = ["/HD37061/RED_564/HD37061_w564_redl_20190101_O13.fits", "/HD37061/RED_564/HD37061_w564_redl_20160912_O13.fits", "/HD153919/RED_564/HD153919_w564_redl_20160707_O13.fits", "/HD185859/RED_564/HD185859_w564_redl_20160814_O13.fits", "/HD185859/RED_564/HD185859_w564_redl_20150920_O13.fits", "/HD185859/RED_564/HD185859_w564_redl_20160813_O13.fits", "/HD148937/RED_564/HD148937_w564_redl_20170424_O13.fits", "/HD148937/RED_564/HD148937_w564_redl_20150817_O13.fits", "/HD152424/RED_564/HD152424_w564_redl_20160411_O13.fits", "/HD149404/RED_564/HD149404_w564_redl_20180630_O13.fits", "/HD149404/RED_564/HD149404_w564_redl_20180623_O13.fits", "/HD149404/RED_564/HD149404_w564_redl_20170418_O13.fits", "/HD149404/RED_564/HD149404_w564_redl_20170507_O13.fits", "/HD154043/RED_564/HD154043_w564_redl_20170502_O13.fits", "/HD168076/RED_564/HD168076_w564_redl_20180912_O13.fits", "/HD168076/RED_564/HD168076_w564_redl_20180911_O13.fits", "/HD156201/RED_564/HD156201_w564_redl_20180830_O13.fits", "/HD75860/RED_564/HD75860_w564_redl_20170429_O13.fits", "/HD75860/RED_564/HD75860_w564_redl_20170430_O13.fits", "/HD186841/RED_564/HD186841_w564_redl_20160910_O13.fits", "/HD186841/RED_564/HD186841_w564_redl_20160909_O13.fits", "/HD112272/RED_564/HD112272_w564_redl_20170430_O13.fits", "/HD112272/RED_564/HD112272_w564_redl_20170501_O13.fits", "/HD170938/RED_564/HD170938_w564_redl_20180911_O13.fits", "/HD167971/RED_564/HD167971_w564_redl_20140921_O13.fits", "/HD169454/RED_564/HD169454_w564_redl_20160714_O13.fits", "/HD169454/RED_564/HD169454_w564_redl_20160724_O13.fits", "/HD147889/RED_564/HD147889_w564_redl_20140928_O13.fits"]

SelectedListBands1516 = [x for x in List1516_O13x if "101065" not in x and "147084" not in x and "50820" not in x and "157978" not in x and "133518" not in x]

High_Extinction_SelectedListBands1516 = [x for x in SelectedListBands1516 if "183143" in x or "169454" in x or "167971" in x or "170938" in x or "186841" in x or "164740" in x or "186745" in x or "75860" in x or "156201" in x or "112272" in x or "168076" in x or "154043" in x or "73882" in x or "149404" in x or "152424" in x or "148937" in x or "161056" in x or "80558" in x or "185859" in x or "167838" in x or "153919" in x or "43384" in x or "37061" in x]
Low_Extinction_SelectedListBands1516 = [x for x in SelectedListBands1516 if "23016" in x or "157246" in x or "38771" in x or "148605" in x or "93030" in x or "66811" in x or "37128" in x or "158926" in x or "81188" in x]
Low_Extinction_SelectedListBands1516_33AddedTargets = [x for x in SelectedListBands1516 if "23016" in x or "157246" in x or "38771" in x or "148605" in x or "93030" in x or "66811" in x or "37128" in x or "158926" in x or "81188" in x or "40111" in x or "180554" in x or "57061" in x or "55879" in x or "36861" in x or "133518" in x or "34748" in x or "36822" in x or "75309" in x or "143275" in x or "49787" in x or "53975" in x or "94493" in x or "113904" in x or "135591" in x or "144470" in x or "37041" in x or "22951" in x or "104705" in x or "184915" in x or "167264" in x or "103779" in x or "91824" in x or "109399" in x or "93843" in x or "171957" in x or "145502" in x or "145502" in x or "54439" in x or "155806" in x]
"""


#C6H Band17 Strong HeI Line (for stellar frame), O12 Range:
"""
plotrange = [5005, 5030]
#plotrange = [5013, 5019]
#Stellar Line Check:
#SELECTED, EXTINCTION ORDERED (Increasing E(B-V)) TARGETS WITHOUT STELLAR FEATURES
ControlledTargetListNoStellarBands17 = ["/HD37061/RED_564/HD37061_w564_redl_20190101_O12.fits", "/HD37061/RED_564/HD37061_w564_redl_20160912_O12.fits", "/HD153919/RED_564/HD153919_w564_redl_20160707_O12.fits", "/HD185859/RED_564/HD185859_w564_redl_20160814_O12.fits", "/HD185859/RED_564/HD185859_w564_redl_20150920_O12.fits", "/HD185859/RED_564/HD185859_w564_redl_20160813_O12.fits", "/HD148937/RED_564/HD148937_w564_redl_20170424_O12.fits", "/HD148937/RED_564/HD148937_w564_redl_20150817_O12.fits", "/HD152424/RED_564/HD152424_w564_redl_20160411_O12.fits", "/HD149404/RED_564/HD149404_w564_redl_20180630_O12.fits", "/HD149404/RED_564/HD149404_w564_redl_20180623_O12.fits", "/HD149404/RED_564/HD149404_w564_redl_20170418_O12.fits", "/HD149404/RED_564/HD149404_w564_redl_20170507_O12.fits", "/HD154043/RED_564/HD154043_w564_redl_20170502_O12.fits", "/HD168076/RED_564/HD168076_w564_redl_20180912_O12.fits", "/HD168076/RED_564/HD168076_w564_redl_20180911_O12.fits", "/HD156201/RED_564/HD156201_w564_redl_20180830_O12.fits", "/HD75860/RED_564/HD75860_w564_redl_20170429_O12.fits", "/HD75860/RED_564/HD75860_w564_redl_20170430_O12.fits", "/HD186841/RED_564/HD186841_w564_redl_20160910_O12.fits", "/HD186841/RED_564/HD186841_w564_redl_20160909_O12.fits", "/HD112272/RED_564/HD112272_w564_redl_20170430_O12.fits", "/HD112272/RED_564/HD112272_w564_redl_20170501_O12.fits", "/HD170938/RED_564/HD170938_w564_redl_20180911_O12.fits", "/HD167971/RED_564/HD167971_w564_redl_20140921_O12.fits", "/HD169454/RED_564/HD169454_w564_redl_20160714_O12.fits", "/HD169454/RED_564/HD169454_w564_redl_20160724_O12.fits", "/HD147889/RED_564/HD147889_w564_redl_20140928_O12.fits"]
print(len(ControlledTargetListNoStellarBands17))
print(ControlledTargetListNoStellarBands17)
"""

#C6H BAND 18 LISTS:
"""
plotrange = [4820,4860]
#plotrange = [4840,4850]
pythia = EdiblesOracle()
ListBand18 = pythia.GetObsListByWavelength(4845)
ListBand18_O7 = [x for x in ListBand18 if not "O8" in x]
ListBand18_O7x = [x for x in ListBand18_O7 if not "O28" in x]
ListBand18_O7xx = [x for x in ListBand18_O7x if not "O29" in x]
print(ListBand18_O7xx)
print(len(ListBand18_O7xx))
#DE MEGA FLUCTUATING TARGETS ERUIT HALEN (HD101065, HD147084, HD50820, HD157978, 133518)!!!
SelectedListBand18 = [x for x in ListBand18_O7xx if "101065" not in x and "147084" not in x and "50820" not in x and "157978" not in x and "133518" not in x]

High_Extinction_SelectedListBand18 = [x for x in SelectedListBand18 if "183143" in x or "169454" in x or "186841" in x or "164740" in x or "186745" in x or "75860" in x or "156201" in x or "112272" in x or "168076" in x or "154043" in x or "149404" in x or "152424" in x or "148937" in x or "161056" in x or "80558" in x or "185859" in x or "153919" in x or "43384" in x or "37061" in x]
Low_Extinction_SelectedListBand18 = [x for x in SelectedListBand18 if "23016" in x or "157246" in x or "38771" in x or "148605" in x or "93030" in x or "66811" in x or "37128" in x or "158926" in x or "81188" in x]
Low_Extinction_SelectedListBand18_33AddedTargets = [x for x in SelectedListBand18 if "23016" in x or "157246" in x or "38771" in x or "148605" in x or "93030" in x or "66811" in x or "37128" in x or "158926" in x or "81188" in x or "40111" in x or "180554" in x or "57061" in x or "55879" in x or "36861" in x or "133518" in x or "34748" in x or "36822" in x or "75309" in x or "143275" in x or "49787" in x or "53975" in x or "94493" in x or "113904" in x or "135591" in x or "144470" in x or "37041" in x or "22951" in x or "104705" in x or "184915" in x or "167264" in x or "103779" in x or "91824" in x or "109399" in x or "93843" in x or "171957" in x or "145502" in x or "145502" in x or "54439" in x or "155806" in x]
Low_Extinction_SelectedListBand18_33AddedTargets_FilteredAgain = [x for x in SelectedListBand18 if "157246" in x or "37128" in x or "158926" in x or "81188" in x or "180554" in x or "133518" in x or "34748" in x or "143275" in x or "49787" in x or "135591" in x or "144470" in x or "22951" in x or "109399" in x or "93843" in x or "171957" in x or "54439" in x or "155806" in x]

#MIMICKED LIST TO BANDS1,2, SELECTED, EXTINCTION ORDERED (Increasing E(B-V)) TARGETS WITHOUT STELLAR FEATURES
ControlledTargetListNoStellarBand18 = ["/HD37061/RED_564/HD37061_w564_redl_20190101_O7.fits", "/HD37061/RED_564/HD37061_w564_redl_20160912_O7.fits", "/HD153919/RED_564/HD153919_w564_redl_20160707_O7.fits", "/HD185859/RED_564/HD185859_w564_redl_20160814_O7.fits", "/HD185859/RED_564/HD185859_w564_redl_20150920_O7.fits", "/HD185859/RED_564/HD185859_w564_redl_20160813_O7.fits", "/HD148937/RED_564/HD148937_w564_redl_20170424_O7.fits", "/HD148937/RED_564/HD148937_w564_redl_20150817_O7.fits", "/HD152424/RED_564/HD152424_w564_redl_20160411_O7.fits", "/HD149404/RED_564/HD149404_w564_redl_20180630_O7.fits", "/HD149404/RED_564/HD149404_w564_redl_20180623_O7.fits", "/HD149404/RED_564/HD149404_w564_redl_20170418_O7.fits", "/HD149404/RED_564/HD149404_w564_redl_20170507_O7.fits", "/HD154043/RED_564/HD154043_w564_redl_20170502_O7.fits", "/HD168076/RED_564/HD168076_w564_redl_20180912_O7.fits", "/HD168076/RED_564/HD168076_w564_redl_20180911_O7.fits", "/HD156201/RED_564/HD156201_w564_redl_20180830_O7.fits", "/HD75860/RED_564/HD75860_w564_redl_20170429_O7.fits", "/HD75860/RED_564/HD75860_w564_redl_20170430_O7.fits", "/HD186841/RED_564/HD186841_w564_redl_20160910_O7.fits", "/HD186841/RED_564/HD186841_w564_redl_20160909_O7.fits", "/HD112272/RED_564/HD112272_w564_redl_20170430_O7.fits", "/HD112272/RED_564/HD112272_w564_redl_20170501_O7.fits", "/HD170938/RED_564/HD170938_w564_redl_20180911_O7.fits", "/HD167971/RED_564/HD167971_w564_redl_20140921_O7.fits", "/HD169454/RED_564/HD169454_w564_redl_20160714_O7.fits", "/HD169454/RED_564/HD169454_w564_redl_20160724_O7.fits", "/HD147889/RED_564/HD147889_w564_redl_20140928_O7.fits"]
print(len(ControlledTargetListNoStellarBand18))
print(ControlledTargetListNoStellarBand18)
"""

#C6H BAND 19 LISTS:
"""
plotrange = [4730, 4760]
plotrange = [4740, 4755]
#plotrange = [4740, 4750]
pythia = EdiblesOracle()
ListBand19 = pythia.GetObsListByWavelength(4745)
ListBand19_O4 = [x for x in ListBand19 if not "O5" in x]
ListBand19_O4x = [x for x in ListBand19_O4 if not "O26" in x]
ListBand19_O4xx = [x for x in ListBand19_O4x if not "O27" in x]

#DE MEGA FLUCTUATING TARGETS ERUIT HALEN (HD101065, HD147084, HD50820, HD157978, 133518)!!!
SelectedListBand19 = [x for x in ListBand19_O4xx if "101065" not in x and "147084" not in x and "50820" not in x and "157978" not in x and "133518" not in x]
print(SelectedListBand19)
print(len(ListBand19_O4xx))
print(len(SelectedListBand19))

High_Extinction_SelectedListBand19 = [x for x in SelectedListBand19 if "183143" in x or "169454" in x or "186841" in x or "164740" in x or "186745" in x or "75860" in x or "156201" in x or "112272" in x or "168076" in x or "154043" in x or "149404" in x or "152424" in x or "148937" in x or "161056" in x or "80558" in x or "185859" in x or "153919" in x or "43384" in x or "37061" in x]
Low_Extinction_SelectedListBand19 = [x for x in SelectedListBand19 if "23016" in x or "157246" in x or "38771" in x or "148605" in x or "93030" in x or "66811" in x or "37128" in x or "158926" in x or "81188" in x]
Low_Extinction_SelectedListBand19_33AddedTargets = [x for x in SelectedListBand19 if "23016" in x or "157246" in x or "38771" in x or "148605" in x or "93030" in x or "66811" in x or "37128" in x or "158926" in x or "81188" in x or "40111" in x or "180554" in x or "57061" in x or "55879" in x or "36861" in x or "133518" in x or "34748" in x or "36822" in x or "75309" in x or "143275" in x or "49787" in x or "53975" in x or "94493" in x or "113904" in x or "135591" in x or "144470" in x or "37041" in x or "22951" in x or "104705" in x or "184915" in x or "167264" in x or "103779" in x or "91824" in x or "109399" in x or "93843" in x or "171957" in x or "145502" in x or "145502" in x or "54439" in x or "155806" in x]
Low_Extinction_SelectedListBand19_33AddedTargets_FilteredAgain = [x for x in SelectedListBand19 if "157246" in x or "37128" in x or "158926" in x or "81188" in x or "180554" in x or "133518" in x or "34748" in x or "143275" in x or "49787" in x or "135591" in x or "144470" in x or "22951" in x or "109399" in x or "93843" in x or "171957" in x or "54439" in x or "155806" in x]


#MIMICKED LIST TO BANDS1,2, SELECTED, EXTINCTION ORDERED (Increasing E(B-V)) TARGETS WITHOUT STELLAR FEATURES
ControlledTargetListNoStellarBand19 = ["/HD37061/RED_564/HD37061_w564_redl_20190101_O4.fits", "/HD37061/RED_564/HD37061_w564_redl_20160912_O4.fits", "/HD153919/RED_564/HD153919_w564_redl_20160707_O4.fits", "/HD185859/RED_564/HD185859_w564_redl_20160814_O4.fits", "/HD185859/RED_564/HD185859_w564_redl_20150920_O4.fits", "/HD185859/RED_564/HD185859_w564_redl_20160813_O4.fits", "/HD148937/RED_564/HD148937_w564_redl_20170424_O4.fits", "/HD148937/RED_564/HD148937_w564_redl_20150817_O4.fits", "/HD152424/RED_564/HD152424_w564_redl_20160411_O4.fits", "/HD149404/RED_564/HD149404_w564_redl_20180630_O4.fits", "/HD149404/RED_564/HD149404_w564_redl_20180623_O4.fits", "/HD149404/RED_564/HD149404_w564_redl_20170418_O4.fits", "/HD149404/RED_564/HD149404_w564_redl_20170507_O4.fits", "/HD154043/RED_564/HD154043_w564_redl_20170502_O4.fits", "/HD168076/RED_564/HD168076_w564_redl_20180912_O4.fits", "/HD168076/RED_564/HD168076_w564_redl_20180911_O4.fits", "/HD156201/RED_564/HD156201_w564_redl_20180830_O4.fits", "/HD75860/RED_564/HD75860_w564_redl_20170429_O4.fits", "/HD75860/RED_564/HD75860_w564_redl_20170430_O4.fits", "/HD186841/RED_564/HD186841_w564_redl_20160910_O4.fits", "/HD186841/RED_564/HD186841_w564_redl_20160909_O4.fits", "/HD112272/RED_564/HD112272_w564_redl_20170430_O4.fits", "/HD112272/RED_564/HD112272_w564_redl_20170501_O4.fits", "/HD170938/RED_564/HD170938_w564_redl_20180911_O4.fits", "/HD167971/RED_564/HD167971_w564_redl_20140921_O4.fits", "/HD169454/RED_564/HD169454_w564_redl_20160714_O4.fits", "/HD169454/RED_564/HD169454_w564_redl_20160724_O4.fits", "/HD147889/RED_564/HD147889_w564_redl_20140928_O4.fits"]
ControlledTargetListNoStellarBand19_O5 = ["/HD37061/RED_564/HD37061_w564_redl_20190101_O5.fits", "/HD37061/RED_564/HD37061_w564_redl_20160912_O5.fits", "/HD153919/RED_564/HD153919_w564_redl_20160707_O5.fits", "/HD185859/RED_564/HD185859_w564_redl_20160814_O5.fits", "/HD185859/RED_564/HD185859_w564_redl_20150920_O5.fits", "/HD185859/RED_564/HD185859_w564_redl_20160813_O5.fits", "/HD148937/RED_564/HD148937_w564_redl_20170424_O5.fits", "/HD148937/RED_564/HD148937_w564_redl_20150817_O5.fits", "/HD152424/RED_564/HD152424_w564_redl_20160411_O5.fits", "/HD149404/RED_564/HD149404_w564_redl_20180630_O5.fits", "/HD149404/RED_564/HD149404_w564_redl_20180623_O5.fits", "/HD149404/RED_564/HD149404_w564_redl_20170418_O5.fits", "/HD149404/RED_564/HD149404_w564_redl_20170507_O5.fits", "/HD154043/RED_564/HD154043_w564_redl_20170502_O5.fits", "/HD168076/RED_564/HD168076_w564_redl_20180912_O5.fits", "/HD168076/RED_564/HD168076_w564_redl_20180911_O5.fits", "/HD156201/RED_564/HD156201_w564_redl_20180830_O5.fits", "/HD75860/RED_564/HD75860_w564_redl_20170429_O5.fits", "/HD75860/RED_564/HD75860_w564_redl_20170430_O5.fits", "/HD186841/RED_564/HD186841_w564_redl_20160910_O5.fits", "/HD186841/RED_564/HD186841_w564_redl_20160909_O5.fits", "/HD112272/RED_564/HD112272_w564_redl_20170430_O5.fits", "/HD112272/RED_564/HD112272_w564_redl_20170501_O5.fits", "/HD170938/RED_564/HD170938_w564_redl_20180911_O5.fits", "/HD167971/RED_564/HD167971_w564_redl_20140921_O5.fits", "/HD169454/RED_564/HD169454_w564_redl_20160714_O5.fits", "/HD169454/RED_564/HD169454_w564_redl_20160724_O5.fits", "/HD147889/RED_564/HD147889_w564_redl_20140928_O5.fits"]
print(len(ControlledTargetListNoStellarBand19))
print(ControlledTargetListNoStellarBand19)
ControlledTargetListNoStellarBand19ExtraFilter = ["/HD37061/RED_564/HD37061_w564_redl_20190101_O4.fits", "/HD37061/RED_564/HD37061_w564_redl_20160912_O4.fits", "/HD153919/RED_564/HD153919_w564_redl_20160707_O4.fits", "/HD185859/RED_564/HD185859_w564_redl_20160814_O4.fits", "/HD185859/RED_564/HD185859_w564_redl_20150920_O4.fits", "/HD185859/RED_564/HD185859_w564_redl_20160813_O4.fits", "/HD148937/RED_564/HD148937_w564_redl_20170424_O4.fits", "/HD148937/RED_564/HD148937_w564_redl_20150817_O4.fits", "/HD152424/RED_564/HD152424_w564_redl_20160411_O4.fits", "/HD149404/RED_564/HD149404_w564_redl_20180630_O4.fits", "/HD149404/RED_564/HD149404_w564_redl_20180623_O4.fits", "/HD149404/RED_564/HD149404_w564_redl_20170418_O4.fits", "/HD149404/RED_564/HD149404_w564_redl_20170507_O4.fits", "/HD154043/RED_564/HD154043_w564_redl_20170502_O4.fits", "/HD168076/RED_564/HD168076_w564_redl_20180912_O4.fits", "/HD168076/RED_564/HD168076_w564_redl_20180911_O4.fits", "/HD156201/RED_564/HD156201_w564_redl_20180830_O4.fits", "/HD75860/RED_564/HD75860_w564_redl_20170429_O4.fits", "/HD75860/RED_564/HD75860_w564_redl_20170430_O4.fits", "/HD186841/RED_564/HD186841_w564_redl_20160910_O4.fits", "/HD186841/RED_564/HD186841_w564_redl_20160909_O4.fits", "/HD112272/RED_564/HD112272_w564_redl_20170430_O4.fits", "/HD112272/RED_564/HD112272_w564_redl_20170501_O4.fits", "/HD170938/RED_564/HD170938_w564_redl_20180911_O4.fits", "/HD169454/RED_564/HD169454_w564_redl_20160714_O4.fits", "/HD169454/RED_564/HD169454_w564_redl_20160724_O4.fits"]
print(len(ControlledTargetListNoStellarBand19ExtraFilter))
"""


#C6H BANDS 1,2 LISTS:
#"""
#plotrange = [5230, 5300]
#plotrange = [5230, 5275]
#plotrange = [5250,5280]
#plotrange = [5255, 5275]
#plotrange = [5230,5240]
#plotrange = [5240,5250]
#plotrange = [5250,5260]
plotrange = [5260,5270]
#plotrange = [5270,5280]
#plotrange = [5280,5290]
#plotrange = [5290,5300]
#plotrange = [5262, 5267]
#plotrange = [5265, 5267]
pythia = EdiblesOracle()
List = pythia.GetObsListByWavelength(5265.5)
ReverseList = List[::-1]
ReverseListO17 = [x for x in ReverseList if not "O18" in x]

#DE MEGA FLUCTUATING TARGETS ERUIT HALEN (HD101065, HD147084, HD50820, HD157978, 133518)!!!
SelectedListO17 = [x for x in ReverseListO17 if "101065" not in x and "147084" not in x and "50820" not in x and "157978" not in x and "133518" not in x and "152408" not in x]
print(len(ReverseListO17))
print(len(SelectedListO17))

#Oscillating ones filtered out In General, ook HD164073
SelectedListO17x = [x for x in SelectedListO17 if "/HD73882/RED_564/HD73882_w564_redl_20170425_O17.fits" not in x and "/HD149404/RED_564/HD149404_w564_redl_20170417_O17.fits" not in x and "164073" not in x]
print(len(SelectedListO17x))


#HIGH EXTINCTION ONLY TARGET LIST:
High_Extinction_SelectedListO17 = [x for x in SelectedListO17 if "183143" in x or "147889" in x or "169454" in x or "167971" in x or "170938" in x or "186841" in x or "164740" in x or "186745" in x or "75860" in x or "156201" in x or "112272" in x or "61827" in x or "165319" in x or "168076" in x or "154043" in x or "73882" in x or "149404" in x or "152424" in x or "148937" in x or "161056" in x or "80558" in x or "185859" in x or "167838" in x or "153919" in x or "43384" in x or "37061" in x]
#SelectedListO17 = [x for x in ListO17 if "183143" in x or "147889" in x or "169454" in x or "167971" in x or "170938" in x or "186841" in x or "164740" in x or "186745" in x or "75860" in x or "156201" in x or "112272" in x or "61827" in x or "165319" in x or "168076" in x or "154043" in x or "73882" in x or "149404" in x or "152424" in x or "148937" in x or "157978" in x or "161056" in x or "80558" in x or "185859" in x or "167838" in x or "153919" in x or "43384" in x or "37061" in x]
#"Too Variable" Targetlist: HD101065, HD147084, HD50820, HD157978, HD133518, HD152408
print(len(High_Extinction_SelectedListO17))

#FILTERING OUT THE WEIRD OSCILATTING ONES:
ExtraFiltered_High_Extinction_SelectedListO17 = [x for x in High_Extinction_SelectedListO17 if "/HD149404/RED_564/HD149404_w564_redl_20170417_O17.fits" not in x and "73882" not in x]
print(len(ExtraFiltered_High_Extinction_SelectedListO17))

#FILTERING OUT THE BAD S/N RATIO ONES:
ExtraFiltered2_High_Extinction_SelectedListO17 = [x for x in ExtraFiltered_High_Extinction_SelectedListO17 if "/HD154043/RED_564/HD154043_w564_redl_20170718_O17.fits" not in x and "/HD154043/RED_564/HD154043_w564_redl_20170719_O17.fits" not in x and "/HD147889/RED_564/HD147889_w564_redl_20140927_O17.fits" not in x and "/HD147889/RED_564/HD147889_w564_redl_20150816_O17.fits" not in x and "/HD61827/RED_564/HD61827_w564_redl_20180325_O17.fits" not in x and "165319" not in x]
print(len(ExtraFiltered2_High_Extinction_SelectedListO17))
print(ExtraFiltered2_High_Extinction_SelectedListO17)

#SELECTED, EXTINCTION ORDERED (Increasing E(B-V)) TARGETS WITH STELLAR FEATURES
ControlledTargetListStellar = ["/HD43384/RED_564/HD43384_w564_redl_20150930_O17.fits", "/HD167838/RED_564/HD167838_w564_redl_20180902_O17.fits", "/HD167838/RED_564/HD167838_w564_redl_20160808_O17.fits", "/HD80558/RED_564/HD80558_w564_redl_20160221_O17.fits", "/HD80558/RED_564/HD80558_w564_redl_20150603_O17.fits", "/HD61827/RED_564/HD61827_w564_redl_20180324_O17.fits", "/HD186745/RED_564/HD186745_w564_redl_20150920_O17.fits", "/HD186745/RED_564/HD186745_w564_redl_20160909_O17.fits", "/HD183143/RED_564/HD183143_w564_redl_20180903_O17.fits", "/HD183143/RED_564/HD183143_w564_redl_20180907_O17.fits"]
print(len(ControlledTargetListStellar))

#SELECTED, EXTINCTION NON-ORDERED (Increasing E(B-V)) TARGETS WITHOUT STELLAR FEATURES
ExtraFiltered3_High_Extinction_SelectedListO17 = [x for x in ExtraFiltered2_High_Extinction_SelectedListO17 if "/HD43384/RED_564/HD43384_w564_redl_20150930_O17.fits" not in x and "/HD167838/RED_564/HD167838_w564_redl_20180902_O17.fits" not in x and "/HD167838/RED_564/HD167838_w564_redl_20160808_O17.fits" not in x and "/HD80558/RED_564/HD80558_w564_redl_20150603_O17.fits" not in x and "/HD186745/RED_564/HD186745_w564_redl_20150920_O17.fits" not in x and "/HD186745/RED_564/HD186745_w564_redl_20160909_O17.fits" not in x and "/HD183143/RED_564/HD183143_w564_redl_20180903_O17.fits" not in x and "/HD183143/RED_564/HD183143_w564_redl_20180907_O17.fits" not in x and "/HD80558/RED_564/HD80558_w564_redl_20160221_O17.fits" not in x and "/HD61827/RED_564/HD61827_w564_redl_20180324_O17.fits" not in x]
print(len(ExtraFiltered3_High_Extinction_SelectedListO17))
print(ExtraFiltered3_High_Extinction_SelectedListO17)

#SELECTED, EXTINCTION ORDERED (Increasing E(B-V)) TARGETS WITHOUT STELLAR FEATURES
ControlledTargetListNoStellar = ["/HD37061/RED_564/HD37061_w564_redl_20190101_O17.fits", "/HD37061/RED_564/HD37061_w564_redl_20160912_O17.fits", "/HD153919/RED_564/HD153919_w564_redl_20160707_O17.fits", "/HD185859/RED_564/HD185859_w564_redl_20160814_O17.fits", "/HD185859/RED_564/HD185859_w564_redl_20150920_O17.fits", "/HD185859/RED_564/HD185859_w564_redl_20160813_O17.fits", "/HD148937/RED_564/HD148937_w564_redl_20170424_O17.fits", "/HD148937/RED_564/HD148937_w564_redl_20150817_O17.fits", "/HD152424/RED_564/HD152424_w564_redl_20160411_O17.fits", "/HD149404/RED_564/HD149404_w564_redl_20180630_O17.fits", "/HD149404/RED_564/HD149404_w564_redl_20180623_O17.fits", "/HD149404/RED_564/HD149404_w564_redl_20170418_O17.fits", "/HD149404/RED_564/HD149404_w564_redl_20170507_O17.fits", "/HD154043/RED_564/HD154043_w564_redl_20170502_O17.fits", "/HD168076/RED_564/HD168076_w564_redl_20180912_O17.fits", "/HD168076/RED_564/HD168076_w564_redl_20180911_O17.fits", "/HD156201/RED_564/HD156201_w564_redl_20180830_O17.fits", "/HD75860/RED_564/HD75860_w564_redl_20170429_O17.fits", "/HD75860/RED_564/HD75860_w564_redl_20170430_O17.fits", "/HD186841/RED_564/HD186841_w564_redl_20160910_O17.fits", "/HD186841/RED_564/HD186841_w564_redl_20160909_O17.fits", "/HD112272/RED_564/HD112272_w564_redl_20170430_O17.fits", "/HD112272/RED_564/HD112272_w564_redl_20170501_O17.fits", "/HD170938/RED_564/HD170938_w564_redl_20180911_O17.fits", "/HD167971/RED_564/HD167971_w564_redl_20140921_O17.fits", "/HD169454/RED_564/HD169454_w564_redl_20160714_O17.fits", "/HD169454/RED_564/HD169454_w564_redl_20160724_O17.fits", "/HD147889/RED_564/HD147889_w564_redl_20140928_O17.fits"]
#ControlledTargetListNoStellar_O18 = ["/HD37061/RED_564/HD37061_w564_redl_20190101_O18.fits", "/HD37061/RED_564/HD37061_w564_redl_20160912_O18.fits", "/HD153919/RED_564/HD153919_w564_redl_20160707_O18.fits", "/HD185859/RED_564/HD185859_w564_redl_20160814_O18.fits", "/HD185859/RED_564/HD185859_w564_redl_20150920_O18.fits", "/HD185859/RED_564/HD185859_w564_redl_20160813_O18.fits", "/HD148937/RED_564/HD148937_w564_redl_20170424_O18.fits", "/HD148937/RED_564/HD148937_w564_redl_20150817_O18.fits", "/HD152424/RED_564/HD152424_w564_redl_20160411_O18.fits", "/HD149404/RED_564/HD149404_w564_redl_20180630_O18.fits", "/HD149404/RED_564/HD149404_w564_redl_20180623_O18.fits", "/HD149404/RED_564/HD149404_w564_redl_20170418_O18.fits", "/HD149404/RED_564/HD149404_w564_redl_20170507_O18.fits", "/HD154043/RED_564/HD154043_w564_redl_20170502_O18.fits", "/HD168076/RED_564/HD168076_w564_redl_20180912_O18.fits", "/HD168076/RED_564/HD168076_w564_redl_20180911_O18.fits", "/HD156201/RED_564/HD156201_w564_redl_20180830_O18.fits", "/HD75860/RED_564/HD75860_w564_redl_20170429_O18.fits", "/HD75860/RED_564/HD75860_w564_redl_20170430_O18.fits", "/HD186841/RED_564/HD186841_w564_redl_20160910_O18.fits", "/HD186841/RED_564/HD186841_w564_redl_20160909_O18.fits", "/HD112272/RED_564/HD112272_w564_redl_20170430_O18.fits", "/HD112272/RED_564/HD112272_w564_redl_20170501_O18.fits", "/HD170938/RED_564/HD170938_w564_redl_20180911_O18.fits", "/HD167971/RED_564/HD167971_w564_redl_20140921_O18.fits", "/HD169454/RED_564/HD169454_w564_redl_20160714_O18.fits", "/HD169454/RED_564/HD169454_w564_redl_20160724_O18.fits", "/HD147889/RED_564/HD147889_w564_redl_20140928_O18.fits"]
print(len(ControlledTargetListNoStellar))
print(ControlledTargetListNoStellar)
#Just the 112272´s out:
ControlledTargetListNoStellar1ExtraFilter = ["/HD37061/RED_564/HD37061_w564_redl_20190101_O17.fits", "/HD37061/RED_564/HD37061_w564_redl_20160912_O17.fits", "/HD153919/RED_564/HD153919_w564_redl_20160707_O17.fits", "/HD185859/RED_564/HD185859_w564_redl_20160814_O17.fits", "/HD185859/RED_564/HD185859_w564_redl_20150920_O17.fits", "/HD185859/RED_564/HD185859_w564_redl_20160813_O17.fits", "/HD148937/RED_564/HD148937_w564_redl_20170424_O17.fits", "/HD148937/RED_564/HD148937_w564_redl_20150817_O17.fits", "/HD152424/RED_564/HD152424_w564_redl_20160411_O17.fits", "/HD149404/RED_564/HD149404_w564_redl_20180630_O17.fits", "/HD149404/RED_564/HD149404_w564_redl_20180623_O17.fits", "/HD149404/RED_564/HD149404_w564_redl_20170418_O17.fits", "/HD149404/RED_564/HD149404_w564_redl_20170507_O17.fits", "/HD154043/RED_564/HD154043_w564_redl_20170502_O17.fits", "/HD168076/RED_564/HD168076_w564_redl_20180912_O17.fits", "/HD168076/RED_564/HD168076_w564_redl_20180911_O17.fits", "/HD156201/RED_564/HD156201_w564_redl_20180830_O17.fits", "/HD75860/RED_564/HD75860_w564_redl_20170429_O17.fits", "/HD75860/RED_564/HD75860_w564_redl_20170430_O17.fits", "/HD186841/RED_564/HD186841_w564_redl_20160910_O17.fits", "/HD186841/RED_564/HD186841_w564_redl_20160909_O17.fits", "/HD170938/RED_564/HD170938_w564_redl_20180911_O17.fits", "/HD167971/RED_564/HD167971_w564_redl_20140921_O17.fits", "/HD169454/RED_564/HD169454_w564_redl_20160714_O17.fits", "/HD169454/RED_564/HD169454_w564_redl_20160724_O17.fits", "/HD147889/RED_564/HD147889_w564_redl_20140928_O17.fits"]
#Also 152424 out, Target List of 14:
ControlledTargetListNoStellar1ExtraFilter2 = ["/HD37061/RED_564/HD37061_w564_redl_20190101_O17.fits", "/HD37061/RED_564/HD37061_w564_redl_20160912_O17.fits", "/HD153919/RED_564/HD153919_w564_redl_20160707_O17.fits", "/HD185859/RED_564/HD185859_w564_redl_20160814_O17.fits", "/HD185859/RED_564/HD185859_w564_redl_20150920_O17.fits", "/HD185859/RED_564/HD185859_w564_redl_20160813_O17.fits", "/HD148937/RED_564/HD148937_w564_redl_20170424_O17.fits", "/HD148937/RED_564/HD148937_w564_redl_20150817_O17.fits", "/HD149404/RED_564/HD149404_w564_redl_20180630_O17.fits", "/HD149404/RED_564/HD149404_w564_redl_20180623_O17.fits", "/HD149404/RED_564/HD149404_w564_redl_20170418_O17.fits", "/HD149404/RED_564/HD149404_w564_redl_20170507_O17.fits", "/HD154043/RED_564/HD154043_w564_redl_20170502_O17.fits", "/HD168076/RED_564/HD168076_w564_redl_20180912_O17.fits", "/HD168076/RED_564/HD168076_w564_redl_20180911_O17.fits", "/HD156201/RED_564/HD156201_w564_redl_20180830_O17.fits", "/HD75860/RED_564/HD75860_w564_redl_20170429_O17.fits", "/HD75860/RED_564/HD75860_w564_redl_20170430_O17.fits", "/HD186841/RED_564/HD186841_w564_redl_20160910_O17.fits", "/HD186841/RED_564/HD186841_w564_redl_20160909_O17.fits", "/HD170938/RED_564/HD170938_w564_redl_20180911_O17.fits", "/HD167971/RED_564/HD167971_w564_redl_20140921_O17.fits", "/HD169454/RED_564/HD169454_w564_redl_20160714_O17.fits", "/HD169454/RED_564/HD169454_w564_redl_20160724_O17.fits", "/HD147889/RED_564/HD147889_w564_redl_20140928_O17.fits"]

#EXTRA EXTINCTION TARGETS ADDED 0.4-0.5E, TargetList of 13!!!
High_Extinction_SelectedListO17_Newlist = [x for x in SelectedListO17 if"169454" in x or "170938" in x or "186841" in x or "164740" in x or "75860" in x or "156201" in x or "112272" in x or "168076" in x or "148937" in x or "161056" in x or "185859" in x or "153919" in x or "37061" in x]
High_Extinction_SelectedListO17_Newlist_12AddedTargets = [x for x in SelectedListO17 if"169454" in x or "170938" in x or "186841" in x or "164740" in x or "75860" in x or "156201" in x or "112272" in x or "168076" in x or "148937" in x or "161056" in x or "185859" in x or "153919" in x or "37061" in x or "124314" in x or "170740" in x or "99953" in x or "164906" in x or "147933" in x or "45314" in x or "303308" in x or "185418" in x or "41117" in x]

#LOW EXTINCTION TARGET LIST (Starting with 9 targets of E <= 0.1):
Low_Extinction_SelectedListO17 = [x for x in SelectedListO17 if "23016" in x or "157246" in x or "38771" in x or "148605" in x or "93030" in x or "66811" in x or "37128" in x or "158926" in x or "81188" in x]
#Lambda Sco Out:
Low_Extinction_SelectedListO17x = [x for x in SelectedListO17 if "23016" in x or "157246" in x or "38771" in x or "148605" in x or "93030" in x or "66811" in x or "37128" in x or "81188" in x]
#Extra's:
Low_Extinction_SelectedListO17_6AddedTargets = [x for x in SelectedListO17 if "23016" in x or "157246" in x or "38771" in x or "148605" in x or "93030" in x or "66811" in x or "37128" in x or "158926" in x or "81188" in x or "40111" in x or "180554" in x or "57061" in x or "55879" in x or "36861" in x or "133518" in x]
Low_Extinction_SelectedListO17_13AddedTargets = [x for x in SelectedListO17 if "23016" in x or "157246" in x or "38771" in x or "148605" in x or "93030" in x or "66811" in x or "37128" in x or "158926" in x or "81188" in x or "40111" in x or "180554" in x or "57061" in x or "55879" in x or "36861" in x or "133518" in x or "34748" in x or "36822" in x or "75309" in x or "143275" in x or "49787" in x or "53975" in x or "94493" in x]
Low_Extinction_SelectedListO17_33AddedTargets = [x for x in SelectedListO17 if "23016" in x or "157246" in x or "38771" in x or "148605" in x or "93030" in x or "66811" in x or "37128" in x or "158926" in x or "81188" in x or "40111" in x or "180554" in x or "57061" in x or "55879" in x or "36861" in x or "133518" in x or "34748" in x or "36822" in x or "75309" in x or "143275" in x or "49787" in x or "53975" in x or "94493" in x or "113904" in x or "135591" in x or "144470" in x or "37041" in x or "22951" in x or "104705" in x or "184915" in x or "167264" in x or "103779" in x or "91824" in x or "109399" in x or "93843" in x or "171957" in x or "145502" in x or "145502" in x or "54439" in x or "155806" in x]
Low_Extinction_SelectedListO17_33AddedTargets_FilteredAgain = [x for x in SelectedListO17 if "157246" in x or "37128" in x or "158926" in x or "81188" in x or "180554" in x or "133518" in x or "34748" in x or "143275" in x or "49787" in x or "135591" in x or "144470" in x or "22951" in x or "109399" in x or "93843" in x or "171957" in x or "54439" in x or "155806" in x]

#STELLAR FRAME 8 TARGETS WITH STELLAR LINES STACKED (Also the ones without just to check):
ControlledTargetListNoStellar_StellarRestFrameShifted_8Targets = ["/HD185859/RED_564/HD185859_w564_redl_20160814_O17.fits", "/HD185859/RED_564/HD185859_w564_redl_20150920_O17.fits", "/HD185859/RED_564/HD185859_w564_redl_20160813_O17.fits", "/HD154043/RED_564/HD154043_w564_redl_20170502_O17.fits", "/HD156201/RED_564/HD156201_w564_redl_20180830_O17.fits", "/HD75860/RED_564/HD75860_w564_redl_20170429_O17.fits", "/HD75860/RED_564/HD75860_w564_redl_20170430_O17.fits", "/HD186841/RED_564/HD186841_w564_redl_20160910_O17.fits", "/HD186841/RED_564/HD186841_w564_redl_20160909_O17.fits", "/HD112272/RED_564/HD112272_w564_redl_20170430_O17.fits", "/HD112272/RED_564/HD112272_w564_redl_20170501_O17.fits", "/HD170938/RED_564/HD170938_w564_redl_20180911_O17.fits", "/HD169454/RED_564/HD169454_w564_redl_20160714_O17.fits", "/HD169454/RED_564/HD169454_w564_redl_20160724_O17.fits"]
ControlledTargetListNoStellar_TheRest_Also8Targets = ["/HD37061/RED_564/HD37061_w564_redl_20190101_O17.fits", "/HD37061/RED_564/HD37061_w564_redl_20160912_O17.fits", "/HD153919/RED_564/HD153919_w564_redl_20160707_O17.fits", "/HD148937/RED_564/HD148937_w564_redl_20170424_O17.fits", "/HD148937/RED_564/HD148937_w564_redl_20150817_O17.fits", "/HD152424/RED_564/HD152424_w564_redl_20160411_O17.fits", "/HD149404/RED_564/HD149404_w564_redl_20180630_O17.fits", "/HD149404/RED_564/HD149404_w564_redl_20180623_O17.fits", "/HD149404/RED_564/HD149404_w564_redl_20170418_O17.fits", "/HD149404/RED_564/HD149404_w564_redl_20170507_O17.fits", "/HD168076/RED_564/HD168076_w564_redl_20180912_O17.fits", "/HD168076/RED_564/HD168076_w564_redl_20180911_O17.fits", "/HD167971/RED_564/HD167971_w564_redl_20140921_O17.fits", "/HD147889/RED_564/HD147889_w564_redl_20140928_O17.fits"]

#NEW LOW EXTINCTION FOLLOW UP LISTS, NOW STARTING WITH A LIST WITH ~20TARGETS ASWELL:
#E<0.2 targets.
#HD148605... en HD34748.. eruit, slechtste S/N
###15 TARGETS:
LowExList1 = [x for x in SelectedListO17 if "23016" in x or "157246" in x or "38771" in x or "93030" in x or "66811" in x or "37128" in x or "158926" in x or "81188" in x or "40111" in x or "180554" in x or "57061" in x or "55879" in x or "36861" in x or "133518" in x or "36822" in x or "75309" in x or "143275" in x]
#Nu alles E<0.3:
#Filtered Out: HD116852, bad S/N, HD166937, HD164073, weird oscillating
#Verder ook nog HD23180 eruit ivm 1 slechte S/N van de 5 en verder ook nog HD147683 eruit
###39 TARGETS: 
LowExList2 = [x for x in SelectedListO17 if "23016" in x or "157246" in x or "38771" in x or "93030" in x or "66811" in x or "37128" in x or "158926" in x or "81188" in x or "40111" in x or "180554" in x or "57061" in x or "55879" in x or "36861" in x or "133518" in x or "36822" in x or "75309" in x or "143275" in x or "49787" in x or "53975" in x or "94493" in x or "113904" in x or "135591" in x or "144470" in x or "37041" in x or "22951" in x or "104705" in x or "184915" in x or "167264" in x or "103779" in x or "91824" in x or "109399" in x or "93843" in x or "171957" in x or "145502" in x or "145502" in x or "54439" in x or "155806" in x or "114886" in x or "38087" in x or "79186" in x or "37020" in x or "24398" in x]
#Nu alles E<0.4:
#Filtered Out: HD36982, HD39680, HD151804, HD149038, bad S/N
###48 TARGETS: 
LowExList3 = [x for x in SelectedListO17 if "23016" in x or "157246" in x or "38771" in x or "93030" in x or "66811" in x or "37128" in x or "158926" in x or "81188" in x or "40111" in x or "180554" in x or "57061" in x or "55879" in x or "36861" in x or "133518" in x or "36822" in x or "75309" in x or "143275" in x or "49787" in x or "53975" in x or "94493" in x or "113904" in x or "135591" in x or "144470" in x or "37041" in x or "22951" in x or "104705" in x or "184915" in x or "167264" in x or "103779" in x or "91824" in x or "109399" in x or "93843" in x or "171957" in x or "145502" in x or "145502" in x or "54439" in x or "155806" in x or "114886" in x or "38087" in x or "79186" in x or "37020" in x or "24398" in x or "37021" in x or "149757" in x or "54662" in x or "111934" in x or "122879" in x  or "37903" in x or "93222" in x or "27778" in x or "37023" in x or "37367" in x or "147165" in x]
#Nu alles E<0.6:
#Filtered Out: HD147888, HD150136, HD80558, HD43384, HD167838, HD170740
###59 TARGETS:
LowExList4 = [x for x in SelectedListO17 if "23016" in x or "157246" in x or "38771" in x or "93030" in x or "66811" in x or "37128" in x or "158926" in x or "81188" in x or "40111" in x or "180554" in x or "57061" in x or "55879" in x or "36861" in x or "133518" in x or "36822" in x or "75309" in x or "143275" in x or "49787" in x or "53975" in x or "94493" in x or "113904" in x or "135591" in x or "144470" in x or "37041" in x or "22951" in x or "104705" in x or "184915" in x or "167264" in x or "103779" in x or "91824" in x or "109399" in x or "93843" in x or "171957" in x or "145502" in x or "145502" in x or "54439" in x or "155806" in x or "114886" in x or "38087" in x or "79186" in x or "37020" in x or "24398" in x or "37021" in x or "149757" in x or "54662" in x or "111934" in x or "122879" in x  or "37903" in x or "93222" in x or "27778" in x or "37023" in x or "37367" in x or "147165" in x or "41117" in x or "185418" in x or "303308" in x or "45314" in x or "147933" in x or "164906" in x or "152408" in x or "99953" in x or "124314" in x or "37061" in x or "153919" in x or "185859" in x]
#Nu alles E<1.1, alles dus. Co-add je main High Extinction Filtered List hier gwn:
#Nog extra filtered toch: HD154043, HD147889
###70 TARGETS:
LowExList5 = [x for x in SelectedListO17 if "23016" in x or "157246" in x or "38771" in x or "93030" in x or "66811" in x or "37128" in x or "158926" in x or "81188" in x or "40111" in x or "180554" in x or "57061" in x or "55879" in x or "36861" in x or "133518" in x or "36822" in x or "75309" in x or "143275" in x or "49787" in x or "53975" in x or "94493" in x or "113904" in x or "135591" in x or "144470" in x or "37041" in x or "22951" in x or "104705" in x or "184915" in x or "167264" in x or "103779" in x or "91824" in x or "109399" in x or "93843" in x or "171957" in x or "145502" in x or "145502" in x or "54439" in x or "155806" in x or "114886" in x or "38087" in x or "79186" in x or "37020" in x or "24398" in x or "37021" in x or "149757" in x or "54662" in x or "111934" in x or "122879" in x  or "37903" in x or "93222" in x or "27778" in x or "37023" in x or "37367" in x or "147165" in x or "41117" in x or "185418" in x or "303308" in x or "45314" in x or "147933" in x or "164906" in x or "152408" in x or "99953" in x or "124314" in x or "37061" in x or "153919" in x or "185859" in x or "148937" in x or "152424" in x or "149404" in x or "168076" in x or "156201" in x or "75860" in x or "186841" in x or "112272" in x or "170938" in x or "167971" in x or "169454" in x]

#Telluric Lists:
#Filtered Out List:
telluric = PYTHONDIR + "/edibles/data/Telluric_Subset_Main_HDList.dat"
t = pd.read_csv(telluric, delim_whitespace=True)
print(t)
tx = t.target
print(t.target)
Telluric_TargetList = tx.values.tolist()
print(Telluric_TargetList)

Less_Telluric = [x for x in SelectedListO17x if x not in Telluric_TargetList]
print(Less_Telluric)
print(len(Less_Telluric))

#16Target High Extinction Target List, aiming to filter out the most prominent telluric emission feature ones:
ControlledTargetListNoStellarxxx = ["/HD37061/RED_564/HD37061_w564_redl_20160912_O17.fits", "/HD153919/RED_564/HD153919_w564_redl_20160707_O17.fits", "/HD148937/RED_564/HD148937_w564_redl_20170424_O17.fits", "/HD148937/RED_564/HD148937_w564_redl_20150817_O17.fits", "/HD152424/RED_564/HD152424_w564_redl_20160411_O17.fits", "/HD168076/RED_564/HD168076_w564_redl_20180912_O17.fits", "/HD168076/RED_564/HD168076_w564_redl_20180911_O17.fits", "/HD156201/RED_564/HD156201_w564_redl_20180830_O17.fits", "/HD186841/RED_564/HD186841_w564_redl_20160910_O17.fits", "/HD186841/RED_564/HD186841_w564_redl_20160909_O17.fits", "/HD170938/RED_564/HD170938_w564_redl_20180911_O17.fits", "/HD167971/RED_564/HD167971_w564_redl_20140921_O17.fits", "/HD169454/RED_564/HD169454_w564_redl_20160714_O17.fits", "/HD169454/RED_564/HD169454_w564_redl_20160724_O17.fits", "/HD147889/RED_564/HD147889_w564_redl_20140928_O17.fits"]
#"""




#ISM CLOUD VELOCITIES:
cloudvelocitiesdffilename = PYTHONDIR + '/edibles/data/ISM_Cloud_Velocities_inKms_Approximate_Database_withtargetsasindex.dat'
cloudvelocitiesdf = pd.read_csv(cloudvelocitiesdffilename, delim_whitespace=True)

#Stellar Rest Frames Database:
stellarrestframedffilename = PYTHONDIR + '/edibles/data/Stellar_Rest_Frame_Velocities_OII5160_withtargetsasindex.txt'
stellarrestframedf = pd.read_csv(stellarrestframedffilename, delim_whitespace=True)
print(stellarrestframedf.stellar_rest_velocity)

#FOR-LOOP QUICK STACKING:
totalplotwave = np.linspace((plotrange[0]+0.01), (plotrange[1]-0.01), ((plotrange[1]-plotrange[0])*50))
#print(totalplotwave)
totalnormfluxlistism = []
totalnormfluxlistgeo = []
totalnormfluxlistbary = []
totalnormfluxliststellar = []
counter = 0
trial = 0.01
trial = 0
for filename in ControlledTargetListNoStellar: 
	counter = counter + 1
	print(counter)
	sp = EdiblesSpectrum(filename)
	print(sp.target)
	print(filename)
	
	if sp.target not in cloudvelocitiesdf.index:
		cloudvelocityoftarget = 0
		print("De sp.target naam komt niet overeen met je naam in de Cloud Velocity Database (bij verschillende sp.target names per 1 Sightline), hier kan je niks mee, is nu in Barycentric frame")
	else:
		cloudvelocityoftarget = cloudvelocitiesdf.cloud_velocity[sp.target]
		print("Cloud Velocity of {} (in km/s) =".format(sp.target), cloudvelocityoftarget)
	
	if sp.target not in stellarrestframedf.index:
		stellarrestvelocityoftarget = 0
		print("Deze Had geen Stellar Line, dus onverschoven")
	else:
		stellarrestvelocityoftarget = stellarrestframedf.stellar_rest_velocity[sp.target]
		print("Stellar Rest Velocity of {} (in km/s) =".format(sp.target), stellarrestvelocityoftarget)
		
	wave = sp.wave
	#print("geo_wave", wave[400])
	flux = np.clip(sp.flux, 0, None) 
	bool_keep_geo = (wave > plotrange[0]) & (wave < plotrange[1])
	plotgeo_wave = wave[bool_keep_geo]
	plotgeoflux = flux[bool_keep_geo]
	#print("Comparison Geo Flux Values", plotgeoflux[:1])
	normgeoflux = plotgeoflux / np.median(plotgeoflux)
	trialplotnormgeoflux = normgeoflux +trial
	
	bary_wave = sp.bary_wave
	#print("bary_wave", bary_wave[400])
	bool_keep_bary = (bary_wave > plotrange[0]) & (bary_wave < plotrange[1])
	plotbary_wave = bary_wave[bool_keep_bary]
	plotbaryflux = flux[bool_keep_bary]
	#print ("Comparison Bary Flux Values", plotbaryflux[:1])
	normbaryflux = plotbaryflux / np.median(plotbaryflux)
	trialplotnormbaryflux = normbaryflux + trial
	
	ism_wave = bary_wave + (bary_wave * (cloudvelocityoftarget/(3*10**5)))
	#print("ism_wave", ism_wave[400])
	bool_keep_ism = (ism_wave > plotrange[0]) & (ism_wave < plotrange[1])
	plotism_wave = ism_wave[bool_keep_ism]
	plotismflux = flux[bool_keep_ism]
	#print("ISM Flux Values", plotismflux[:1])
	normismflux = plotismflux / np.median(plotismflux)
	trialplotnormismflux = normismflux + trial

	stellar_rest_wave = ism_wave + (ism_wave * (stellarrestvelocityoftarget/(3*10**5)))
	#print("stellar_rest_wave", stellar_rest_wave)
	bool_keep_stellar_rest = (stellar_rest_wave > plotrange[0]) & (stellar_rest_wave < plotrange[1])
	plotstellar_rest_wave = stellar_rest_wave[bool_keep_stellar_rest]
	plotstellar_rest_flux = flux[bool_keep_stellar_rest]
	normstellarrestflux = plotstellar_rest_flux / np.median(plotstellar_rest_flux)
	trialplotnormstellarrestflux = normstellarrestflux + trial
	
	#trial = trial + 0.01

	#plt.plot(plotbary_wave, trialplotnormbaryflux, alpha=1, label="Barycentric Frame {}".format(sp.target))
	#plt.plot(plotism_wave, trialplotnormismflux, alpha=1.0, label="ISM Frame {}".format(sp.target))
	#plt.plot(plotgeo_wave, trialplotnormgeoflux, alpha=1.0, label="Geocentric Frame {}".format(sp.target))
	#plt.plot(plotstellar_rest_wave, trialplotnormstellarrestflux, alpha=1, label="Stellar Rest Frame {}".format(sp.target))
	
	#plt.xlim(plotrange[0], plotrange[1])
	#ylim = [np.min(trialplotnormismflux), np.max(trialplotnormismflux)]
	#dynrange = ylim[1]-ylim[0]
	#plt.plot(labspec.wavelength, labspec.norm * dynrange/2 + 1, color='lime', alpha=1, linewidth = 3.0, label="C6H Lab Spectrum Bands 1,2")
	#plt.legend()
	#plt.show()
	
	totalnormfluxliststellar.append(normstellarrestflux)
	totalnormfluxlistism.append(normismflux)
	totalnormfluxlistgeo.append(normgeoflux)
	totalnormfluxlistbary.append(normbaryflux)

totalsumfluxliststellar = [sum(i) for i in zip(*totalnormfluxliststellar)]
normalizedtotalsumfluxliststellar = totalsumfluxliststellar / np.median(totalsumfluxliststellar)
#print(normalizedtotalsumfluxliststellar)
print(len(normalizedtotalsumfluxliststellar))
print(len(totalplotwave))

totalsumfluxlistism = [sum(i) for i in zip(*totalnormfluxlistism)]
normalizedtotalsumfluxlistism = totalsumfluxlistism / np.median(totalsumfluxlistism)
#Saving the Data:
#np.savetxt("High_Telluric_ISM", normalizedtotalsumfluxlistism)

totalsumfluxlistgeo = [sum(i) for i in zip(*totalnormfluxlistgeo)]
normalizedtotalsumfluxlistgeo = totalsumfluxlistgeo / np.median(totalsumfluxlistgeo)
#Saving the Data:
#np.savetxt("High_Telluric_GEO", normalizedtotalsumfluxlistgeo)

totalsumfluxlistbary = [sum(i) for i in zip(*totalnormfluxlistbary)]
normalizedtotalsumfluxlistbary = totalsumfluxlistbary / np.median(totalsumfluxlistbary)
#Saving the Data:
#np.savetxt("High_Telluric_BARY", normalizedtotalsumfluxlistbary)

plt.xlim(plotrange[0], plotrange[1])
plotdifference = np.max(normalizedtotalsumfluxlistism) - np.min(normalizedtotalsumfluxlistism)
print("plotdifference (voor ander soort plotje)", plotdifference)
ylim = [(np.min(normalizedtotalsumfluxlistism)-plotdifference/5), (np.max(normalizedtotalsumfluxlistism)+plotdifference/5)]
#ylim = [np.min(normalizedtotalsumfluxlistism), np.max(normalizedtotalsumfluxlistism)]
plt.ylim(ylim)
#plt.ylim(0.85, 2.5)
#plt.ylim(0.9, 2.2)
#plt.ylim(0.87, 1.5)
#plt.ylim(0.9, 1.3)
#plt.ylim(0.98, 1.13)
#plt.ylim(0.93, 1.04)
#plt.ylim(0.996, 1.010)
plt.ylim(0.995, 1.0252)
#plt.plot(totalplotwave, normalizedtotalsumfluxlistgeo, color='b', alpha=0.9, linewidth=2, label="Stacked Edibles Data Normalized Geocentric Frame")
#plt.plot(totalplotwave, normalizedtotalsumfluxlistbary, color='c', alpha=0.8, linewidth=2, label="Stacked Edibles Data Normalized Barycentric Frame")
#plt.plot(totalplotwave, normalizedtotalsumfluxlistism, color='k', alpha=1, linewidth=3.0, label="Stacked Edibles Data Normalized ISM Frame")
#plt.plot(totalplotwave, normalizedtotalsumfluxliststellar, color='darkviolet', alpha=1, linewidth=3.0, label="Stacked Edibles Data Normalized Stellar Rest Frame")



#"""
#Crafted Stacks in all Frames for Low_Extinction vs High Extinction Showal:
#1:
Low_Extinction_Stack_GEO_1 = np.loadtxt(fname = "Low_Extinction_Stack__GEO_1.txt")
Low_Extinction_Stack_BARY_1 = np.loadtxt(fname = "Low_Extinction_Stack_BARY_1.txt")
Low_Extinction_Stack_ISM_1 = np.loadtxt(fname = "Low_Extinction_Stack__ISM_1.txt")
plt.plot(totalplotwave, np.asarray(Low_Extinction_Stack_GEO_1)+0.018, color='b', alpha=0.6, linewidth=1.5)
plt.plot(totalplotwave, np.asarray(Low_Extinction_Stack_BARY_1)+0.018, color='c', alpha=0.4, linewidth=1.5)
plt.plot(totalplotwave, np.asarray(Low_Extinction_Stack_ISM_1)+0.018, color='darkkhaki', alpha=1, linewidth=3.0, label="Stacked E(B-V)<0.2 15 Targets in ISM Frame")
#2:
Low_Extinction_Stack_GEO_2 = np.loadtxt(fname = "Low_Extinction_Stack__GEO_2.txt")
Low_Extinction_Stack_BARY_2 = np.loadtxt(fname = "Low_Extinction_Stack_BARY_2.txt")
Low_Extinction_Stack_ISM_2 = np.loadtxt(fname = "Low_Extinction_Stack__ISM_2.txt")
plt.plot(totalplotwave, np.asarray(Low_Extinction_Stack_GEO_2)+0.015, color='b', alpha=0.6, linewidth=1.5)
plt.plot(totalplotwave, np.asarray(Low_Extinction_Stack_BARY_2)+0.015, color='c', alpha=0.4, linewidth=1.5)
plt.plot(totalplotwave, np.asarray(Low_Extinction_Stack_ISM_2)+0.015, color='gold', alpha=1, linewidth=3.0, label="Stacked E(B-V)<0.3 39 Targets in ISM Frame")
#3:
Low_Extinction_Stack_GEO_3 = np.loadtxt(fname = "Low_Extinction_Stack__GEO_3.txt")
Low_Extinction_Stack_BARY_3 = np.loadtxt(fname = "Low_Extinction_Stack_BARY_3.txt")
Low_Extinction_Stack_ISM_3 = np.loadtxt(fname = "Low_Extinction_Stack__ISM_3.txt")
plt.plot(totalplotwave, np.asarray(Low_Extinction_Stack_GEO_3)+0.012, color='b', alpha=0.6, linewidth=1.5)
plt.plot(totalplotwave, np.asarray(Low_Extinction_Stack_BARY_3)+0.012, color='c', alpha=0.4, linewidth=1.5)
plt.plot(totalplotwave, np.asarray(Low_Extinction_Stack_ISM_3)+0.012, color='orange', alpha=1, linewidth=3.0, label="Stacked E(B-V)<0.4 48 Targets in ISM Frame")
#4:
Low_Extinction_Stack_GEO_4 = np.loadtxt(fname = "Low_Extinction_Stack__GEO_4.txt")
Low_Extinction_Stack_BARY_4 = np.loadtxt(fname = "Low_Extinction_Stack_BARY_4.txt")
Low_Extinction_Stack_ISM_4 = np.loadtxt(fname = "Low_Extinction_Stack__ISM_4.txt")
plt.plot(totalplotwave, np.asarray(Low_Extinction_Stack_GEO_4)+0.009, color='b', alpha=0.6, linewidth=1.5)
plt.plot(totalplotwave, np.asarray(Low_Extinction_Stack_BARY_4)+0.009, color='c', alpha=0.4, linewidth=1.5)
plt.plot(totalplotwave, np.asarray(Low_Extinction_Stack_ISM_4)+0.009, color='red', alpha=1, linewidth=3.0, label="Stacked E(B-V)<0.6 59 Targets in ISM Frame")
#5:
Low_Extinction_Stack_GEO_5 = np.loadtxt(fname = "Low_Extinction_Stack__GEO_5.txt")
Low_Extinction_Stack_BARY_5 = np.loadtxt(fname = "Low_Extinction_Stack_BARY_5.txt")
Low_Extinction_Stack_ISM_5 = np.loadtxt(fname = "Low_Extinction_Stack__ISM_5.txt")
plt.plot(totalplotwave, np.asarray(Low_Extinction_Stack_GEO_5)+0.006, color='b', alpha=0.6, linewidth=1.5)
plt.plot(totalplotwave, np.asarray(Low_Extinction_Stack_BARY_5)+0.006, color='c', alpha=0.4, linewidth=1.5)
plt.plot(totalplotwave, np.asarray(Low_Extinction_Stack_ISM_5)+0.006, color='indianred', alpha=1, linewidth=3.0, label="Stacked E(B-V)<1.1 70 Targets in ISM Frame")
#6:
Low_Extinction_Stack_GEO_6 = np.loadtxt(fname = "Low_Extinction_Stack__GEO_6.txt")
Low_Extinction_Stack_BARY_6 = np.loadtxt(fname = "Low_Extinction_Stack_BARY_6.txt")
Low_Extinction_Stack_ISM_6 = np.loadtxt(fname = "Low_Extinction_Stack__ISM_6.txt")
plt.plot(totalplotwave, np.asarray(Low_Extinction_Stack_BARY_6)+0.003, color='c', alpha=0.4, linewidth=1.5)
plt.plot(totalplotwave, np.asarray(Low_Extinction_Stack_ISM_6)+0.003, color='darkred', alpha=1, linewidth=3.0, label="Stacked 0.5<E(B-V)<1.1 16 Targets in ISM Frame")
plt.plot(totalplotwave, np.asarray(Low_Extinction_Stack_GEO_6)+0.003, color='b', alpha=0.6, linewidth=1.5, label="Stacked Targets in Respective Geocentric Frames")

plt.axvline(5264.94, color='k', linewidth = 3.0)
plt.axvline(5267.10, color='k', linewidth = 3.0)
#plt.axvline(5264.4, color='darkred', linewidth=3.0, label="Stellar Lines")
#"""

"""
#Similar, but then for Band19:
#1:
Low_Extinction_Stack_GEO_19 = np.loadtxt(fname= "Low_Extinction_Stack_GEO_19.txt")
Low_Extinction_Stack_ISM_19 = np.loadtxt(fname= "Low_Extinction_Stack_ISM_19.txt")
plt.plot(totalplotwave, np.asarray(Low_Extinction_Stack_GEO_19)+0.006, color='b', alpha=0.6, linewidth=1.5)
plt.plot(totalplotwave, np.asarray(Low_Extinction_Stack_ISM_19)+0.006, color='orange', alpha=1, linewidth=2.5, label="Stacked E(B-V)<0.3 16 Targets in ISM Frame")
#2:
High_Extinction_Stack_GEO_19 = np.loadtxt(fname= "High_Extinction_Stack_GEO_19.txt")
High_Extinction_Stack_ISM_19 = np.loadtxt(fname= "High_Extinction_Stack_ISM_19.txt")
plt.plot(totalplotwave, np.asarray(High_Extinction_Stack_GEO_19)+0.003, color='b', alpha=0.6, linewidth=1.5)
plt.plot(totalplotwave, np.asarray(High_Extinction_Stack_ISM_19)+0.003, color='darkred', alpha=1, linewidth=2.5, label="Stacked 0.5<E(B-V)<1.1 17 Targets in ISM Frame")

plt.axvline(4743.95, color='k', linewidth = 3.0)
plt.axvline(4746.00, color='k', linewidth = 3.0)
"""

"""
#Now for the Telluric Test:
dynrange = ylim[1]-ylim[0]
#1
Low_T_Geo = np.loadtxt(fname = "Less_Telluric_GEO.txt")
Low_T_Bary = np.loadtxt(fname = "Less_Telluric_BARY.txt")
Low_T_ISM = np.loadtxt(fname = "Less_Telluric_ISM.txt")
plt.plot(totalplotwave, np.asarray(Low_T_Geo)+0.007, color='b', alpha=0.9, linewidth=2.0)
plt.plot(totalplotwave, np.asarray(Low_T_Bary)+0.007, color='c', alpha=0.8, linewidth=2.0)
plt.plot(totalplotwave, np.asarray(Low_T_ISM)+0.007, color='darkred', alpha=1, linewidth=3.0, label="Filtered Emission Stack 133 targets")
plt.plot(labspec.wavelength, (labspec.norm * dynrange/2 + 1)+0.007, color='lime', alpha=1, linewidth = 3.0)
#2
High_T_Geo = np.loadtxt(fname = "High_Telluric_GEO.txt")
High_T_Bary = np.loadtxt(fname = "High_Telluric_BARY.txt")
High_T_ISM = np.loadtxt(fname = "High_Telluric_ISM.txt")
plt.plot(totalplotwave, np.asarray(High_T_Geo)-0.008, color='b', alpha=0.9, linewidth=2.0)
plt.plot(totalplotwave, np.asarray(High_T_Bary)-0.008, color='c', alpha=0.8, linewidth=2.0)
plt.plot(totalplotwave, np.asarray(High_T_ISM)-0.008, color='orange', alpha=1, linewidth=3.0, label="High Emission Stack 74 targets")
plt.plot(labspec.wavelength, (labspec.norm * dynrange/2 + 1)-0.008, color='lime', alpha=1, linewidth = 3.0)

plt.axvline(5264.94, color='k', linewidth = 3.0)
plt.axvline(5267.10, color='k', linewidth = 3.0)
"""

dynrange = ylim[1]-ylim[0]
#Band1,2:
plt.plot(labspec.wavelength, labspec.norm * dynrange/2 +1, color='lime', alpha=1, linewidth = 3.0, label="C6H Spectrum Bands 1,2")
#Band19:
#plt.plot(labspec.wavelength, labspec.norm * dynrange/2 + 1, color='lime', alpha=1, linewidth = 3.0, label="C6H Lab Spectrum Band 19")
#Band18:
#plt.plot(labspec.wavelength, labspec.norm * dynrange/2 + 1, color='lime', alpha=1, linewidth = 3.0, label="C6H Lab Spectrum Band 18")
#C8H:
#plt.plot(labspec.wavelength, labspec.norm * dynrange/2 + 1.001, color='lime', alpha=1, linewidth = 3.0, label="C8H Spectrum Main 2 Bands")
#C10H:
#plt.plot(labspec.wavelength, labspec.norm * dynrange/2 + 1.003, color='lime', alpha=1, linewidth = 3.0, label="C10H Spectrum Main Band")

#9x Plots averaging all temperatures:
#plt.plot(labspec1.wavelength, labspec1.norm * dynrange/2 +1, color='lime', alpha=1, linewidth = 3.0, label="C6H Spectrum Bands 1,2, 2K")
#plt.plot(labspec2.wavelength, labspec2.norm * dynrange/2 +1, color='lime', alpha=1, linewidth = 3.0, label="C6H Spectrum Bands 1,2, 5K")
#plt.plot(labspec3.wavelength, labspec3.norm * dynrange/2 +1, color='lime', alpha=1, linewidth = 3.0, label="C6H Spectrum Bands 1,2, 10K")
#plt.plot(labspec4.wavelength, labspec4.norm * dynrange/2 +1, color='lime', alpha=1, linewidth = 3.0, label="C6H Spectrum Bands 1,2, 15K")
#plt.plot(labspec5.wavelength, labspec5.norm * dynrange/2 +1, color='lime', alpha=1, linewidth = 3.0, label="C6H Spectrum Bands 1,2, 20K")
#plt.plot(labspec6.wavelength, labspec6.norm * dynrange/2 +1, color='lime', alpha=1, linewidth = 3.0, label="C6H Spectrum Bands 1,2, 25K")
#plt.plot(labspec7.wavelength, labspec7.norm * dynrange/2 +1, color='lime', alpha=1, linewidth = 3.0, label="C6H Spectrum Bands 1,2, 50K")
#plt.plot(labspec8.wavelength, labspec8.norm * dynrange/2 +1, color='lime', alpha=1, linewidth = 3.0, label="C6H Spectrum Bands 1,2, 200K")
#plt.plot(labspec9.wavelength, labspec9.norm * dynrange/2 +1, color='lime', alpha=1, linewidth = 3.0, label="C6H Spectrum Bands 1,2, 2000K")


#DIBS and Stellar Features:
#Band1,2 Range:
plt.axvline(5236.27, color='m', linewidth = 3.0)
plt.axvline(5245.44, color='m', linewidth = 3.0)
plt.axvline(5247.39, color='m', linewidth = 3.0)
plt.axvline(5251.72, color='m', linewidth = 3.0)
plt.axvline(5257.44, color='m', linewidth = 3.0)
plt.axvline(5262.44, color='m', linewidth = 3.0, label="Haoyu DIB Catalog Known DIBS")
plt.axvline(5285.53, color='m', linewidth = 3.0)
plt.axvline(5298.05, color='m', linewidth = 3.0)
plt.axvline(5299.37, color='m', linewidth = 3.0)
#Bands345678910111213141516171819 Range:
plt.axvline(5217.85, color='m', linewidth = 3.0)
plt.axvline(5176.00, color='m', linewidth = 3.0)
plt.axvline(5170.49, color='m', linewidth = 3.0)
plt.axvline(5117.62, color='m', linewidth = 3.0)
plt.axvline(5110.77, color='m', linewidth = 3.0)
plt.axvline(5100.96, color='m', linewidth = 3.0)
plt.axvline(5092.17, color='m', linewidth = 3.0)
plt.axvline(5074.48, color='m', linewidth = 3.0)
plt.axvline(5061.50, color='m', linewidth = 3.0)
plt.axvline(5054.85, color='m', linewidth = 3.0)
plt.axvline(4817.62, color='m', linewidth = 3.0)
plt.axvline(4734.77, color='m', linewidth = 3.0)
plt.axvline(4761, color='m', linewidth = 3.0)
plt.axvline(4662.44, color='m', linewidth = 3.0)
#C8H range:
plt.axvline(6259.65, color='m', linewidth = 3.0)
#Stellar Lines
#plt.axvline(5282, color='darkred', linewidth = 3.0, label="Stellar Lines from DIB Surveys")
plt.axvline(5276, color='darkred', linewidth = 3.0)
plt.axvline(5272.5, color='darkred', linewidth = 3.0)
plt.axvline(5242.5, color='darkred', linewidth = 3.0)
plt.axvline(5235.5, color='darkred', linewidth = 3.0)
#plt.axvline(5159.7, color='goldenrod', linewidth = 4.0, label="Strong Stellar Line O-II for Stellar Rest Frame")
#plt.axvline(5155.7, color='darkgoldenrod', linewidth = 4.0, label="Strong Stellar Line Fe-III for Stellar Rest Frame")
#plt.axvline(5015.678, color='grey', linewidth = 3.0, label="Stronger Stellar Line HeI for Stellar Rest Frame")
plt.axvline(4751.5, color='darkred', linewidth = 3.0)
plt.axvline(4741.5, color='darkred', linewidth = 3.0)
plt.axvline(4843.0, color='darkred', linewidth = 3.0)
plt.xlabel("Air Wavelength in Angstrom")
plt.ylabel("Normalized Total Flux")
plt.legend(bbox_to_anchor=(0, 1.06, 0, 0.1), loc = 'upper left', prop={"size":5.95})
#plt.legend(framealpha=0.8, loc='upper left', ncol=4, fontsize= 'x-small')
plt.grid()
plt.show()

