import pandas as pd
import time
import numpy as np
import PyCO2SYS as pyco2

# Main script calling to
import os
import time


# Store relative path to this script
#script_dir = os.path.dirname(os.path.realpath(__file__))
#script_dir = os.getwd()
input_files_dir = '/Users/rpr061/Documents/'
files_path_remote = ('https://www.ncei.noaa.gov'
                   '/data/oceans/ncei/ocads/data/0237935/')
output_files_dir = '/Users/rpr061/Documents/CMEMS_INSTAC_PRODUCT/'

GLODAP_files=['GLODAPv2.2021_Merged_Master_File.csv','EXPOCODES.txt','Dataset_DOIs.txt']

# Create and store output directory
if not os.path.isdir(output_files_dir):
    os.mkdir(output_files_dir)
# Check if input data files are there

if not all([os.path.isfile(os.path.join(input_files_dir,f)) for f in GLODAP_files]):
    input_files_dir  = files_path_remote

headerlines = 0
separator = ','
ddtype = None
sourcefile = [i for i in GLODAP_files if 'GLODAP' in i][0]
doifile = [i for i in GLODAP_files if 'DOI' in i][0]
expocodefile = [i for i in GLODAP_files if 'EXPOCODES' in i][0]

# Read files
tempdf = pd.read_csv(os.path.join(input_files_dir, sourcefile), sep = separator,
                     skiprows = headerlines, dtype = ddtype, na_values = -9999,
                     on_bad_lines = 'skip')
dois = pd.read_csv(os.path.join(input_files_dir, doifile), sep ='\t', header = None,
                   names = ['G2cruise', 'DOI'], dtype = ddtype, on_bad_lines = 'skip',
                   encoding = 'utf_16_le')
expocodes = pd.read_csv(os.path.join(input_files_dir, expocodefile), sep ='\t',
                        header = None, names = ['G2cruise', 'EXPOCODE'], dtype = ddtype,
                        on_bad_lines = 'skip')

# Identify the nominal depths
tempdf['UNICAST'] = tempdf.set_index(['G2cruise', 'G2station',
    'G2cast']).index.factorize()[0]
alldepths=tempdf['G2depth']
allcasts=tempdf['UNICAST']

alldepths_onemeter=alldepths.round()

nominal_depths=tuple(itertools.chain(range(0,105,5),range(100,1025,25),range(1000,10100,100)))

### Extract from each platform
MeteorCruises=expocodes.loc[expocodes['EXPOCODE'].str.contains('06AQ')]
MeteorData=tempdf.loc[tempdf['G2cruise'].isin(MeteorCruises['G2cruise'])]


### The matlab way of creating the depth dimension was oh boy SO WRONGG



###############
# Rename G2fco2 to G2fco2_20_0 (in GLODAP, it's given at 20 dg, 0dbar).
tempdf.rename(columns = {'G2fco2': 'G2fco2_20_0'}, inplace = True)

# Subset for surface and reset indices
# Upper 10 m, with fCO2 measurement (measured and calculated; # Flag 2 for Good,
# MEASURED. Flag 0 is for calculated)
surfilt = (tempdf["G2depth"] <= 10.) & (~pd.isna(tempdf["G2fco2_20_0"])) & (
        (tempdf["G2fco2f"] == 2) | (tempdf["G2fco2f"] == 0))
dfsurf = tempdf[surfilt]
dfsurf.reset_index(drop = True, inplace = True)

# Filter for only the uppermost measurement at each unique cast (if, e.g.
# samples at 2 and 10 m)
dfsurf['UNICAST'] = dfsurf.set_index(['G2cruise', 'G2station',
    'G2cast']).index.factorize()[0] + 1
surfacemostind = []
for x in np.unique(dfsurf['UNICAST']):
    surfacemostind.append(dfsurf['G2depth'].iloc[np.where(dfsurf[
        'UNICAST'] == x)].idxmin())
dfsurfgood = dfsurf.iloc[surfacemostind].copy()

# Change to UNIX date format and create Python datevector (to work internally)
# Some hours and minutes are NA: change to 0.0
dfsurfgood['G2hour'].iloc[np.where(dfsurfgood['G2hour'].isna())] = 0.0
dfsurfgood['G2minute'].iloc[np.where(dfsurfgood['G2minute'].isna())] = 0.0
tempdtframe = pd.DataFrame(
    {'year': dfsurfgood['G2year'], 'month': dfsurfgood['G2month'],
    'day': dfsurfgood['G2day'], 'hour': dfsurfgood['G2hour'],
    'minute': dfsurfgood['G2minute']})
dfsurfgood['DATEVECTOR1'] = pd.to_datetime(tempdtframe, utc = True)
dfsurfgood[vardict['unixd']] = dfsurfgood['DATEVECTOR1'].astype(
    'int64') // 10 ** 9
dfsurfgood[vardict['datevec']] = dfsurfgood['DATEVECTOR1'].dt.strftime(
    '%Y-%m-%dT%H:%M:%SZ')

# Subset location and time
morefilt = (dfsurfgood['G2latitude'] >= minlat) & (
    dfsurfgood['G2latitude'] <= maxlat) & (
    dfsurfgood['G2longitude'] >= minlon) & (
    dfsurfgood['G2longitude'] <= maxlon) & (
    dfsurfgood['DATEVECTOR1'] >= pd.to_datetime(mindate)) & (
    dfsurfgood['DATEVECTOR1'] <= pd.to_datetime(maxdate))

dfsurfgood = dfsurfgood[morefilt]
dfsurf.reset_index(drop = True, inplace = True)

# Assign Expocodes and DOIs
# Transform expocode dataframes into dictionaries (easier to lookup)
expocodesdict = expocodes.set_index('G2cruise')['EXPOCODE'].to_dict()
doisdict = dois.set_index('G2cruise')['DOI'].to_dict()

for cruises in dfsurfgood['G2cruise'].unique():
    dfsurfgood.loc[dfsurfgood.G2cruise == cruises,
     'EXPOCODE'] = expocodesdict[cruises]
    dfsurfgood.loc[dfsurfgood.G2cruise == cruises,
     'DOI'] = doisdict[cruises].rsplit('.org/', 1)[1]

### CARBON STUFF
### Calculate fCO2 in situ (in GLODAP it's at 20 C, 0dbar)
# Define input and output conditions
kwargs = dict(
    par1_type = 1,  # The first parameter is of type "1". meaning "alkalinity"
    par1 = dfsurfgood['G2talk'],  # value of the first parameter
    par2_type = 5,  # The second parameter is of type "5", meaning "fCO2"
    par2 = dfsurfgood['G2fco2_20_0'],  # value of the second parameter
    salinity = dfsurfgood['G2salinity'],  # Salinity of the sample
    temperature = 20,  # T at input conditions
    temperature_out = dfsurfgood['G2temperature'],  # T at output conditions
    pressure = 0,  # Pressure    at input conditions
    pressure_out = dfsurfgood['G2pressure'],  # Pressure at output conditions
    total_silicate = dfsurfgood['G2silicate'],  # Silicate in sample [umol/kg]
    total_phosphate = dfsurfgood['G2phosphate'], # Phosphate in sample [umol/kg]
    opt_pH_scale = 1,  # pH scale of the input ("1" means "Total Scale")
    opt_k_carbonic = 10,  # Choice of H2CO3 and HCO3- dissociation constants K1
    # and K2 ("10" means "Lueker 2000")
    opt_k_bisulfate = 1,  # Choice of HSO4- dissociation constant KSO4 ("1"
    #means "Dickson")
    opt_total_borate = 1,  # Choice of boron:sal ("1" means "Uppstrom")
)
start_time = time.time()
print('CO2SYS Conditions have been defined!')
results = pyco2.sys(**kwargs)
print("--- %s seconds ---" % (time.time() - start_time))
dfsurfgood['G2fco2'] = results['fCO2_out']

# Assign calculation method (based on what data is available) 0: measured,
# 1:f(Alk,DIC), 2ALK,pH, 3 DIC ph
dfsurfgood[vardict['fco2wc']] = 9 # Data not available
#dfsurfgood[vardict['dicc']]=0
#dfsurfgood[vardict['alkc']]=0
#dfsurfgood[vardict['phc']]=0

calcindex = dfsurfgood.index[dfsurfgood[vardict['fco2wc']] == 0]
# '0' means calculated
for ind in calcindex:
    if (dfsurfgood['G2fco2f'][ind] == 2):
        dfsurfgood[vardict['fco2wc']][ind] = 0
    elif (dfsurfgood['G2tco2f'][ind] == 1) & (dfsurfgood['G2talkf'][ind] == 1):
                                              # Check if it should be 2 instead!
        dfsurfgood[vardict['fco2wc']][ind] = 1
    elif (dfsurfgood['G2phtsinsitutpf'][ind] == 1) & (
        dfsurfgood['G2talkf'][ind] == 1):
        dfsurfgood[vardict['fco2wc']][ind] = 2
    elif (dfsurfgood['G2tco2f'][ind] == 1) & (
        dfsurfgood['G2phtsinsitutpf'][ind] == 1):
        dfsurfgood[vardict['fco2wc']][ind] = 3


# Rename columns
dfsurfgood.rename(
    columns = {'EXPOCODE': vardict['id'], 'DOI': vardict['doi'],
      'G2latitude': vardict['lat'], 'G2longitude': vardict['lon'],
      'G2depth': vardict['dep'], 'G2temperature': vardict['temp'],
      'G2salinity': vardict['sal'], 'G2salinityf': vardict['salf'],
      'G2tco2': vardict['dic'], 'G2tco2qc': vardict['dicf'],
      'G2talk': vardict['alk'], 'G2talkqc': vardict['alkf'],
      'G2phtsinsitutp': vardict['ph'], 'G2phtsqc': vardict['phf'],
      'G2fco2': vardict['fco2w'], 'G2fco2f': vardict['fco2wf']},
    inplace = True)

print(tempdf.shape, dfsurf.shape, dfsurfgood.shape)

# Add source (SOCAT, GLODAP, ARGO, etc...)
dfsurfgood['SOURCE'] = source

# Rename and reset indices
printdf = dfsurfgood
printdf.reset_index(drop = True, inplace = True)

print('GLODAP frame size is ')
print(printdf.shape)