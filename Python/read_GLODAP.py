import pandas as pd
import time
import numpy as np
import PyCO2SYS as pyco2
import itertools

# Main script calling to
import os
import time

# Store relative path to this script
# script_dir = os.path.dirname(os.path.realpath(__file__))
# script_dir = os.getwd()
input_files_dir = '/Users/rpr061/Documents/'
files_path_remote = ('https://www.ncei.noaa.gov'
                     '/data/oceans/ncei/ocads/data/0237935/')
output_files_dir = '/Users/rpr061/Documents/CMEMS_INSTAC_PRODUCT/'

GLODAP_files = ['GLODAPv2.2021_Merged_Master_File.csv', 'EXPOCODES.txt', 'Dataset_DOIs.txt']
info_file = 'GLODAPv22021CMEMS.tsv'

# Create and store output directory
if not os.path.isdir(output_files_dir):
    os.mkdir(output_files_dir)
# Check if input data files are there

if not all([os.path.isfile(os.path.join(input_files_dir, f)) for f in GLODAP_files]):
    input_files_dir = files_path_remote

headerlines = 0
separator = ','
ddtype = None
sourcefile = [i for i in GLODAP_files if 'GLODAP' in i][0]
doifile = [i for i in GLODAP_files if 'DOI' in i][0]
expocodefile = [i for i in GLODAP_files if 'EXPOCODES' in i][0]

# Read files
tempdf = pd.read_csv(os.path.join(input_files_dir, sourcefile), sep=separator,
                     skiprows=headerlines, dtype=ddtype, na_values=-9999,
                     on_bad_lines='skip')
dois = pd.read_csv(os.path.join(input_files_dir, doifile), sep='\t', header=None,
                   names=['G2cruise', 'DOI'], dtype=ddtype, on_bad_lines='skip',
                   encoding='utf_16_le')
expocodes = pd.read_csv(os.path.join(input_files_dir, expocodefile), sep='\t',
                        header=None, names=['G2cruise', 'EXPOCODE'], dtype=ddtype,
                        on_bad_lines='skip')
GLODAP_info = pd.read_csv(os.path.join('/Users/rpr061/Downloads', info_file), sep='\t',
                          skiprows=headerlines, dtype=ddtype, na_values=-9999,
                          on_bad_lines='skip')

# Some parameters / unchanging values / dictionaries
nominal_depths = tuple(itertools.chain(range(0, 105, 5), range(100, 1025, 25), range(1000, 10100, 100)))
input_vars = ('G2bottomdepth', 'G2pressure', 'G2temperature', 'G2salinity', 'G2oxygen', 'G2nitrate',
              'G2nitrite', 'G2phosphate', 'G2silicate', 'G2phtsinsitutp', 'G2phts25p0', 'G2tco2',
              'G2talk', 'G2doc', 'G2don', 'G2tdn', 'G2chla')
input_flag_vars = ('G2bottomdepthf', 'G2pressure', 'G2temperaturef', 'G2salinityf', 'G2oxygenf', 'G2nitratef',
                   'G2nitritef', 'G2phosphatef', 'G2silicatef', 'G2phtsinsitutpf', 'G2phts25p0f', 'G2tco2f',
                   'G2talkf', 'G2docf', 'G2donf', 'G2tdnf', 'G2chlaf')
SDN_var_names = ('BATH', '',
                 'TEMP', 'PSAL', 'DOX2', 'NTAW',
                 'NTIW', 'PHOW', 'SLCW', 'PHPH', 'PH25', 'TICW',
                 'ALKW', 'CORG', 'NODW', 'NT1D', 'CPHL')
units = ('m', '',
         'degrees_C', '0.001', 'µmol kg-1', 'µmol kg-1',
         'µmol kg-1', 'µmol kg-1', 'µmol kg-1', '1', '1', 'µmol kg-1',
         'µmol kg-1', 'µmol kg-1', 'µmol kg-1', 'µmol kg-1', 'mg m-3')
long_name = ('Bathymetric depth', '',
             'Sea temperature', 'Practical salinity',
             'Dissolved oxygen', 'Nitrate (NO3-N)',
             'Nitrite (NO2-N)', 'Phosphate (PO4-P)', 'Silicate (SIO4-SI)',
             'Ph', 'Ph at 25 °C and 0 dbar', 'Dissolved inorganic carbon',
             'Total alkalinity', 'Dissolved organic carbon', 'Dissolved organic nitrogen',
             'Total dissolved nitrogen', 'Chlorophyll-a')
CF_standard_name = ('sea_floor_depth_below_sea_surface', '',
                    'sea_water_temperature',
                    'sea_water_practical_salinity', 'moles_of_oxygen_per_unit_mass_in_sea_water',
                    'moles_of_nitrate_per_unit_mass_in_sea_water', 'moles_of_nitrite_per_unit_mass_in_sea_water',
                    'moles_of_phosphate_per_unit_mass_in_sea_water', 'moles_of_silicate_per_unit_mass_in_sea_water',
                    'sea_water_ph_reported_on_total_scale', ' ',
                    ' ', 'sea_water_alkalinity_per_unit_mass',
                    'moles_of_dissolved_organic_carbon_per_unit_mass_in_sea_water', ' ',
                    'moles_of_dissolved_total_nitrogen_per_unit_mass_in_sea_water',
                    'mass_concentration_of_chlorophyll_a_in_sea_water')
variables_dict = {}
variables_dict['flag'] = dict(zip(input_vars, input_flag_vars))
variables_dict['SDN'] = dict(zip(input_vars, input_flag_vars))
variables_dict['unit'] = dict(zip(input_vars, units))
variables_dict['long'] = dict(zip(input_vars, long_name))
variables_dict['CF'] = dict(zip(input_vars, CF_standard_name))

# Calculate the nominal depths series
nominal_depths_array = np.array(nominal_depths)
actual_depths_array = np.array(tempdf['G2depth'])
depths_as_nominal = nominal_depths_array[
    abs(actual_depths_array[None, :] - nominal_depths_array[:, None]).argmin(axis=0)]
tempdf['G2depthnominal'] = depths_as_nominal

# Change NaN hour/minute to 0
tempdf.loc[np.isnan(tempdf['G2hour']), 'G2hour'] = 0
tempdf.loc[np.isnan(tempdf['G2minute']), 'G2minute'] = 0

### QC flags
# If temperature and bottom depth exist, they're assumed good (flag 1)
# If value not recorded/measured, then flag 0
tempdf['G2bottomdepthf'] = np.ones(tempdf['G2bottomdepth'].__len__())
tempdf.loc[np.isnan(tempdf['G2bottomdepth']), 'G2bottomdepthf'] = 9
tempdf['G2temperaturef'] = np.ones(tempdf['G2temperature'].__len__())
tempdf.loc[np.isnan(tempdf['G2temperature']), 'G2temperaturef'] = 9

# Remap flags to OceanSITES NaN -> 0; 0->8; 1->9; 2->1; 3->2; 4->4;
# % 5->9; 6->8; 7->0; 8->0; 9->9 as int8
GLODAP_flags = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
OceanSITES_flags = [8, 9, 1, 2, 4, 9, 8, 0, 0, 9]
G2OS_flag_dict = dict(zip(GLODAP_flags, OceanSITES_flags))

for fv in variables_dict['flag']:
    if fv:
        print(fv)
        tempdf.loc[np.isnan(tempdf[fv]), [fv]] = 9
        for fvv in G2OS_flag_dict.keys():
            print(fvv, G2OS_flag_dict[fvv])
            tempdf.loc[tempdf[fv] == fvv, [fv]] = G2OS_flag_dict[fvv]

### Extract from each platform
# create the platform codes unique list, merging Callsign and Name (if callsign missing)
GLODAP_info['PlatformCode']=GLODAP_info['CallSign_WMO']
GLODAP_info.loc[GLODAP_info['PlatformCode'].isna(),'PlatformCode']=\
    GLODAP_info.loc[GLODAP_info['PlatformCode'].isna(),'Name']
GLODAP_info['PlatformCode']=GLODAP_info['PlatformCode'].str.replace(' ','')
unique_platform_codes=GLODAP_info['PlatformCode'].unique()



# loop through the platforms (the basis of INSTAC files)
for pc in unique_platform_codes:
    cruises=expocodes.loc[expocodes['EXPOCODE'].str.contains('06AQ')]


MeteorCruises = expocodes.loc[expocodes['EXPOCODE'].str.contains('06AQ')]
MeteorData = tempdf.loc[tempdf['G2cruise'].isin(MeteorCruises['G2cruise'])]

### The matlab way of creating the depth dimension was oh boy SO WRONGG


###############
# Rename G2fco2 to G2fco2_20_0 (in GLODAP, it's given at 20 dg, 0dbar).
tempdf.rename(columns={'G2fco2': 'G2fco2_20_0'}, inplace=True)

# Subset for surface and reset indices
# Upper 10 m, with fCO2 measurement (measured and calculated; # Flag 2 for Good,
# MEASURED. Flag 0 is for calculated)
surfilt = (tempdf["G2depth"] <= 10.) & (~pd.isna(tempdf["G2fco2_20_0"])) & (
        (tempdf["G2fco2f"] == 2) | (tempdf["G2fco2f"] == 0))
dfsurf = tempdf[surfilt]
dfsurf.reset_index(drop=True, inplace=True)

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
dfsurfgood['DATEVECTOR1'] = pd.to_datetime(tempdtframe, utc=True)
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
dfsurf.reset_index(drop=True, inplace=True)

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
    par1_type=1,  # The first parameter is of type "1". meaning "alkalinity"
    par1=dfsurfgood['G2talk'],  # value of the first parameter
    par2_type=5,  # The second parameter is of type "5", meaning "fCO2"
    par2=dfsurfgood['G2fco2_20_0'],  # value of the second parameter
    salinity=dfsurfgood['G2salinity'],  # Salinity of the sample
    temperature=20,  # T at input conditions
    temperature_out=dfsurfgood['G2temperature'],  # T at output conditions
    pressure=0,  # Pressure    at input conditions
    pressure_out=dfsurfgood['G2pressure'],  # Pressure at output conditions
    total_silicate=dfsurfgood['G2silicate'],  # Silicate in sample [umol/kg]
    total_phosphate=dfsurfgood['G2phosphate'],  # Phosphate in sample [umol/kg]
    opt_pH_scale=1,  # pH scale of the input ("1" means "Total Scale")
    opt_k_carbonic=10,  # Choice of H2CO3 and HCO3- dissociation constants K1
    # and K2 ("10" means "Lueker 2000")
    opt_k_bisulfate=1,  # Choice of HSO4- dissociation constant KSO4 ("1"
    # means "Dickson")
    opt_total_borate=1,  # Choice of boron:sal ("1" means "Uppstrom")
)
start_time = time.time()
print('CO2SYS Conditions have been defined!')
results = pyco2.sys(**kwargs)
print("--- %s seconds ---" % (time.time() - start_time))
dfsurfgood['G2fco2'] = results['fCO2_out']

# Assign calculation method (based on what data is available) 0: measured,
# 1:f(Alk,DIC), 2ALK,pH, 3 DIC ph
dfsurfgood[vardict['fco2wc']] = 9  # Data not available
# dfsurfgood[vardict['dicc']]=0
# dfsurfgood[vardict['alkc']]=0
# dfsurfgood[vardict['phc']]=0

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
    columns={'EXPOCODE': vardict['id'], 'DOI': vardict['doi'],
             'G2latitude': vardict['lat'], 'G2longitude': vardict['lon'],
             'G2depth': vardict['dep'], 'G2temperature': vardict['temp'],
             'G2salinity': vardict['sal'], 'G2salinityf': vardict['salf'],
             'G2tco2': vardict['dic'], 'G2tco2qc': vardict['dicf'],
             'G2talk': vardict['alk'], 'G2talkqc': vardict['alkf'],
             'G2phtsinsitutp': vardict['ph'], 'G2phtsqc': vardict['phf'],
             'G2fco2': vardict['fco2w'], 'G2fco2f': vardict['fco2wf']},
    inplace=True)

print(tempdf.shape, dfsurf.shape, dfsurfgood.shape)

# Add source (SOCAT, GLODAP, ARGO, etc...)
dfsurfgood['SOURCE'] = source

# Rename and reset indices
printdf = dfsurfgood
printdf.reset_index(drop=True, inplace=True)

print('GLODAP frame size is ')
print(printdf.shape)
