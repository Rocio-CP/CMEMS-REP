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

all_glodap_input_variables=['G2cruise','G2region','G2station','G2cast',
                            'G2year','G2month','G2day','G2hour','G2minute',
                            'G2latitude','G2longitude','G2bottomdepth','G2maxsampdepth',
                            'G2bottle','G2pressure','G2depth','G2temperature','G2theta',
                            'G2salinity','G2salinityf','G2salinityqc',
                            'G2sigma0','G2sigma1','G2sigma2','G2sigma3','G2sigma4','G2gamma',
                            'G2oxygen','G2oxygenf','G2oxygenqc','G2aou','G2aouf',
                            'G2nitrate','G2nitratef','G2nitrateqc','G2nitrite','G2nitritef',
                            'G2silicate','G2silicatef','G2silicateqc',
                            'G2phosphate','G2phosphatef','G2phosphateqc',
                            'G2tco2','G2tco2f','G2tco2qc','G2talk','G2talkf','G2talkqc',
                            'G2fco2','G2fco2f','G2fco2temp',
                            'G2phts25p0','G2phts25p0f','G2phtsinsitutp','G2phtsinsitutpf','G2phtsqc',
                            'G2cfc11','G2pcfc11','G2cfc11f','G2cfc11qc',
                            'G2cfc12','G2pcfc12','G2cfc12f','G2cfc12qc',
                            'G2cfc113','G2pcfc113','G2cfc113f','G2cfc113qc',
                            'G2ccl4','G2pccl4','G2ccl4f','G2ccl4qc',
                            'G2sf6','G2psf6','G2sf6f',
                            'G2c13','G2c13f','G2c13qc','G2c14','G2c14f','G2c14err',
                            'G2h3','G2h3f','G2h3err','G2he3','G2he3f','G2he3err',
                            'G2he','G2hef','G2heerr','G2neon','G2neonf','G2neonerr',
                            'G2o18','G2o18f','G2toc','G2tocf','G2doc','G2docf',
                            'G2don','G2donf','G2tdn','G2tdnf','G2chla','G2chlaf']
# pd.Int16Dtype() is the type of pandas integer that allows for na
all_glodap_input_variables_dtype=['UInt64']*all_glodap_input_variables.__len__()
for n,gi in enumerate(all_glodap_input_variables):
    if not gi.endswith('f') and 'qc' not in gi and gi not in ['G2cruise','G2region','G2cast','G2year','G2month','G2day','G2hour','G2minute','G2bottle']:
           all_glodap_input_variables_dtype[n]='float'
dtype_dict=dict(zip(all_glodap_input_variables,all_glodap_input_variables_dtype))
from datetime import datetime
custom_date_parser = lambda x: datetime.strptime(x, "%Y %m %d %H:%M")

# Read files
# parse_dates does not work because some hours and minutes are nan
tempdf = pd.read_csv(os.path.join(input_files_dir, sourcefile), sep=separator,
                     skiprows=headerlines, na_values=-9999,
                     dtype=dtype_dict,
                     on_bad_lines='skip')
# Change NaN hour/minute to 0
tempdf.loc[np.isnan(tempdf['G2hour']), 'G2hour'] = 0
tempdf.loc[np.isnan(tempdf['G2minute']), 'G2minute'] = 0


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


# Rename G2fco2 to G2fco2_20_0 (in GLODAP, it's given at 20 dg, 0dbar).
tempdf.rename(columns={'G2fco2': 'G2fco2_20_0'}, inplace=True)

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

### Assign Expocodes and DOIs
# Transform expocode dataframes into dictionaries (easier to lookup)
expocodesdict = expocodes.set_index('G2cruise')['EXPOCODE'].to_dict()
doisdict = dois.set_index('G2cruise')['DOI'].to_dict()

for cruise in tempdf['G2cruise'].unique():
    tempdf.loc[tempdf.G2cruise == cruise,
                   'EXPOCODE'] = expocodesdict[cruise]
    tempdf.loc[tempdf.G2cruise == cruise,
                   'DOI'] = doisdict[cruise].rsplit('.org/', 1)[1]


### Extract from each platform
# create the platform codes unique list, merging Callsign and Name (if callsign missing)
GLODAP_info['PlatformCode']=GLODAP_info['CallSign_WMO']
GLODAP_info.loc[GLODAP_info['PlatformCode'].isna(),'PlatformCode']=\
    GLODAP_info.loc[GLODAP_info['PlatformCode'].isna(),'Name']
GLODAP_info['PlatformCode']=GLODAP_info['PlatformCode'].str.replace(' ','')
unique_platform_codes=GLODAP_info['PlatformCode'].unique()



# loop through the platforms (the basis of INSTAC files)

filter_variables=['G2year','G2month','G2day','G2hour','G2minute',
                   'G2latitude','G2longitude', 'G2depthnominal',
                   *input_vars]

for pc in unique_platform_codes:
    # Create boolean filters for
    filter_expocodes=GLODAP_info['PlatformCode']==pc
    current_expocodes=GLODAP_info.loc[filter_expocodes,'EXPOCODE']
    print(pc, current_expocodes)
    filter_expocodes_data=tempdf['EXPOCODE'].isin(current_expocodes)

    # Platform-specific global attributes
    # Platform attributes
    platform_code= pc # only alphanumeric characters
    platform_name= GLODAP_info.loc[filter_expocodes, 'Name'].unique().item()
    wmo=GLODAP_info.loc[filter_expocodes, 'CallSign_WMO'].unique().item()
    ices=GLODAP_info.loc[filter_expocodes, 'ICEScode'].unique().item()
    if not wmo : wmo= ' '
    if not ices : ices= ' '
    platform_category = str(GLODAP_info.loc[filter_expocodes, 'PlatformType'].unique().item())
    if platform_category == '31' :
        source = 'research vessel'
        output_file_full_name='VESSEL/GL_PR_BO_'+platform_code+'-GLODAPv22022.nc'
    elif platform_category == '21' :
        source = 'propelled manned submersible'
        output_file_full_name='ETC/GL_PR_BO_'+platform_code+'-GLODAPv22022.nc'
    elif platform_category == '62' :
        source = 'aeroplane'
        output_file_full_name='ETC/GL_PR_BO_'+platform_code+'-GLODAPv22022.nc'

    # Institution attributes
    institution_name=GLODAP_info.loc[filter_expocodes, 'Institution'].unique()
    edmo=GLODAP_info.loc[filter_expocodes, 'EDMO'].unique()
    pi_institution_name=GLODAP_info.loc[filter_expocodes, 'PI_Institution'].unique()
    pi_edmo=GLODAP_info.loc[filter_expocodes, 'PI_EDMO'].unique()
    all_institution_name="/".join(np.unique([*institution_name, *pi_institution_name]))
    all_edmo=",".join(np.unique([*edmo, *pi_edmo]))

    # PI attributes
    pi=GLODAP_info.loc[filter_expocodes, 'PI'].unique()
    pi=";".join(pi)

    # Extract sub-dataframe
    current_dataframe = tempdf.loc[filter_expocodes_data,filter_variables]

    # Create date variable
    tempdtframe = pd.DataFrame(
        {'year': current_dataframe['G2year'], 'month': current_dataframe['G2month'],
         'day': current_dataframe['G2day'], 'hour': current_dataframe['G2hour'],
         'minute': current_dataframe['G2minute']})
    #tempdf['DATEVECTOR1'] = pd.to_datetime(tempdtframe, utc=True)

    # Create days-from-1950 numeric date
    datetime_diff_1950 = tempdtframe - pd.Timestamp('1950-01-01T00:00:00', tz='UTC')
    datetime_numeric_1950 = datetime_diff_1950.dt.total_seconds() / 86400
    time_variable = datetime_numeric_1950[0:datetime_numeric_1950.idxmax() + 1]

    # Order by increasing depth values



    # Average at same depths and casts

    break


### The matlab way of creating the depth dimension was oh boy SO WRONGG


############
#
# # Subset for surface and reset indices
# # Upper 10 m, with fCO2 measurement (measured and calculated; # Flag 2 for Good,
# # MEASURED. Flag 0 is for calculated)
# surfilt = (tempdf["G2depth"] <= 10.) & (~pd.isna(tempdf["G2fco2_20_0"])) & (
#         (tempdf["G2fco2f"] == 2) | (tempdf["G2fco2f"] == 0))
# dfsurf = tempdf[surfilt]
# dfsurf.reset_index(drop=True, inplace=True)
#
# # Filter for only the uppermost measurement at each unique cast (if, e.g.
# # samples at 2 and 10 m)
# dfsurf['UNICAST'] = dfsurf.set_index(['G2cruise', 'G2station',
#                                       'G2cast']).index.factorize()[0] + 1
# surfacemostind = []
# for x in np.unique(dfsurf['UNICAST']):
#     surfacemostind.append(dfsurf['G2depth'].iloc[np.where(dfsurf[
#                                                               'UNICAST'] == x)].idxmin())
# dfsurfgood = dfsurf.iloc[surfacemostind].copy()
#
# # Change to UNIX date format and create Python datevector (to work internally)
# # Some hours and minutes are NA: change to 0.0
# dfsurfgood['G2hour'].iloc[np.where(dfsurfgood['G2hour'].isna())] = 0.0
# dfsurfgood['G2minute'].iloc[np.where(dfsurfgood['G2minute'].isna())] = 0.0
# tempdtframe = pd.DataFrame(
#     {'year': dfsurfgood['G2year'], 'month': dfsurfgood['G2month'],
#      'day': dfsurfgood['G2day'], 'hour': dfsurfgood['G2hour'],
#      'minute': dfsurfgood['G2minute']})
# dfsurfgood['DATEVECTOR1'] = pd.to_datetime(tempdtframe, utc=True)
# dfsurfgood[vardict['unixd']] = dfsurfgood['DATEVECTOR1'].astype(
#     'int64') // 10 ** 9
# dfsurfgood[vardict['datevec']] = dfsurfgood['DATEVECTOR1'].dt.strftime(
#     '%Y-%m-%dT%H:%M:%SZ')
#
#
#
#
# # Rename columns
# dfsurfgood.rename(
#     columns={'EXPOCODE': vardict['id'], 'DOI': vardict['doi'],
#              'G2latitude': vardict['lat'], 'G2longitude': vardict['lon'],
#              'G2depth': vardict['dep'], 'G2temperature': vardict['temp'],
#              'G2salinity': vardict['sal'], 'G2salinityf': vardict['salf'],
#              'G2tco2': vardict['dic'], 'G2tco2qc': vardict['dicf'],
#              'G2talk': vardict['alk'], 'G2talkqc': vardict['alkf'],
#              'G2phtsinsitutp': vardict['ph'], 'G2phtsqc': vardict['phf'],
#              'G2fco2': vardict['fco2w'], 'G2fco2f': vardict['fco2wf']},
#     inplace=True)
#
# print(tempdf.shape, dfsurf.shape, dfsurfgood.shape)
#
# # Add source (SOCAT, GLODAP, ARGO, etc...)
# dfsurfgood['SOURCE'] = source
#
# # Rename and reset indices
# printdf = dfsurfgood
# printdf.reset_index(drop=True, inplace=True)
#
# print('GLODAP frame size is ')
# print(printdf.shape)
