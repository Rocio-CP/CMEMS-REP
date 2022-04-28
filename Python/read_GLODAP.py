import itertools
import json
# Main script calling to
import os

import numpy as np
import pandas as pd
from netCDF4 import Dataset
import netCDF4

# Some parameters / unchanging values / dictionaries
dimensions=('G2datetime', 'G2depthnominal','G2latitude','G2longitude')
dimension_flags=('G2datetimef','G2depthnominalf','G2positionf','G2positionf')
SDN_dimension_names=('TIME','DEPH','LATITUDE','LONGITUDE')
SDN_dimension_variable_names=('TIME','DEPTH','LATITUDE','LONGITUDE')
dimension_units=("days since 1950-01-01T00:00:00Z","m","degree_north","degree_east")
dimension_long_name=("Time","Depth","Latitude of each location","Longitude of each location")
dimension_CF_standard_name=("time","depth","latitude","longitude")
dimension_axis=("T","Z","Y","X")
dimension_valid_max=(90000.0,12000.0,90.0,180.0)
dimension=valid_min=(-90000.0,-12000.0,-90.0,-180.0)
dimension_dict = {}
dimension_dict['flag'] = dict(zip(dimensions, dimension_flags))
dimension_dict['SDN'] = dict(zip(dimensions, SDN_dimension_names))
dimension_dict['unit'] = dict(zip(dimensions, dimension_units))
dimension_dict['long'] = dict(zip(dimensions, dimension_long_name))
dimension_dict['CF'] = dict(zip(dimensions, dimension_CF_standard_name))



input_vars = ('G2bottomdepth', 'G2pressure', 'G2temperature', 'G2salinity', 'G2oxygen', 'G2nitrate',
              'G2nitrite', 'G2phosphate', 'G2silicate', 'G2phtsinsitutp', 'G2phts25p0', 'G2tco2',
              'G2talk', 'G2doc', 'G2don', 'G2tdn', 'G2chla')
input_flag_vars = ('G2bottomdepthf', 'G2pressuref', 'G2temperaturef', 'G2salinityf', 'G2oxygenf', 'G2nitratef',
                   'G2nitritef', 'G2phosphatef', 'G2silicatef', 'G2phtsinsitutpf', 'G2phts25p0f', 'G2tco2f',
                   'G2talkf', 'G2docf', 'G2donf', 'G2tdnf', 'G2chlaf')
SDN_var_names = ('BATH', 'PRES','TEMP', 'PSAL', 'DOX2', 'NTAW',
                 'NTIW', 'PHOW', 'SLCW', 'PHPH', 'PH25', 'TICW',
                 'ALKW', 'CORG', 'NODW', 'NT1D', 'CPHL')
units = ('m', 'dbar','degrees_C', '0.001', 'µmol kg-1', 'µmol kg-1',
         'µmol kg-1', 'µmol kg-1', 'µmol kg-1', '1', '1', 'µmol kg-1',
         'µmol kg-1', 'µmol kg-1', 'µmol kg-1', 'µmol kg-1', 'mg m-3')
long_name = ('Bathymetric depth', 'Sea pressure',
             'Sea temperature', 'Practical salinity',
             'Dissolved oxygen', 'Nitrate (NO3-N)',
             'Nitrite (NO2-N)', 'Phosphate (PO4-P)', 'Silicate (SIO4-SI)',
             'Ph', 'Ph at 25 °C and 0 dbar', 'Dissolved inorganic carbon',
             'Total alkalinity', 'Dissolved organic carbon', 'Dissolved organic nitrogen',
             'Total dissolved nitrogen', 'Chlorophyll-a')
CF_standard_name = ('sea_floor_depth_below_sea_surface', 'sea_water_pressure','sea_water_temperature',
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
variables_dict['SDN'] = dict(zip(input_vars, SDN_var_names))
variables_dict['unit'] = dict(zip(input_vars, units))
variables_dict['long'] = dict(zip(input_vars, long_name))
variables_dict['CF'] = dict(zip(input_vars, CF_standard_name))

GLODAP_flags = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
OceanSITES_flags = [8, 9, 1, 2, 4, 9, 8, 0, 0, 9]
G2OS_flag_dict = dict(zip(GLODAP_flags, OceanSITES_flags))



input_files_dir = '/Users/rpr061/Documents/'
files_path_remote = ('https://www.ncei.noaa.gov'
                     '/data/oceans/ncei/ocads/data/0237935/')
product_dir='/Users/rpr061/Documents/INSITU_GLO_BGC_CARBON_DISCRETE_MY_013_050/'
# Create and store output directory.
if not os.path.isdir(product_dir):
    os.makedirs(os.path.join(product_dir,'DNT/'))
    os.makedirs(os.path.join(product_dir,'cmems_obs-ins_glo_bgc-car_my_glodap-gridded_irr/'))
    os.makedirs(os.path.join(product_dir,'cmems_obs-ins_glo_bgc-car_my_glodap-obs_irr/VESSEL'))
    os.makedirs(os.path.join(product_dir,'cmems_obs-ins_glo_bgc-car_my_glodap-obs_irr/ETC'))
    os.makedirs(os.path.join(product_dir,'cmems_obs-ins_glo_bgc-car_my_socat-gridded_irr/'))
    os.makedirs(os.path.join(product_dir,'cmems_obs-ins_glo_bgc-car_my_socat-obs_irr/VESSEL/'))
    os.makedirs(os.path.join(product_dir,'cmems_obs-ins_glo_bgc-car_my_socat-obs_irr/MOORING/'))
    os.makedirs(os.path.join(product_dir,'cmems_obs-ins_glo_bgc-car_my_socat-obs_irr/DRIFTER/'))
    os.makedirs(os.path.join(product_dir,'cmems_obs-ins_glo_bgc-car_my_socat-obs_irr/AUV/'))

output_files_dir = os.path.join(product_dir,'cmems_obs-ins_glo_bgc-car_my_glodap-obs_irr/')

GLODAP_files = ['GLODAPv2.2021_Merged_Master_File.csv', 'EXPOCODES.txt', 'Dataset_DOIs.txt']
info_file = 'GLODAPv22021CMEMS.tsv'


# Check if input data files are there; if not, fetch data from the internet
if not all([os.path.isfile(os.path.join(input_files_dir, f)) for f in GLODAP_files]):
    input_files_dir = files_path_remote

sourcefile = [i for i in GLODAP_files if 'GLODAP' in i][0]
doifile = [i for i in GLODAP_files if 'DOI' in i][0]
expocodefile = [i for i in GLODAP_files if 'EXPOCODES' in i][0]

all_glodap_input_variables=pd.read_csv(os.path.join(input_files_dir, sourcefile), index_col=0, nrows=0).columns.tolist()
# pd.Int16Dtype() is the type of pandas integer that allows for na
all_glodap_input_variables_dtype=['UInt64']*all_glodap_input_variables.__len__()
for n,gi in enumerate(all_glodap_input_variables):
    if not gi.endswith('f') and 'qc' not in gi and gi not in ['G2cruise','G2region','G2cast','G2year','G2month','G2day','G2hour','G2minute','G2bottle']:
           all_glodap_input_variables_dtype[n]='float'
data_dtype_dict=dict(zip(all_glodap_input_variables, all_glodap_input_variables_dtype))

glodap_info_dtype=['int',*['str']*12]
info_dtype_dict=dict(zip(range(0, 13), glodap_info_dtype))

# Read files
# parse_dates does not work because some hours and minutes are nan
tempdf = pd.read_csv(os.path.join(input_files_dir, sourcefile), sep=',',
                     skiprows=0, na_values=-9999, # dtype=data_dtype_dict,
                     on_bad_lines='skip')

### Create datetime from 1950 series
# Change NaN hour/minute to 0. Calculate 1950-time
tempdf.loc[np.isnan(tempdf['G2hour']), 'G2hour'] = 0
tempdf.loc[np.isnan(tempdf['G2minute']), 'G2minute'] = 0
# Create days-from-1950 numeric datetime
tempdtframe = pd.to_datetime({'year': tempdf['G2year'], 'month': tempdf['G2month'],
     'day': tempdf['G2day'], 'hour': tempdf['G2hour'],
     'minute': tempdf['G2minute']}, utc=True)
timedelta_1950 = tempdtframe - pd.Timestamp('1950-01-01T00:00:00', tz='UTC')
tempdf['G2datetime'] = timedelta_1950.dt.total_seconds() / 86400

# Calculate the nominal depths series
nominal_depths = tuple(itertools.chain(range(0, 105, 5), range(100, 1025, 25), range(1000, 10100, 100)))
nominal_depths_array = np.array(nominal_depths)
actual_depths_array = np.array(tempdf['G2depth'])
depths_as_nominal = nominal_depths_array[
    abs(actual_depths_array[None, :] - nominal_depths_array[:, None]).argmin(axis=0)]
tempdf['G2depthnominal'] = depths_as_nominal

### QC flags
# If temperature and bottom depth exist, they're assumed good (flag 1)
# If value not recorded/measured, then flag 9
tempdf['G2bottomdepthf'] = np.ones(tempdf['G2bottomdepth'].__len__())
tempdf.loc[np.isnan(tempdf['G2bottomdepth']), 'G2bottomdepthf'] = 9
tempdf['G2temperaturef'] = np.ones(tempdf['G2temperature'].__len__())
tempdf.loc[np.isnan(tempdf['G2temperature']), 'G2temperaturef'] = 9
# Needed for CMEMS: time, depth and position flags.
tempdf['G2datetimef']=np.ones(tempdf['G2datetime'].__len__())
tempdf['G2positionf']=np.ones(tempdf['G2latitude'].__len__())
tempdf['G2depthnominalf']=np.ones(tempdf['G2depthnominal'].__len__())*7
tempdf.loc[np.isnan(tempdf['G2depthnominal']), 'G2depthnominalf'] = 9
tempdf['G2pressuref']=np.ones(tempdf['G2pressure'].__len__())
tempdf.loc[np.isnan(tempdf['G2pressure']), 'G2pressuref'] = 9

# Remap flags to OceanSITES NaN -> 0; 0->8; 1->9; 2->1; 3->2; 4->4;
# % 5->9; 6->8; 7->0; 8->0; 9->9 as int8 (equivalent to byte)
for fv in variables_dict['flag']:
    if fv:
        tempdf.loc[np.isnan(tempdf[fv]), [fv]] = 9
        for fvv in G2OS_flag_dict.keys():
            tempdf.loc[tempdf[fv] == fvv, [fv]] = G2OS_flag_dict[fvv]

# Rename G2fco2 to G2fco2_20_0 (in GLODAP, it's given at 20 dg, 0dbar).
tempdf.rename(columns={'G2fco2': 'G2fco2_20_0'}, inplace=True)

# Read extra info files
dois = pd.read_csv(os.path.join(input_files_dir, doifile), sep='\t', header=None,
                   names=['G2cruise', 'doi'], dtype={'G2cruise':'int','doi':'str'}, on_bad_lines='skip',
                   encoding='utf_16_le')
expocodes = pd.read_csv(os.path.join(input_files_dir, expocodefile), sep='\t',
                        header=None, names=['G2cruise', 'expocode'], dtype={'G2cruise':'int','expocode':'str'},
                        on_bad_lines='skip')

### Assign Expocodes and DOIs to data frame
# Transform expocode dataframes into dictionaries (easier to lookup)
expocodesdict = expocodes.set_index('G2cruise')['expocode'].to_dict()
doisdict = dois.set_index('G2cruise')['doi'].to_dict()

for cruise in tempdf['G2cruise'].unique():
    tempdf.loc[tempdf.G2cruise == cruise,
                   'expocode'] = expocodesdict[cruise]
    tempdf.loc[tempdf.G2cruise == cruise,
                   'doi'] = doisdict[cruise].rsplit('.org/', 1)[1]

# GLODAP info: platform names, codes, EDMO, PI...
GLODAP_info = pd.read_csv(os.path.join('/Users/rpr061/Downloads', info_file), sep='\t',
                          skiprows=0, dtype=info_dtype_dict, na_values=-9999,
                          on_bad_lines='skip')

# Create platform_code series: Callsign if exists, name if it doesn't
GLODAP_info['PlatformCode']=GLODAP_info['CallSign_WMO']
GLODAP_info.loc[GLODAP_info['PlatformCode'].isna(),'PlatformCode']=\
    GLODAP_info.loc[GLODAP_info['PlatformCode'].isna(),'Name']
GLODAP_info['PlatformCode']=GLODAP_info['PlatformCode'].str.replace(' ','')

### Global and attribute netcdf variables dictionary
with open("CMEMS_INSTAC_metadata.json", "r") as template:
    all_attributes = json.load(template)
template.close()

### Extract from each platform
# loop through the platforms (the basis of INSTAC files)
unique_platform_codes=GLODAP_info['PlatformCode'].unique()
in_CMEMS_variables=['G2datetime', 'G2latitude', 'G2longitude', 'G2depthnominal', *input_vars,
                    'G2datetimef','G2positionf','G2depthnominalf',*input_flag_vars]
loop_dimensions=['G2datetime', 'G2latitude', 'G2longitude', 'G2depthnominal']
loop_variables=input_vars

for pc in unique_platform_codes:
    # Create boolean filters for
    filter_expocodes=GLODAP_info['PlatformCode']==pc
    current_expocodes=GLODAP_info.loc[filter_expocodes,'EXPOCODE']
    filter_expocodes_data=tempdf['expocode'].isin(current_expocodes)

    # Extract sub-dataframe
    current_dataframe = tempdf.loc[filter_expocodes_data, in_CMEMS_variables]

    all_attributes['globalatt'][0]['id']='GL_PR_BO_'+pc+'-GLODAPv22022.nc'
    platform_category = str(GLODAP_info.loc[filter_expocodes, 'PlatformType'].unique().item())
    if platform_category == '31' :
        all_attributes['globalatt'][0]['source'] = 'research vessel'
        output_file_full_name=output_files_dir+'VESSEL/'+all_attributes['globalatt'][0]['id']
    elif platform_category == '21' :
        all_attributes['globalatt'][0]['source'] = 'propelled manned submersible'
        output_file_full_name=output_files_dir+'ETC/'+all_attributes['globalatt'][0]['id']
    elif platform_category == '62' :
        all_attributes['globalatt'][0]['source'] = 'aeroplane'
        output_file_full_name=output_files_dir+'ETC/'+all_attributes['globalatt'][0]['id']
    all_attributes['globalatt'][0]['source_platform_category_code'] = platform_category

    # Create NetCDF file
    nc_filename = output_file_full_name
    nc = Dataset(nc_filename, format="NETCDF4_CLASSIC", mode="w")
    timedim = nc.createDimension("TIME", current_dataframe['G2datetime'].sort_values().unique().__len__())
    depthdim = nc.createDimension("DEPTH", current_dataframe['G2depthnominal'].sort_values().unique().__len__())
    depthvar = nc.createVariable("DEPH", "float", "DEPTH")
    depthvar[:]=current_dataframe['G2depthnominal'].sort_values().unique()
    timevar = nc.createVariable("TIME", "float", "TIME")
    timevar[:]=current_dataframe['G2datetime'].sort_values().unique()

    # Generate variables and attach the attributes
    for variable_ind, variable_name in enumerate(loop_variables):
        variable = nc.createVariable(variables_dict['SDN'][variable_name],
                                     all_attributes['variable_varatt'][0]["datatype"],
                                     all_attributes['variable_varatt'][0]["dimensions"])
        value_table = pd.pivot_table(current_dataframe, index="G2datetime", columns="G2depthnominal",
                                     values=variable_name, aggfunc="mean")
        value_table_nobs = pd.pivot_table(current_dataframe, index="G2datetime", columns="G2depthnominal",
                                     values=variable_name, aggfunc="count")
        value_table[value_table.isna()] = netCDF4.default_fillvals['i4']
        variable[:] = value_table*1000

        # Variable attributes:
        all_attributes['variable_varatt'][1]["standard_name"]=variables_dict['CF'][variable_name]
        all_attributes['variable_varatt'][1]["units"]=variables_dict['unit'][variable_name]
        all_attributes['variable_varatt'][1]["long_name"]=variables_dict['long'][variable_name]
        all_attributes['variable_varatt'][1]["ancillary_variables"]=variable_name+"_QC"

        #variable_attributes = attributes[variable_name + "_varatt"][1]
        for key, value in all_attributes['variable_varatt'][1].items():
            variable.setncattr(key, value)

        value_qc_table= pd.pivot_table(current_dataframe, index="G2datetime", columns="G2depthnominal",
                                     values=variables_dict['flag'][variable_name], aggfunc="mean")
        value_qc_table[value_table_nobs > 1]= 8
        value_qc_table[value_qc_table.isna()] = 9
        variableqc = nc.createVariable(variable_name+"_QC",
                                     all_attributes['QC_varatt'][0]["datatype"],
                                     all_attributes['QC_varatt'][0]["dimensions"])
        variableqc[:] = value_qc_table


        for key, value in all_attributes['variable_varatt'][1].items():
            variable.setncattr(key, value)

    # Platform-specific global attributes
    # Platform attributes
    all_attributes['globalatt'][0]['platform_code']= pc # only alphanumeric characters
    all_attributes['globalatt'][0]['platform_name']= GLODAP_info.loc[filter_expocodes, 'Name'].unique().item()
    all_attributes['globalatt'][0]['wmo_platform_code']=GLODAP_info.loc[filter_expocodes, 'CallSign_WMO'].unique().item()
    all_attributes['globalatt'][0]['ices_platform_code']=GLODAP_info.loc[filter_expocodes, 'ICEScode'].unique().item()
    #if not wmo : wmo= ''
    #if not ices : ices= ''

    # Institution attributes
    institution_name=GLODAP_info.loc[filter_expocodes, 'Institution'].unique()
    edmo=GLODAP_info.loc[filter_expocodes, 'EDMO'].unique()
    pi_institution_name=GLODAP_info.loc[filter_expocodes, 'PI_Institution'].unique()
    pi_edmo=GLODAP_info.loc[filter_expocodes, 'PI_EDMO'].unique()
    all_institution_name=";".join(np.unique([*institution_name, *pi_institution_name]))
    all_edmo=",".join(np.unique([*edmo, *pi_edmo]))
    all_attributes['globalatt'][0]['institution'] =";".join(np.unique([*institution_name, *pi_institution_name]))
    all_attributes['globalatt'][0]['institution_edmo_code'] = ",".join(np.unique([*edmo, *pi_edmo]))
    # PI attributes
    pi=GLODAP_info.loc[filter_expocodes, 'PI'].unique()
    all_attributes['globalatt'][0]['pi_name'] =";".join(pi)



    for key, value in all_attributes['globalatt'][0].items():
        nc.setncattr(key, value)

    nc.close()
    break


