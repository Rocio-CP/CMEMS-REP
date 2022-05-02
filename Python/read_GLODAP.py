import itertools
import os
import numpy as np
import pandas as pd
import netCDF4
import Python.CMEMS_dictionaries as CMEMSdict

import datetime
tic=datetime.datetime.now()

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

[dimension_dict,variables_dict,G2OS_flag_dict]=CMEMSdict.common_dictionaries()


# It may not be necessary to specify the datatype for reading: it'll transform when writing the netCDF file
all_glodap_input_variables=pd.read_csv(os.path.join(input_files_dir, sourcefile),
                                       index_col=0, nrows=0).columns.tolist()
# pd.Int16Dtype() is the type of pandas integer that allows for na
all_glodap_input_variables_dtype=['UInt64']*all_glodap_input_variables.__len__()
for n,gi in enumerate(all_glodap_input_variables):
    if not gi.endswith('f') and 'qc' not in gi and gi not in ['G2cruise','G2region','G2cast','G2year','G2month','G2day','G2hour','G2minute','G2bottle']:
           all_glodap_input_variables_dtype[n]='float'
data_dtype_dict=dict(zip(all_glodap_input_variables, all_glodap_input_variables_dtype))
# Dtype for GLODAP_info
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
#dois = pd.read_csv(os.path.join(input_files_dir, doifile), sep='\t', header=None,
#                   names=['G2cruise', 'doi'], dtype={'G2cruise':'int','doi':'str'}, on_bad_lines='skip',
#                   encoding='utf_16_le')
expocodes = pd.read_csv(os.path.join(input_files_dir, expocodefile), sep='\t',
                        header=None, names=['G2cruise', 'expocode'], dtype={'G2cruise':'int','expocode':'str'},
                        on_bad_lines='skip')

### Assign Expocodes and DOIs to data frame
# Transform expocode dataframes into dictionaries (easier to lookup)
expocodesdict = expocodes.set_index('G2cruise')['expocode'].to_dict()
#doisdict = dois.set_index('G2cruise')['doi'].to_dict()

for cruise in tempdf['G2cruise'].unique():
    tempdf.loc[tempdf.G2cruise == cruise,
                   'expocode'] = expocodesdict[cruise]
#    tempdf.loc[tempdf.G2cruise == cruise,
#                   'doi'] = doisdict[cruise].rsplit('.org/', 1)[1]

# GLODAP info: platform names, codes, EDMO, PI...
GLODAP_info = pd.read_csv(os.path.join('/Users/rpr061/Downloads', info_file), sep='\t',
                          skiprows=0, dtype=info_dtype_dict, na_values=-9999,
                          on_bad_lines='skip')

# Create platform_code series: Callsign if exists, name if it doesn't
GLODAP_info['PlatformCode']=GLODAP_info['CallSign_WMO']
GLODAP_info.loc[GLODAP_info['PlatformCode'].isna(),'PlatformCode']=\
    GLODAP_info.loc[GLODAP_info['PlatformCode'].isna(),'Name']
GLODAP_info['PlatformCode']=GLODAP_info['PlatformCode'].str.replace(' ','')


### Extract from each platform
# loop through the platforms (the basis of INSTAC files)
in_CMEMS_variables=['G2datetime', 'G2latitude', 'G2longitude', 'G2depthnominal', *variables_dict['flag'].keys(),
                    'G2datetimef','G2positionf','G2depthnominalf',*variables_dict['flag'].values()]
#loop_dimensions=['G2datetime', 'G2latitude', 'G2longitude', 'G2depthnominal']
#loop_variables=input_vars

unique_platform_codes=GLODAP_info['PlatformCode'].unique()

for pc in unique_platform_codes:
    # Create boolean filters for
    filter_expocodes=GLODAP_info['PlatformCode']==pc
    current_expocodes_info=GLODAP_info.loc[filter_expocodes,]
    filter_expocodes_data=tempdf['expocode'].isin(current_expocodes_info['EXPOCODE'])

    # Extract sub-dataframe
    current_dataframe = tempdf.loc[filter_expocodes_data, in_CMEMS_variables]

    global_attributes = CMEMSdict.global_attributes_dictionary(current_expocodes_info, current_dataframe)

    if global_attributes['source_platform_category_code'] == '31' :
        nc_filename=output_files_dir+'VESSEL/'+global_attributes['ID']
    else :
        nc_filename=output_files_dir+'ETC/'+global_attributes['ID']

    # Create NetCDF file
   # nc_filename = output_file_full_name
    nc = netCDF4.Dataset(nc_filename, format="NETCDF4_CLASSIC", mode="w")

    # Create dimensions and their variables
    for d in dimension_dict['dim_name'].keys():
        dimension_attributes = CMEMSdict.dimension_attributes_dictionary(variable_name)


        dim=nc.createDimension(dimension_dict['dim_name'][d],
                           current_dataframe[d].sort_values().unique().__len__())
        dim_variable=nc.CreateVariable(dimension_dict['dim_var_name'][d],
                                       dimension_attributes[dimension_dict['dim_var_name'][d]][0]["datatype"],
                                       dimension_attributes[dimension_dict['dim_var_name'][d]][0]["dimensions"])
        dim_variable=current_dataframe[d].sort_values().unique()

        for key, value in dimension_attributes[1].items():
            dim_variable.setncattr(key, value)


    #timedim = nc.createDimension("TIME", current_dataframe['G2datetime'].sort_values().unique().__len__())
    #depthdim = nc.createDimension("DEPTH", current_dataframe['G2depthnominal'].sort_values().unique().__len__())
#    depthvar = nc.createVariable("DEPH", "float", "DEPTH")
#    depthvar[:]=current_dataframe['G2depthnominal'].sort_values().unique()
#    timevar = nc.createVariable("TIME", "float", "TIME")
#    timevar[:]=current_dataframe['G2datetime'].sort_values().unique()

    # Create variables
    for variable_ind, variable_name in enumerate(variables_dict['SDN'].keys()):

        variable_attributes, QC_attributes = CMEMSdict.variable_attributes_dictionary(variable_name)

        variable = nc.createVariable(variables_dict['SDN'][variable_name],
                                     variable_attributes[0]["datatype"],
                                     variable_attributes[0]["dimensions"])
        value_table = pd.pivot_table(current_dataframe, index="G2datetime", columns="G2depthnominal",
                                     values=variable_name, aggfunc="mean")
        value_table_nobs = pd.pivot_table(current_dataframe, index="G2datetime", columns="G2depthnominal",
                                     values=variable_name, aggfunc="count")
        value_table[value_table.isna()] = netCDF4.default_fillvals['i4']
        variable[:] = value_table*1000

        #variable_attributes = attributes[variable_name + "_varatt"][1]
        for key, value in variable_attributes[1].items():
            variable.setncattr(key, value)

        # QC variables
        value_qc_table= pd.pivot_table(current_dataframe, index="G2datetime", columns="G2depthnominal",
                                     values=variables_dict['flag'][variable_name], aggfunc="mean")
        value_qc_table[value_table_nobs > 1]= 8
        value_qc_table[value_qc_table.isna()] = 9
#        value_qc_table[value_table.isna()] = netCDF4.default_fillvals['i1']
        variableqc = nc.createVariable(variable_name+"_QC",
                                     QC_attributes[0]["datatype"],
                                     QC_attributes[0]["dimensions"])
        variableqc[:] = value_qc_table

        # Variable (and QC variable) attributes
        for key, value in QC_attributes[1].items():
            variableqc.setncattr(key, value)

    # Global attributes
    for key, value in global_attributes.items():
        nc.setncattr(key, value)

    nc.close()
    break

toc = datetime.datetime.now() - tic
print (toc)
