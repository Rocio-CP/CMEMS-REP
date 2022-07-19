import itertools
import os
import numpy as np
import pandas as pd
import netCDF4
import Python.CMEMS_dictionaries as CMEMSdict

# def main():
# Stablish where the files are / should be
input_files_dir = '/Users/rocio/Documents/'
files_path_remote = ('https://www.ncei.noaa.gov'
                     '/data/oceans/ncei/ocads/data/0237935/')
product_dir = '/Users/rocio/Documents/INSITU_GLO_BGC_CARBON_DISCRETE_MY_013_050/'
# Create and store output directory.
if not os.path.isdir(product_dir):
    os.makedirs(os.path.join(product_dir, 'DNT/'))
    os.makedirs(os.path.join(product_dir, 'cmems_obs-ins_glo_bgc-car_my_glodap-gridded_irr/'))
    os.makedirs(os.path.join(product_dir, 'cmems_obs-ins_glo_bgc-car_my_glodap-obs_irr/BO'))
    os.makedirs(os.path.join(product_dir, 'cmems_obs-ins_glo_bgc-car_my_socat-gridded_irr/'))
    os.makedirs(os.path.join(product_dir, 'cmems_obs-ins_glo_bgc-car_my_socat-obs_irr/CO/'))
    os.makedirs(os.path.join(product_dir, 'cmems_obs-ins_glo_bgc-car_my_socat-obs_irr/DB/'))
    os.makedirs(os.path.join(product_dir, 'cmems_obs-ins_glo_bgc-car_my_socat-obs_irr/GL/'))
    os.makedirs(os.path.join(product_dir, 'cmems_obs-ins_glo_bgc-car_my_socat-obs_irr/MO/'))
    os.makedirs(os.path.join(product_dir, 'cmems_obs-ins_glo_bgc-car_my_socat-obs_irr/SD/'))

output_files_dir = os.path.join(product_dir, 'cmems_obs-ins_glo_bgc-car_my_glodap-obs_irr/')

SOCAT_files = ['SOCATv2022.tsv','SOCATv2022_FlagE.tsv']
SOCAT_info_file='SOCATv2021CMEMS.tsv'
GLODAP_files = ['GLODAPv2.2021_Merged_Master_File.csv', 'EXPOCODES.txt', 'Dataset_DOIs.txt']
GLODAP_info_file = 'GLODAPv22021CMEMS.tsv'

# Check if input data files are there; if not, fetch data from the internet
# if not all([os.path.isfile(os.path.join(input_files_dir, f)) for f in GLODAP_files]):
#    input_files_dir = files_path_remote

sourcefile = [i for i in GLODAP_files if 'GLODAP' in i][0]
doifile = [i for i in GLODAP_files if 'DOI' in i][0]
expocodefile = [i for i in GLODAP_files if 'EXPOCODES' in i][0]

variables_dict = CMEMSdict.common_dictionaries('GLODAP')

### Read GLODAP data file(s)
# It may not be necessary to specify the datatype for reading: it'll transform when writing the netCDF file
all_glodap_input_variables = pd.read_csv(os.path.join(input_files_dir, sourcefile),
                                         index_col=0, nrows=0).columns.tolist()

# pd.Int16Dtype() is the type of pandas integer that allows for na
all_glodap_input_variables_dtype = ['UInt64'] * all_glodap_input_variables.__len__()
for n, gi in enumerate(all_glodap_input_variables):
    if not gi.endswith('f') and 'qc' not in gi and gi not in ['G2cruise', 'G2region', 'G2cast', 'G2year', 'G2month',
                                                              'G2day', 'G2hour', 'G2minute', 'G2bottle']:
        all_glodap_input_variables_dtype[n] = 'float'

# dtype for GLODAP data
data_dtype_dict = dict(zip(all_glodap_input_variables, all_glodap_input_variables_dtype))
# Dtype for GLODAP info
glodap_info_dtype = ['int', *['str'] * 12]
info_dtype_dict = dict(zip(range(0, 13), glodap_info_dtype))

# Read GLODAP file (parse_dates does not work because some hours and minutes are nan)
tempdf = pd.read_csv(os.path.join(input_files_dir, sourcefile), sep=',',
                     skiprows=0, na_values=-9999,  # dtype=data_dtype_dict,
                     on_bad_lines='skip')

# Create datetime from 1950 series
# Change NaN hour/minute to 0. Calculate 1950-time
tempdf.loc[np.isnan(tempdf['G2hour']), 'G2hour'] = 0
tempdf.loc[np.isnan(tempdf['G2minute']), 'G2minute'] = 0
# Create days-from-1950 numeric datetime
tempdtframe = pd.to_datetime({'year': tempdf['G2year'], 'month': tempdf['G2month'],
                              'day': tempdf['G2day'], 'hour': tempdf['G2hour'],
                              'minute': tempdf['G2minute']}, utc=True)
timedelta_1950 = tempdtframe - pd.Timestamp('1950-01-01T00:00:00', tz='UTC')
tempdf['TIME'] = timedelta_1950.dt.total_seconds() / 86400

# Calculate the nominal depths series
nominal_depths = tuple(itertools.chain(range(0, 105, 5), range(100, 1025, 25), range(1000, 10100, 100)))
nominal_depths_array = np.array(nominal_depths)
actual_depths_array = np.array(tempdf['G2depth'])
depths_as_nominal = nominal_depths_array[
    abs(actual_depths_array[None, :] - nominal_depths_array[:, None]).argmin(axis=0)]
tempdf['DEPH'] = depths_as_nominal

# Rename latitude/longitude series
tempdf['LATITUDE'] = tempdf['G2latitude']
tempdf['LONGITUDE'] = tempdf['G2longitude']

### QC flags
# Remap existing flags to OceanSITES NaN -> 0; 0->8; 1->9; 2->1; 3->2; 4->4;
# % 5->9; 6->8; 7->0; 8->0; 9->9 as int8 (equivalent to byte)
for fvar in variables_dict['flag'].values():
    if fvar in tempdf.columns:
        # Turn all flag variables into integers
        #tempdf[fvar].astype('int')
        tempdf.loc[np.isnan(tempdf[fvar]), [fvar]] = 9
        for fvval in variables_dict['flag_values'].keys():
            tempdf.loc[tempdf[fvar] == fvval, [fvar]] = variables_dict['flag_values'][fvval]

# If temperature and bottom depth exist, they're assumed good (flag 1)
# If value not recorded/measured, then flag 9
tempdf['G2bottomdepthf'] = np.ones(tempdf['G2bottomdepth'].__len__())
tempdf.loc[np.isnan(tempdf['G2bottomdepth']), 'G2bottomdepthf'] = 9
tempdf['G2bottomdepthf']=tempdf['G2bottomdepthf'].astype('int')
tempdf['G2temperaturef'] = np.ones(tempdf['G2temperature'].__len__())
tempdf.loc[np.isnan(tempdf['G2temperature']), 'G2temperaturef'] = 9
tempdf['G2temperaturef']=tempdf['G2temperaturef'].astype('int')
# Needed for CMEMS: time, depth and position flags.
tempdf['G2pressuref'] = np.ones(tempdf['G2pressure'].__len__())
tempdf.loc[np.isnan(tempdf['G2pressure']), 'G2pressuref'] = 9
tempdf['G2pressuref']=tempdf['G2pressuref'].astype('int')

# Rename G2fco2 to G2fco2_20_0 (in GLODAP, it's given at 20 dg, 0dbar).
tempdf.rename(columns={'G2fco2': 'G2fco2_20_0'}, inplace=True)

# Read extra info files
# dois = pd.read_csv(os.path.join(input_files_dir, doifile), sep='\t', header=None,
#                   names=['G2cruise', 'doi'], dtype={'G2cruise':'int','doi':'str'}, on_bad_lines='skip',
#                   encoding='utf_16_le')
expocodes = pd.read_csv(os.path.join(input_files_dir, expocodefile), sep='\t',
                        header=None, names=['G2cruise', 'expocode'], dtype={'G2cruise': 'int', 'expocode': 'str'},
                        on_bad_lines='skip')

### Assign Expocodes and DOIs to data frame
# Transform expocode dataframes into dictionaries (easier to lookup)
expocodesdict = expocodes.set_index('G2cruise')['expocode'].to_dict()
# doisdict = dois.set_index('G2cruise')['doi'].to_dict()

for cruise in tempdf['G2cruise'].unique():
    tempdf.loc[tempdf.G2cruise == cruise,
               'EXPOCODE'] = expocodesdict[cruise]

### Read GLODAP info: platform names, codes, EDMO, PI...
GLODAP_info = pd.read_csv(os.path.join(input_files_dir, GLODAP_info_file), sep='\t',
                          skiprows=0, dtype=info_dtype_dict, na_values=-9999,
                          on_bad_lines='skip')

# Create platform_code series: Callsign if exists, name if it doesn't
GLODAP_info['PlatformCode'] = GLODAP_info['CallSign_WMO']
GLODAP_info.loc[GLODAP_info['PlatformCode'].isna(), 'PlatformCode'] = \
    GLODAP_info.loc[GLODAP_info['PlatformCode'].isna(), 'Name']
GLODAP_info['PlatformCode'] = GLODAP_info['PlatformCode'].str.replace(' ', '')

# Extract only variables that will go to CMEMS
all_in_CMEMS_variables = ['TIME', 'LATITUDE', 'LONGITUDE', 'DEPH', *variables_dict['flag'].keys(),
                          'G2datetimef', 'G2positionf', 'G2depthnominalf', *variables_dict['flag'].values()]
in_CMEMS_variables = list(set(all_in_CMEMS_variables) & set(tempdf.columns))
in_CMEMS_variables.append('EXPOCODE')
tempdf=tempdf.loc[ : , in_CMEMS_variables]

# Call the netCDF builder
import CMEMS_build_nc as ncbuild
#ncbuild.build_nc(tempdf,GLODAP_info,output_files_dir)
infoframe=GLODAP_info
dataframe=tempdf
output_files_dir=output_files_dir

unique_platform_codes = infoframe['PlatformCode'].unique()

for pc in unique_platform_codes:
    # Create boolean filters for
    filter_expocodes = infoframe['PlatformCode'] == pc
    current_expocodes_info = infoframe.loc[filter_expocodes,]
    filter_expocodes_data = dataframe['EXPOCODE'].isin(current_expocodes_info['EXPOCODE'])

    # Extract sub-dataframe and order chronologically
    current_dataframe = dataframe.loc[filter_expocodes_data, ]
    current_dataframe = current_dataframe.sort_values(by="TIME")

    global_attributes = CMEMSdict.global_attributes_dictionary(current_expocodes_info, current_dataframe)

    # Figure out dataset and platform type folder
    # Check if there is only one platform type and name per platform
    platform_type = current_expocodes_info['PlatformType'].unique()
    platform_name  = current_expocodes_info['Name'].unique()
    if platform_type.__len__() > 1:
        raise Exception("There is more than one PlatformType code for the platform " + pc)
    if platform_name.__len__() > 1:
        raise Exception("There is more than one Name for the platform " + pc)

    if 'glodap-obs' in output_files_dir:
        # All GLODAP data come from bottle samples
        platform_type_folder= 'BO'
        nc_filename='GL_PR_' + platform_type_folder + '_' + pc + '-GLODAPv22022.nc'

    elif 'socat-obs' in output_files_dir:
        if  platform_type in ['31','32', '37']:
            platform_type_folder= 'CO'
        elif platform_type in ['3B']:
            if 'SD' in platform_name:
                platform_type_folder= 'GL'
            elif 'WG' in platform_name:
                platform_type_folder = 'SD'
        elif platform_type in ['41']:
            platform_type_folder= 'MO'
        elif platform_type in ['42']:
            platform_type_folder= 'DB'
        nc_filename='GL_PR_' + platform_type_folder + '_' + pc + '-SOCATv2022.nc'

    nc_filepathname = output_files_dir + '/' + platform_type_folder +'/'+ nc_filename


# Create NetCDF file
    if os.path.isfile(nc_filepathname):
        os.remove(nc_filepathname)

    nc = netCDF4.Dataset(nc_filepathname, format="NETCDF4_CLASSIC", mode="w")

    # Create the dimension variables values
    timedim_df = current_dataframe[['TIME', 'LATITUDE', 'LONGITUDE']]
    timedim_df = timedim_df.drop_duplicates(subset=['TIME'])
    timedim_df = timedim_df.sort_values('TIME')

    # Create dimensions and their variables; the input to the function are the dimension variables AS IS (e.g. for GLODAP, DEPH MUST be already a matrix, monotonically increasing, etc)
    depth_matrix = current_dataframe['DEPH'].sort_values().unique()
    depth_matrix = depth_matrix.reshape((1, -1))
    depth_matrix = np.repeat(depth_matrix, timedim_df['TIME'].__len__(), axis=0)

    # create_dimensions(time, lat, lon, depth, nc)
    CMEMSdict.create_dimensions(timedim_df['TIME'], timedim_df['LATITUDE'],
                                timedim_df['LONGITUDE'], depth_matrix, nc)

    # Create variables
    # Import variables_dict
    variables_dict = CMEMSdict.common_dictionaries('GLODAP')

    for variable_ind, variable_name in enumerate(variables_dict['SDN'].keys()):

        # Create variable
        value_table = pd.pivot_table(current_dataframe, index="TIME", columns="DEPH",
                                     values=variable_name, aggfunc="mean", dropna=False)
        value_table_nobs = pd.pivot_table(current_dataframe, index="TIME", columns="DEPH",
                                          values=variable_name, aggfunc="count", dropna=False)

        variable_attributes, QC_attributes = CMEMSdict.variable_attributes_dictionary(variable_name, 'GLODAP')

        variable = nc.createVariable(variables_dict['SDN'][variable_name],
                                     variable_attributes[0]["datatype"],
                                     variable_attributes[0]["dimensions"],
                                     fill_value=netCDF4.default_fillvals['i4'])
                                     #fill_value=variable_attributes[0]["_FillValue"])


        #value_table=value_table1.copy()
        value_table[np.isnan(value_table)] = netCDF4.default_fillvals['i4']*0.001
        #value_table[value_table1.isna()] = variable_attributes[0]["_FillValue"] * 0.001#netCDF4.default_fillvals['i4']*1000
        #variable.setncattr('scale_factor',0.001)
        for key, value in variable_attributes[1].items():
           if value or key=='add_offset':
               variable.setncattr(key, value)

        variable[:] = value_table

        # QC variables
        value_qc_table = pd.pivot_table(current_dataframe, index="TIME", columns="DEPH",
                                        values=variables_dict['flag'][variable_name], aggfunc="mean", dropna=False)
        value_qc_table[value_table_nobs > 1] = 8
        value_qc_table[value_qc_table.isna()] = 9
        #        value_qc_table[value_table.isna()] = netCDF4.default_fillvals['i1']
        variableqc = nc.createVariable(variables_dict['SDN'][variable_name] + "_QC",
                                       QC_attributes[0]["datatype"],
                                       QC_attributes[0]["dimensions"],
                                       fill_value=QC_attributes[0]["_FillValue"])
        variableqc[:] = value_qc_table.to_numpy()

        # Variable (and QC variable) attributes
        for key, value in QC_attributes[1].items():
            variableqc.setncattr(key, value)


    # Global attributes
    for key, value in global_attributes.items():
        nc.setncattr(key, value)

    nc.close()
    break
