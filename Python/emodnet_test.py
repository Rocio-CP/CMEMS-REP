import datetime
import pandas as pd
import os
import re

def read_emodnet_carbon(input_file):
    wanted_vars = ['Depth [m]', 'time_ISO8601',
                   'ITS-90 water temperature [degrees C]', 'Water body salinity [per mille]',
                   'Water body pH [pH units]', 'Water body dissolved inorganic carbon [umol/l]',
                   'Water body total alkalinity [mEquiv/l]']
    wanted_carbon_vars = ['Water body total alkalinity [mEquiv/l]',
                          'Water body dissolved inorganic carbon [umol/l]',
                          'Water body pH [pH units]']
    header_dict = {}

    line = ''
    headerlines = -1
    headertext = "Cruise\tStation\tType"

    f = open(filepath_emodnet)

    while headertext not in line:
        line = f.readline()

        if (line.startswith('//<DataType>')):
            header_dict['DataType'] = re.search('>([a-zA-Z]+?)<', line).group(1)

        # Split the info on each variable into a dictionary key:value
        if (line.startswith('//<MetaVariable>')) | (line.startswith('//<DataVariable>')):
            temp_header_dict = {}
            temp_header_dict['type'] = re.search('<(.+?)>', line).group(1)
            # regex key:val
            reg_search = re.compile(r'([a-z_]+?)="([A-Za-z_0-9 µ:/%^\-\[\]\(\)\.]*?)"')
            reg_search_match = reg_search.findall(line)
            for l in reg_search_match:
                temp_header_dict[l[0]] = l[1]
            header_dict[temp_header_dict['label']] = temp_header_dict

        headerlines = headerlines + 1

    f.close()

    # Lists with the metadata and the wanted variables available in the file
    metadata = []
    for metad in header_dict.keys():
        if metad != 'DataType':
            if header_dict[metad]['type'] == 'MetaVariable':
                metadata.append(metad)
    vars_available = list(set(header_dict.keys()) & set(wanted_vars))
    carbon_vars_available = list(set(header_dict.keys()) & set(wanted_carbon_vars))

    # If no carbon variables are available, skip to next file
    if len(carbon_vars_available) == 0:
        print('No carbon variables')
        return

    # Create dtype dictionary to speed up pd.read_csv. For all variables in the file
    dtype_dict_og = dict(zip([a for a in header_dict.keys() if a != 'DataType'],
                             [header_dict[a]['value_type'] for a in header_dict.keys() if a != 'DataType']))
    dtype_dict = dtype_dict_og.copy()
    # For some reason yyyy-mm-ddThh:mm:ss.sss is not declared in the header info; give the same value_type as time_ISO8601
    dtype_dict['yyyy-mm-ddThh:mm:ss.sss'] = 'DOUBLE'
    # Read the file column header to get the QC column names and add them to the dtype dictionary
    colnames_from_file = pd.read_csv(filepath_emodnet, sep='\t', skiprows=headerlines,
                                     index_col=0, nrows=0)
    qv_vars = [a for a in colnames_from_file if a.__contains__('QV:SEADATANET')]
    dtype_dict_qv = {key: 'TEXT' for key in qv_vars}  # L20 vocab
    dtype_dict.update(dtype_dict_qv)
    # Map to python dtype values
    for key, val in dtype_dict.items():
        if (key.__contains__('yyyy-mm-ddThh:mm:ss')) | (key.__contains__('time')):
            dtype_dict[key] = 'str'
        elif val.__contains__('TEXT'):
            dtype_dict[key] = 'str'
        elif (val.__contains__('FLOAT')) | (val.__contains__('DOUBLE')):
            dtype_dict[key] = 'float'
        elif val.__contains__('INT'):
            dtype_dict[key] = 'Int64'

    # Read the file into a dataframe
    dates_to_parse = [a for a in ['yyyy-mm-ddThh:mm:ss.sss', 'time_ISO8601'] if a in dtype_dict.keys()]
    df = pd.read_csv(filepath_emodnet, sep='\t', skiprows=headerlines,
                     dtype=dtype_dict, parse_dates=dates_to_parse)

    # Identify columns that have pH, DIC, Alk, Temp, Sal and their corresponding QC flag
    # Find indices for variables and their QV columns to rename
    vars_ind = [df.columns.get_loc(varname) for varname in vars_available]
    qcvars_ind = [a + 1 for a in vars_ind]
    # Rename the QCvars with the variable they QC
    qcvars = [var + '_QV:SEADATANET' for var in vars_available]
    df.rename(columns=dict(zip(df.iloc[:, qcvars_ind].columns, qcvars)), inplace=True)

    # Subset only the columns we need (all metadata + wanted variables)
    metadata_ind = [df.columns.get_loc(metad) for metad in metadata]
    df2 = df.loc[:, [*metadata, *vars_available, *qcvars]].copy()

    # Fill the metadata gaps BEFORE selecting the rows with valid carbon data
    df3 = df2.copy()
    df3.loc[:, metadata] = df3.loc[:, metadata].fillna(method='ffill')

    # Extract data points with carbon data (regardless of quality, for now)
    # create the filter expression (remove lines without carbon measurements)
    carbon_filter_expression = " | ".join(
        f"(df3['{carbon_qvvar}_QV:SEADATANET'] != '9')" for carbon_qvvar in carbon_vars_available)
    dum = df3.loc[eval(carbon_filter_expression)].copy()

    # Column with file source
    dum['source file']=filepath_emodnet.split('/')[-1]

    return dum

# Script here
start_time = datetime.datetime.now()

emodnet_filespath = '/Users/rocio/Documents/templocal/emodnetchem/'
emodnet_files = [f for f in os.listdir(emodnet_filespath) if f.startswith('Eutrophication_')]

# Med_profiles and timeseries has several lines with extra tab (109 cols vs 108 / 105 vs 104) Run a quick bash script to find those lines with 109 cols. Remove extra tabs manually
# Atlantic_profiles has 106 columns instead of 107
# for f in Eutrophication*.txt; do
# cat ${f} | awk 'BEGIN {FS="\t"}; {print NF}' > number_fields_${f}.txt
# grep -n "109" number_fields_${f}
# done

df4 = pd.DataFrame()

for file in emodnet_files:
    print(file)
    filepath_emodnet = os.path.join(emodnet_filespath, file)
    df3 = read_emodnet_carbon(filepath_emodnet)
    df4 = pd.concat([df4, df3])

df4.reset_index(drop=True, inplace=True)


lapsed = datetime.datetime.now() - start_time
print(lapsed.total_seconds() / 60)

# Create a datetime column; if yyyy-mm-ddThh:mm:ss.sss and time_ISO8601, pick the latter.



# Reorder columns to have the QV flags next to their variables
# Create NetCDF?? Use header_dict (now within the read_emodnet function) for attributes


