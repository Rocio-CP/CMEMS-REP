# SOCAT import. From local .tsv file or ERDDAP

import os
import pandas as pd
import time
import numpy as np
from erddapy import ERDDAP

input_dir='/Users/rocio/Downloads/SOCATv2022_synthesis_files'
source='SOCATv2022_Indian'
datafrom ='local'
headerlines = 0
# Link to server "https://www.ncei.noaa.gov/data/oceans/ncei/ocads/data/0235360/"

# Find where the files are and retrieve them (synthesis and FlagE)
input_dir='/Users/rocio/Downloads/SOCATv2022_synthesis_files'
file='SOCATv2022_Indian.tsv'
#SOCAT_files = os.listdir(os.path.join(input_dir))
SOCAT_files=[os.path.join(input_dir,file)]

##### -----------------------------
# Don't touch from here
SOCATcolheadertextshort = 'Expocode\tversion\tSource_DOI\tQC_Flag'
SOCATmetacolheadertextshort = 'Expocode\tversion\tDataset'


for file in SOCAT_files:
    filepath = os.path.join(input_dir, file)

    # Read metadata header for Cruise Flags and find the number of
    # headerlines before the data
    separator = '\t'
    line = ''
    headerlines = -1

    f = open(filepath)
    while SOCATcolheadertextshort not in line:
        line = f.readline()
        headerlines = headerlines + 1
    f.close()

    # Read SOCAT data in dataframe
    print(file + " file has " + str(headerlines) + " header lines")
    start_time = time.time()  # Time the script
    ddtype = {0: str, 2: str}  # add type str to columns 0 and 2
    # Read the SOCAT file into a pandas dataframe
    tempdf1 = pd.read_csv(filepath, sep = separator, skiprows = headerlines,
        na_values='NaN',on_bad_lines='skip', dtype=ddtype)
    print("--- %s seconds ---" % (time.time() - start_time))
    print(file + " data frame has " + str(len(tempdf1)) + " lines")

    ### Create datetime from 1950 series
    # Create days-from-1950 numeric datetime
    tempdtframe = pd.to_datetime({'year': tempdf1['yr'], 'month': tempdf1['mon'],
        'day': tempdf1['day'], 'hour': tempdf1['hh'],
        'minute': tempdf1['mm'], 'seconds': tempdf1['ss']}, utc=True)
    timedelta_1950 = tempdtframe - pd.Timestamp('1950-01-01T00:00:00', tz='UTC')
    tempdf1['datetime'] = timedelta_1950.dt.total_seconds() / 86400

    # Transform longitude to +-180
    tempdf1.loc[tempdf1['longitude [dec.deg.E]'] > 180, 'longitude [dec.deg.E]'] = \
        tempdf1['longitude [dec.deg.E]'] - 360

    # Merge synthesis and FlagE dataframes in one
    if 'tempdf' not in globals():
        tempdf = tempdf1
    else:
        tempdf = tempdf.append(tempdf1)


# Rename and reset indices
printdf = tempdf
printdf.reset_index(drop=True, inplace=True)

print('SOCAT frame size is ')
print(printdf.shape)

if global_attributes['source_platform_category_code'] == '31':
    nc_filename = output_files_dir + 'VESSEL/' + global_attributes['id']
else:
    nc_filename = output_files_dir + 'ETC/' + global_attributes['id']