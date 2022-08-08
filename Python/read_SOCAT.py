def read_SOCAT_obs(input_files_dir, output_files_dir, SOCAT_files, SOCAT_info_file):
    import os
    import pandas as pd
    import time
    import numpy as np
    import Python.CMEMS_dictionaries as CMEMSdict

    # Read SOCAT files, skipping the header information lines
    SOCATcolheadertextshort = 'Expocode\tversion\tSource_DOI\tQC_Flag'
    SOCATmetacolheadertextshort = 'Expocode\tversion\tDataset'

    tempdf = pd.DataFrame()

    for file in SOCAT_files:
        filepath = os.path.join(input_files_dir, file)

        separator = '\t'
        line = ''
        headerlines = -1

        f = open(filepath)
        while SOCATcolheadertextshort not in line:
            line = f.readline()
            headerlines = headerlines + 1
        f.close()

        # Read SOCAT data in a dataframe
        #print(file + " file has " + str(headerlines) + " header lines")
        start_time = time.time()  # Time the script
        ddtype = {0: str, 2: str}  # add type str to columns 0 and 2
        # Read the SOCAT file into a pandas dataframe
        tempdf1 = pd.read_csv(filepath, sep=separator, skiprows=headerlines,
                              na_values='NaN', on_bad_lines='skip', dtype=ddtype)
        #print("--- %s seconds ---" % (time.time() - start_time))
        #print(file + " data frame has " + str(len(tempdf1)) + " lines")

        # Create days-from-1950 numeric datetime
        tempdtframe = pd.to_datetime({'year': tempdf1['yr'], 'month': tempdf1['mon'],
                                      'day': tempdf1['day'], 'hour': tempdf1['hh'],
                                      'minute': tempdf1['mm'], 'seconds': tempdf1['ss']}, utc=True)
        timedelta_1950 = tempdtframe - pd.Timestamp('1950-01-01T00:00:00', tz='UTC')
        tempdf1['datetime'] = timedelta_1950.dt.total_seconds() / 86400

        # Transform longitude to +-180
        tempdf1.loc[tempdf1['longitude [dec.deg.E]'] > 180, 'longitude [dec.deg.E]'] = \
            tempdf1['longitude [dec.deg.E]'] - 360

        # Merge synthesis and FlagE dataframes in one. Use an empty dataframe as starter
        tempdf = pd.concat([tempdf, tempdf1])

    # Reset indices(?)
    # tempdf.reset_index(drop=True, inplace=True)
    #print('SOCAT frame size is ')
    #print(tempdf.shape)

    # Create days-from-1950 numeric datetime
    tempdtframe = pd.to_datetime({'year': tempdf['yr'], 'month': tempdf['mon'],
                                  'day': tempdf['day'], 'hour': tempdf['hh'],
                                  'minute': tempdf['mm'], 'second': tempdf['ss']},
                                 utc=True)
    timedelta_1950 = tempdtframe - pd.Timestamp('1950-01-01T00:00:00', tz='UTC')
    tempdf['TIME'] = timedelta_1950.dt.total_seconds() / 86400

    # Ignore the sample depths, and create a 5-m nominal depth series.
    # In the future, create a PRES variable and fill with nominal values if needed??
    tempdf['DEPH'] = np.ones(tempdf['sample_depth [m]'].__len__()) * 5.0

    # Rename latitude/longitude/expocode series
    tempdf['LATITUDE'] = tempdf['latitude [dec.deg.N]'].copy()
    tempdf['LONGITUDE'] = tempdf['longitude [dec.deg.E]'].copy()
    tempdf['EXPOCODE'] = tempdf['Expocode'].copy()

    # Remap existing flags to OceanSITES NaN -> 2->1 as int8 (equivalent to byte)
    tempdf.loc[tempdf['fCO2rec_flag'] == 2, 'fCO2rec_flag'] = 1
    tempdf['fCO2rec_flag'] = tempdf['fCO2rec_flag'].astype('int8')

    # Recreate flags for some variables (Dimension variables flags are created in the create_dimensions function)
    # Temperature and salinity are not QCed per se
    tempdf['SST_flag'] = np.zeros(tempdf.shape[0]).astype('int8')
    tempdf['sal_flag'] = np.zeros(tempdf.shape[0]).astype('int8')

    # Read info files
    SOCAT_info = pd.read_csv(os.path.join(input_files_dir, SOCAT_info_file), sep='\t',
                             skiprows=0, dtype='str', na_values='',
                             on_bad_lines='skip')

    # Create platform_code series: Callsign if exists, name if it doesn't
    SOCAT_info['PlatformCode'] = SOCAT_info['CallSign_WMO']
    SOCAT_info.loc[SOCAT_info['PlatformCode'].isna(), 'PlatformCode'] = \
        SOCAT_info.loc[SOCAT_info['PlatformCode'].isna(), 'Name']
    SOCAT_info['PlatformCode'] = SOCAT_info['PlatformCode'].str.replace(' ', '')

    # Extract only variables that will go to CMEMS
    variables_dict = CMEMSdict.generate_variables_dictionary(output_files_dir)
    all_in_CMEMS_variables = ['TIME', 'LATITUDE', 'LONGITUDE', 'DEPH', *variables_dict['flag'].keys(),
                              *variables_dict['flag'].values()]
    in_CMEMS_variables = list(set(all_in_CMEMS_variables) & set(tempdf.columns))
    in_CMEMS_variables.append('EXPOCODE')
    tempdf = tempdf.loc[:, in_CMEMS_variables]

    # Call the netCDF builder
    import CMEMS_build_nc as cmems_build
    cmems_build.build_nc(tempdf, SOCAT_info, output_files_dir)
