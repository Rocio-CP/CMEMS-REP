def read_GLODAP_obs(input_files_dir, output_files_dir, GLODAP_files, GLODAP_info_file):

    import itertools
    import os
    import numpy as np
    import pandas as pd
    import Python.CMEMS_dictionaries as CMEMSdict

    # Check if input data files are there; if not, fetch data from the internet

    ### Read GLODAP data file(s)
    # # It may not be necessary to specify the datatype for reading: it'll transform when writing the netCDF file
    # all_glodap_input_variables = pd.read_csv(os.path.join(input_files_dir, sourcefile),
    #                                          index_col=0, nrows=0).columns.tolist()
    #
    # # pd.Int16Dtype() is the type of pandas integer that allows for na
    # all_glodap_input_variables_dtype = ['UInt64'] * all_glodap_input_variables.__len__()
    # for n, gi in enumerate(all_glodap_input_variables):
    #     if not gi.endswith('f') and 'qc' not in gi and gi not in ['G2cruise', 'G2region', 'G2cast', 'G2year', 'G2month',
    #                                                               'G2day', 'G2hour', 'G2minute', 'G2bottle']:
    #         all_glodap_input_variables_dtype[n] = 'float'
    #
    # # dtype for GLODAP data
    # data_dtype_dict = dict(zip(all_glodap_input_variables, all_glodap_input_variables_dtype))

    if any(['DOI' in i for i in GLODAP_files]):
        sourcefile = [i for i in GLODAP_files if 'GLODAP' in i][0]
        doifile = [i for i in GLODAP_files if 'DOI' in i][0]
        expocodefile = [i for i in GLODAP_files if 'EXPOCODES' in i][0]

        # Read GLODAP file (parse_dates does not work because some hours and minutes are nan)
        tempdf = pd.read_csv(os.path.join(input_files_dir, sourcefile), sep=',',
                             skiprows=0, na_values='-9999',  # dtype=data_dtype_dict,
                             on_bad_lines='skip')

        # Read extra info files
        dois = pd.read_csv(os.path.join(input_files_dir, doifile), sep='\t', header=None,
                          names=['G2cruise', 'doi'], dtype={'G2cruise':'int','doi':'str'}, on_bad_lines='skip',
                          encoding='utf_16_le')
        doisdict = dois.set_index('G2cruise')['doi'].to_dict()

        expocodes = pd.read_csv(os.path.join(input_files_dir, expocodefile), sep='\t',
                                    header=None, names=['G2cruise', 'expocode'],
                                    dtype={'G2cruise': 'int', 'expocode': 'str'},
                                    on_bad_lines='skip')

            ### Assign Expocodes and DOIs to data frame
            # Transform expocode dataframes into dictionaries (easier to lookup)
        expocodesdict = expocodes.set_index('G2cruise')['expocode'].to_dict()

        for cruise in tempdf['G2cruise'].unique():
                tempdf.loc[tempdf.G2cruise == cruise,
                           'EXPOCODE'] = expocodesdict[cruise]
                tempdf.loc[tempdf.G2cruise == cruise,
                    'G2doi'] = doisdict[cruise]


    else:
        sourcefile=GLODAP_files[0]
        # Read GLODAP file (parse_dates does not work because some hours and minutes are nan)
        tempdf = pd.read_csv(os.path.join(input_files_dir, sourcefile), sep=',',
                             skiprows=0, na_values='-9999',  # dtype=data_dtype_dict,
                             on_bad_lines='skip')
        tempdf['EXPOCODE'] = tempdf['G2expocode'].copy()

    variables_dict = CMEMSdict.generate_variables_dictionary(output_files_dir)

    # Create datetime from 1950 series
    # Change NaN hour/minute to 0. Calculate 1950-time
    tempdf.loc[np.isnan(tempdf['G2hour']), 'G2hour'] = 0
    tempdf.loc[np.isnan(tempdf['G2minute']), 'G2minute'] = 0
    # Create days-from-1950 numeric datetime
    tempdtframe = pd.to_datetime({'year': tempdf['G2year'], 'month': tempdf['G2month'],
                                  'day': tempdf['G2day'], 'hour': tempdf['G2hour'],
                                  'minute': tempdf['G2minute']},
                                 utc=True)
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
    tempdf['LATITUDE'] = tempdf['G2latitude'].copy()
    tempdf['LONGITUDE'] = tempdf['G2longitude'].copy()

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

    # Recreate flags for some variables (Dimension variables flags are created in the create_dimensions function)
    # If temperature, pressure, and bottom depth exist, they're assumed good (flag 1)
    tempdf['G2bottomdepthf'] = np.ones(tempdf['G2bottomdepth'].__len__())
    tempdf.loc[np.isnan(tempdf['G2bottomdepth']), 'G2bottomdepthf'] = 9
    tempdf['G2bottomdepthf']=tempdf['G2bottomdepthf'].astype('int')
    tempdf['G2temperaturef'] = np.ones(tempdf['G2temperature'].__len__())
    tempdf.loc[np.isnan(tempdf['G2temperature']), 'G2temperaturef'] = 9
    tempdf['G2temperaturef']=tempdf['G2temperaturef'].astype('int')
    tempdf['G2pressuref'] = np.ones(tempdf['G2pressure'].__len__())
    tempdf.loc[np.isnan(tempdf['G2pressure']), 'G2pressuref'] = 9
    tempdf['G2pressuref']=tempdf['G2pressuref'].astype('int')

    # Rename G2fco2 to G2fco2_20_0 (in GLODAP, it's given at 20 dg, 0dbar).
    # tempdf.rename(columns={'G2fco2': 'G2fco2_20_0'}, inplace=True)

    ### Read GLODAP info: platform names, codes, EDMO, PI...
    # # Dtype for GLODAP info
    glodap_info_dtype = ['int', *['str'] * 12]
    info_dtype_dict = dict(zip(range(0, 13), glodap_info_dtype))
    GLODAP_info = pd.read_csv(os.path.join(input_files_dir, GLODAP_info_file), sep='\t',
                              skiprows=0, dtype=info_dtype_dict, na_values='',
                              on_bad_lines='skip')

    # Create platform_code series: Callsign if exists, name if it doesn't
    GLODAP_info['PlatformCode'] = GLODAP_info['CallSign_WMO']
    GLODAP_info.loc[GLODAP_info['PlatformCode'].isna(), 'PlatformCode'] = \
        GLODAP_info.loc[GLODAP_info['PlatformCode'].isna(), 'Name']
    GLODAP_info['PlatformCode'] = GLODAP_info['PlatformCode'].str.replace(' ', '')

    # Extract only variables that will go to CMEMS
    all_in_CMEMS_variables = ['TIME', 'LATITUDE', 'LONGITUDE', 'DEPH', *variables_dict['flag'].keys(),
                              *variables_dict['flag'].values()]
    in_CMEMS_variables = list(set(all_in_CMEMS_variables) & set(tempdf.columns))
    in_CMEMS_variables.append('EXPOCODE')
    tempdf=tempdf.loc[ : , in_CMEMS_variables]

    # Call the netCDF builder
    import CMEMS_build_nc as cmems_build
    print(tempdf.columns)
    cmems_build.build_nc(tempdf,GLODAP_info,output_files_dir)
