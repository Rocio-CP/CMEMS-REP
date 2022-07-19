def build_nc(dataframe, infoframe, output_files_dir):
    import os
    import numpy as np
    import pandas as pd
    import netCDF4
    import CMEMS_dictionaries as CMEMSdict

    # Loop through the platforms (the basis of INSTAC files)
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
                platform_type_folder= 'CO/' + global_attributes['id']
            elif platform_type in ['3B']:
                if 'SD' in platform_name:
                    platform_type_folder= 'GL/'
                elif 'WG' in platform_name:
                    platform_type_folder = 'SD/'
            elif platform_type in ['41']:
                platform_type_folder= 'MO/'
            elif platform_type in ['42']:
                platform_type_folder= 'DB/'

            nc_filename='GL_PR_' + platform_type_folder + '_' + pc + '-SOCATv2022.nc'

        nc_filepathname = output_files_dir + platform_type_folder + nc_filename


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
                                         values=variable_name, aggfunc="mean")
            value_table_nobs = pd.pivot_table(current_dataframe, index="TIME", columns="DEPH",
                                              values=variable_name, aggfunc="count")
            value_table[value_table.isna()] = netCDF4.default_fillvals['i4']

            variable_attributes, QC_attributes = CMEMSdict.variable_attributes_dictionary(variable_name, 'GLODAP')

            variable = nc.createVariable(variables_dict['SDN'][variable_name],
                                         variable_attributes[0]["datatype"],
                                         variable_attributes[0]["dimensions"],
                                         fill_value=variable_attributes[0]["_FillValue"])
            variable[:] = value_table * 1000

            for key, value in variable_attributes[1].items():
                if value:
                    variable.setncattr(key, value)

            # QC variables
            value_qc_table = pd.pivot_table(current_dataframe, index="TIME", columns="DEPH",
                                            values=variables_dict['flag'][variable_name], aggfunc="mean")
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
