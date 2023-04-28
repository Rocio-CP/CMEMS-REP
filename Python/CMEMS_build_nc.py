def build_nc(dataframe, infoframe, output_files_dir):
    import os
    import numpy as np
    import pandas as pd
    import netCDF4
    import CMEMS_dictionaries as CMEMSdict

    unique_platform_codes = infoframe['PlatformCode'].unique()

    z_counter=0
    
    for pc in unique_platform_codes:
        # Create boolean filters for
        filter_expocodes = infoframe['PlatformCode'] == pc
        current_expocodes_info = infoframe.loc[filter_expocodes,]
        filter_expocodes_data = dataframe['EXPOCODE'].isin(current_expocodes_info['EXPOCODE'])

        # Extract sub-dataframe and order chronologically
        current_dataframe = dataframe.loc[filter_expocodes_data,]
        current_dataframe = current_dataframe.sort_values(by="TIME")

        # Check if there is only one platform type and name per platform
        platform_type = current_expocodes_info['PlatformType'].unique()
        platform_name = current_expocodes_info['Name'].unique()
        if platform_type.__len__() > 1:
            raise Exception("There is more than one PlatformType code for the platform " + pc)
# Some platforms HELLO SAILBOATS change name
        if platform_name.__len__() > 1:
            platform_name=";".join(platform_name)
#        if platform_name.__len__() > 1:
#            raise Exception("There is more than one Name for the platform " + pc)

        # Figure out dataset and platform type folder
        if 'glodap-obs' in output_files_dir:
            # All GLODAP data come from bottle samples
            platform_type_folder = 'BO'
            nc_filename = 'GL_PR_' + platform_type_folder + '_' + pc + '-GLODAPv22022.nc'

        elif 'socat-obs' in output_files_dir:
            if platform_type in ['31', '32', '37']:
                platform_type_folder = 'CO'
            elif platform_type in ['3B']:
                if 'SD' in platform_name.item():
                    platform_type_folder = 'SD'
                elif 'WG' in platform_name.item():
                    platform_type_folder = 'GL'
            elif platform_type in ['41']:
                platform_type_folder = 'MO'
            elif platform_type in ['42']:
                platform_type_folder = 'DB'
            nc_filename = 'GL_TS_' + platform_type_folder + '_' + pc + '-SOCATv2022.nc'

        nc_filepathname = output_files_dir + '/' + platform_type_folder + '/' + nc_filename

        # Create NetCDF file
        # Remove previous version of the file, if it exists
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
        create_dimensions(timedim_df['TIME'], timedim_df['LATITUDE'],
                                    timedim_df['LONGITUDE'], depth_matrix, nc)

        # Create variables
        variables_dict = CMEMSdict.generate_variables_dictionary(output_files_dir)

        for variable_ind, variable_name in enumerate(variables_dict['SDN'].keys()):

            # Variables
            value_table = pd.pivot_table(current_dataframe, index="TIME", columns="DEPH",
                                         values=variable_name, aggfunc="mean", dropna=False)
            value_table_nobs = pd.pivot_table(current_dataframe, index="TIME", columns="DEPH",
                                              values=variable_name, aggfunc="count", dropna=False)

            variable_attributes, QC_attributes = CMEMSdict.variable_attributes_dictionary(variable_name, output_files_dir,platform_type_folder)

            variable = nc.createVariable(variables_dict['SDN'][variable_name],
                                         variable_attributes[0]["datatype"],
                                         variable_attributes[0]["dimensions"],
                                         fill_value=variable_attributes[0]["_FillValue"])

            value_table[np.isnan(value_table)] = variable_attributes[0]["_FillValue"] * \
                                                 variable_attributes[1]['scale_factor']

            for key, value in variable_attributes[1].items():
                if value or key == 'add_offset': # Since add_offset = 0, it takes it as False
                    variable.setncattr(key, value)

            variable[:] = value_table

            # QC variables
            value_qc_table = pd.pivot_table(current_dataframe, index="TIME", columns="DEPH",
                                            values=variables_dict['flag'][variable_name], aggfunc="mean", dropna=False)
            value_qc_table[value_table_nobs > 1] = 8
            value_qc_table[value_qc_table.isna()] = 9
            variableqc = nc.createVariable(variables_dict['SDN'][variable_name] + "_QC",
                                           QC_attributes[0]["datatype"],
                                           QC_attributes[0]["dimensions"],
                                           fill_value=QC_attributes[0]["_FillValue"])

            for key, value in QC_attributes[1].items():
                variableqc.setncattr(key, value)

            variableqc[:] = value_qc_table

        # Global attributes
        global_attributes = CMEMSdict.global_attributes_dictionary(current_expocodes_info, current_dataframe, nc_filename)
        for key, value in global_attributes.items():
            nc.setncattr(key, value)

        nc.close()

        # Python can't store flag_values as bytes. Use ncatted to correct
        nc=netCDF4.Dataset(nc_filepathname, mode='r')
        qcvars=[v for v in nc.variables.keys() if v.__contains__('_QC')]
        nc.close()
        for qc_var in qcvars:
            os.system('ncatted -O -h -a flag_values,' + qc_var + ',m,b,"0, 1, 2, 3, 4, 5, 6, 7, 8, 9" ' + nc_filepathname)


        print(nc_filename + ' ' +
              round(os.path.getsize(nc_filepathname)/(1024*1024),2).__str__() +' MB')


        z_counter=z_counter+1


def create_dimensions(timeval, lat, lon, depth, nc):
    import numpy as np
    import json

    # Import dimensions attributes
    with open("CMEMS_INSTAC_metadata.json", "r") as template:
        all_attributes = json.load(template)
    template.close()
    dimension_attributes = all_attributes['dimatt'][0]
    QC_attributes = all_attributes['QC_varatt']

    # DIMENSIONS
    nc.createDimension('TIME', timeval.__len__())
    nc.createDimension('LONGITUDE', lon.__len__())
    nc.createDimension('LATITUDE', lat.__len__())
    nc.createDimension('POSITION', lon.__len__())
    nc.createDimension('DEPTH', depth.shape[1])

    # Dimension VARIABLES
    # Create dictionary (in Matlab it would be a structure)
    dimensionvar_dict = dict()
    dimensionvar_dict['TIME'] = timeval
    dimensionvar_dict['LONGITUDE'] = lon
    dimensionvar_dict['LATITUDE'] = lat
    dimensionvar_dict['DEPH'] = depth

    for variable_name in dimension_attributes.keys():
        if variable_name == 'DEPH':
            var = nc.createVariable(variable_name,
                                dimension_attributes[variable_name][0]["datatype"],
                                dimension_attributes[variable_name][0]["dimensions"],
                                fill_value=dimension_attributes[variable_name][0]["_FillValue"])
        else:
            var = nc.createVariable(variable_name,
                                dimension_attributes[variable_name][0]["datatype"],
                                dimension_attributes[variable_name][0]["dimensions"])

        var[:] = dimensionvar_dict[variable_name]

        for key, value in dimension_attributes[variable_name][1].items():
            var.setncattr(key, value)

    # Dimension FLAG variables. Assuming DEPH is nominal, and TIME/POSITION are good
    # Create dictionary (in Matlab it would be a structure)
    dimensionQC_dict = dict()
    dimensionQC_dict['TIME_QC'] = dict()
    dimensionQC_dict['POSITION_QC'] = dict()
    dimensionQC_dict['DEPH_QC'] = dict()

    dimensionQC_dict['TIME_QC']["values"] = np.ones(nc.dimensions['TIME'].size).astype('int8')
    dimensionQC_dict['TIME_QC']["dimensions"] = 'TIME'
    dimensionQC_dict['TIME_QC']["long_name"] = "Time quality flag"

    dimensionQC_dict['POSITION_QC']["values"] = np.ones(nc.dimensions['LONGITUDE'].size).astype('int8')
    dimensionQC_dict['POSITION_QC']["dimensions"] = 'POSITION'
    dimensionQC_dict['POSITION_QC']["long_name"] = "Position quality flag"

    dimensionQC_dict['DEPH_QC']["values"] = np.repeat(np.reshape(np.ones(nc.dimensions['DEPTH'].size) * 7, (1, -1)),
                                             nc.dimensions['TIME'].size, axis=0).astype('int8')  # Make matrix!!
    dimensionQC_dict['DEPH_QC']["dimensions"]=['TIME', 'DEPTH']
    dimensionQC_dict['DEPH_QC']["long_name"] = "Depth quality flag"

    for dimQC_name in dimensionQC_dict.keys():

        dimQC_var = nc.createVariable(dimQC_name,
                                      QC_attributes[0]["datatype"],
                                      dimensionQC_dict[dimQC_name]["dimensions"],
                                      fill_value=QC_attributes[0]["_FillValue"])

        dimQC_var[:] = dimensionQC_dict[dimQC_name]["values"]

        QC_attributes[1]["long_name"]=dimensionQC_dict[dimQC_name]["long_name"]
        QC_attributes[1]["valid_min"] = np.int8(QC_attributes[1]["valid_min"])
        QC_attributes[1]["valid_max"] = np.int8(QC_attributes[1]["valid_max"])

        for key, value in QC_attributes[1].items():
            dimQC_var.setncattr(key, value)

    return nc
