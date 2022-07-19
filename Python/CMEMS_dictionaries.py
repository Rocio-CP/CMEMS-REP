def common_dictionaries(data_source):
    # Some parameters / unchanging values / dictionaries
    socat_input_vars = ('', '', 'SST', 'sal', '', '', '', '', '', '', '', '', '', '', '', '', '', 'fCO2rec')

    glodap_input_vars = ('G2bottomdepth', 'G2pressure', 'G2temperature', 'G2salinity', 'G2oxygen', 'G2nitrate',
                         'G2nitrite', 'G2phosphate', 'G2silicate', 'G2phtsinsitutp', 'G2phts25p0', 'G2tco2',
                         'G2talk', 'G2doc', 'G2don', 'G2tdn', 'G2chla', '')
    glodap_input_flag_vars = (
    'G2bottomdepthf', 'G2pressuref', 'G2temperaturef', 'G2salinityf', 'G2oxygenf', 'G2nitratef',
    'G2nitritef', 'G2phosphatef', 'G2silicatef', 'G2phtsinsitutpf', 'G2phts25p0f', 'G2tco2f',
    'G2talkf', 'G2docf', 'G2donf', 'G2tdnf', 'G2chlaf', '')
    SDN_var_names = ('BATH', 'PRES', 'TEMP', 'PSAL', 'DOX2', 'NTAW',
                     'NTIW', 'PHOW', 'SLCW', 'PHPH', 'PH25', 'TICW',
                     'ALKW', 'CORG', 'NODW', 'NT1D', 'CPHL', 'FCO2')
    units = ('m', 'dbar', 'degrees_C', '0.001', 'µmol kg-1', 'µmol kg-1',
             'µmol kg-1', 'µmol kg-1', 'µmol kg-1', '1', '1', 'µmol kg-1',
             'µmol kg-1', 'µmol kg-1', 'µmol kg-1', 'µmol kg-1', 'mg m-3', 'µatm')
    long_name = ('Bathymetric depth', 'Sea pressure',
                 'Sea temperature', 'Practical salinity',
                 'Dissolved oxygen', 'Nitrate (NO3-N)',
                 'Nitrite (NO2-N)', 'Phosphate (PO4-P)', 'Silicate (SIO4-SI)',
                 'Ph', 'Ph at 25 °C and 0 dbar', 'Dissolved inorganic carbon',
                 'Total alkalinity', 'Dissolved organic carbon', 'Dissolved organic nitrogen',
                 'Total dissolved nitrogen', 'Chlorophyll-a', 'CO2 fugacity')
    CF_standard_name = ('sea_floor_depth_below_sea_surface', 'sea_water_pressure', 'sea_water_temperature',
                        'sea_water_practical_salinity', 'moles_of_oxygen_per_unit_mass_in_sea_water',
                        'moles_of_nitrate_per_unit_mass_in_sea_water', 'moles_of_nitrite_per_unit_mass_in_sea_water',
                        'moles_of_phosphate_per_unit_mass_in_sea_water', 'moles_of_silicate_per_unit_mass_in_sea_water',
                        'sea_water_ph_reported_on_total_scale', ' ',
                        ' ', 'sea_water_alkalinity_per_unit_mass',
                        'moles_of_dissolved_organic_carbon_per_unit_mass_in_sea_water', ' ',
                        'moles_of_dissolved_total_nitrogen_per_unit_mass_in_sea_water',
                        'mass_concentration_of_chlorophyll_a_in_sea_water', 'fugacity_of_carbon_dioxide_in_sea_water')
    #valid_min =()
    #valid_max=()
    if 'SOCAT' in data_source:
        variables_dict = {}
        variables_dict['SDN'] = dict(zip(socat_input_vars, SDN_var_names))
        variables_dict['unit'] = dict(zip(socat_input_vars, units))
        variables_dict['long'] = dict(zip(socat_input_vars, long_name))
        variables_dict['CF'] = dict(zip(socat_input_vars, CF_standard_name))

    elif 'GLODAP' in data_source:
        variables_dict = {}
        variables_dict['flag'] = dict(zip(glodap_input_vars, glodap_input_flag_vars))
        variables_dict['SDN'] = dict(zip(glodap_input_vars, SDN_var_names))
        variables_dict['unit'] = dict(zip(glodap_input_vars, units))
        variables_dict['long'] = dict(zip(glodap_input_vars, long_name))
        variables_dict['CF'] = dict(zip(glodap_input_vars, CF_standard_name))
        GLODAP_flags = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        OceanSITES_flags = [8, 9, 1, 2, 4, 9, 8, 0, 0, 9]
        variables_dict['flag_values'] = dict(zip(GLODAP_flags, OceanSITES_flags))

    for keyval in variables_dict:
        if '' in variables_dict[keyval]:
            variables_dict[keyval].pop('')

    return variables_dict

def create_dimensions(timeval, lat, lon, depth, nc):
    import numpy as np

    # Import ATTRIBUTES
    import json
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
                                fill_value=dimension_attributes[variable_name][0]["_FillValue"]    )
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

    dimensionQC_dict['TIME_QC'][ "values"] = np.ones(nc.dimensions['TIME'].size)
    dimensionQC_dict['TIME_QC']["dimensions"] = 'TIME'
    dimensionQC_dict['TIME_QC']["long_name"] = "Time quality flag"

    dimensionQC_dict['POSITION_QC']["values"] = np.ones(nc.dimensions['LONGITUDE'].size)
    dimensionQC_dict['POSITION_QC']["dimensions"] = 'POSITION'
    dimensionQC_dict['POSITION_QC']["long_name"] = "Position quality flag"

    dimensionQC_dict['DEPH_QC']["values"] = np.repeat(np.reshape(np.ones(nc.dimensions['DEPTH'].size) * 7, (1, -1)),
                                             nc.dimensions['TIME'].size, axis=0)  # Make matrix!!
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

def variable_attributes_dictionary(variable_name, data_source):
    import json
    import numpy as np

    variables_dict = common_dictionaries(data_source)
    with open("CMEMS_INSTAC_metadata.json", "r") as template:
        all_attributes = json.load(template)
    template.close()

    variable_attributes=all_attributes['variable_varatt']
    QC_attributes=all_attributes['QC_varatt']

    variable_attributes[1]["standard_name"] = variables_dict['CF'][variable_name]
    variable_attributes[1]["units"] = variables_dict['unit'][variable_name]
    variable_attributes[1]["long_name"] = variables_dict['long'][variable_name]
    variable_attributes[1]["ancillary_variables"] = variables_dict['SDN'][variable_name] + "_QC"
    #variable_attributes[1]["valid_min"] = np.int32(variable_attributes[1]["valid_min"])
    #variable_attributes[1]["valid_max"] = np.int32(variable_attributes[1]["valid_max"])

    QC_attributes[1]["long_name"] = variables_dict['long'][variable_name] + " quality flag"
    QC_attributes[1]["valid_min"] = np.int8(QC_attributes[1]["valid_min"])
    QC_attributes[1]["valid_max"] = np.int8(QC_attributes[1]["valid_max"])

    return variable_attributes, QC_attributes


def global_attributes_dictionary(current_expocodes_info, current_dataframe):
    import datetime
    import json
    import numpy as np

    with open("CMEMS_INSTAC_metadata.json", "r") as template:
        all_attributes = json.load(template)
    template.close()

    all_attributes['globalatt'][0]['platform_code'] = current_expocodes_info[
        'PlatformCode'].unique().item()  # only alphanumeric characters
    all_attributes['globalatt'][0]['platform_name'] = current_expocodes_info['Name'].unique().item()
    all_attributes['globalatt'][0]['id'] = 'GL_PR_BO_' + all_attributes['globalatt'][0][
        'platform_code'] + '-GLODAPv22022.nc'
    all_attributes['globalatt'][0]['wmo_platform_code'] = current_expocodes_info['CallSign_WMO'].unique().item()
    all_attributes['globalatt'][0]['ices_platform_code'] = current_expocodes_info['ICEScode'].unique().item()

    all_attributes['globalatt'][0]['source_platform_category_code'] = str(
        current_expocodes_info['PlatformType'].unique().item())
    if all_attributes['globalatt'][0]['source_platform_category_code'] == '31':
        all_attributes['globalatt'][0]['source'] = 'research vessel'
    elif all_attributes['globalatt'][0]['source_platform_category_code'] == '32':
        all_attributes['globalatt'][0]['source'] = 'vessel of opportunity'
    elif all_attributes['globalatt'][0]['source_platform_category_code'] == '37':
        all_attributes['globalatt'][0]['source'] = 'self-propelled boat'
    elif all_attributes['globalatt'][0]['source_platform_category_code'] == '3B':
        all_attributes['globalatt'][0]['source'] = 'autonomous surface water vehicle'
    elif all_attributes['globalatt'][0]['source_platform_category_code'] == '41':
        all_attributes['globalatt'][0]['source'] = 'moored surface buoy'
    elif all_attributes['globalatt'][0]['source_platform_category_code'] == '42':
        all_attributes['globalatt'][0]['source'] = 'drifting surface float'
    elif all_attributes['globalatt'][0]['source_platform_category_code'] == '21':
        all_attributes['globalatt'][0]['source'] = 'propelled manned submersible'
    elif all_attributes['globalatt'][0]['source_platform_category_code'] == '62':
        all_attributes['globalatt'][0]['source'] = 'aeroplane'

    institution_name = current_expocodes_info['Institution'].unique()
    pi_institution_name = current_expocodes_info['PI_Institution'].unique()
    edmo = current_expocodes_info['EDMO'].unique()
    pi_edmo = current_expocodes_info['PI_EDMO'].unique()
    all_attributes['globalatt'][0]['institution'] = ";".join(np.unique([*institution_name, *pi_institution_name]))
    all_attributes['globalatt'][0]['institution_edmo_code'] = ",".join(np.unique([*edmo, *pi_edmo]))

    # Geospatial and time limits
    all_attributes['globalatt'][0]['geospatial_lat_min'] = current_dataframe['LATITUDE'].min()
    all_attributes['globalatt'][0]['geospatial_lat_max'] = current_dataframe['LATITUDE'].max()
    all_attributes['globalatt'][0]['geospatial_lon_min'] = current_dataframe['LONGITUDE'].min()
    all_attributes['globalatt'][0]['geospatial_lon_max'] = current_dataframe['LONGITUDE'].max()
    all_attributes['globalatt'][0]['geospatial_vertical_min'] = current_dataframe['DEPH'].min().__str__()
    all_attributes['globalatt'][0]['geospatial_vertical_max'] = current_dataframe['DEPH'].max().__str__()
    all_attributes['globalatt'][0]['time_coverage_start'] =date19502string(current_dataframe['TIME'].iloc[0])
    all_attributes['globalatt'][0]['time_coverage_end'] =date19502string(current_dataframe['TIME'].iloc[-1])

    # PI attributes
    pi = current_expocodes_info['PI'].unique()
    all_attributes['globalatt'][0]['pi_name'] = ";".join(pi)

    all_attributes['globalatt'][0]['date_update'] = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%SZ")
    all_attributes['globalatt'][0]['history'] = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%SZ") + ": Creation"

    all_attributes['globalatt'][0]['last_date_observation'] = date19502string(current_dataframe['TIME'].iloc[-1])
    all_attributes['globalatt'][0]['last_latitude_observation'] = current_dataframe['LATITUDE'].iloc[-1]
    all_attributes['globalatt'][0]['last_longitude_observation'] = current_dataframe['LONGITUDE'].iloc[-1]

    return all_attributes['globalatt'][0]


def date19502string(datenumber):
    import datetime
    base_date=datetime.datetime(1950,1,1,0,0,0)
    datedelta=datetime.timedelta(datenumber)
    datestring=(datedelta + base_date).strftime("%Y-%m-%dT%H:%M:%SZ")
    return(datestring)
