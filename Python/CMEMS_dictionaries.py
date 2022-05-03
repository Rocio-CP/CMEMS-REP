def common_dictionaries():
    # Some parameters / unchanging values / dictionaries
    dimensions=('G2datetime', 'G2depthnominal','G2latitude','G2longitude','')
    dimension_flags=('G2datetimef','G2depthnominalf','','','G2positionf')
    SDN_dimension_names=('TIME','DEPTH','LATITUDE','LONGITUDE','POSITION')
    SDN_dimension_variable_names=('TIME','DEPH','LATITUDE','LONGITUDE')
    dimension_dict = {}
    dimension_dict['flag'] = dict(zip(dimensions, dimension_flags))
    dimension_dict['dim_name'] = dict(zip(dimensions, SDN_dimension_names))
    dimension_dict['dim_var_name'] = dict(zip(dimensions, SDN_dimension_variable_names))

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

    return dimension_dict, variables_dict, G2OS_flag_dict

def global_attributes_dictionary(current_expocodes_info, current_dataframe):
    import datetime
    import json
    import numpy as np

    with open("CMEMS_INSTAC_metadata.json", "r") as template:
        all_attributes = json.load(template)
    template.close()

    all_attributes['globalatt'][0]['platform_code']= current_expocodes_info['PlatformCode'].unique().item() # only alphanumeric characters
    all_attributes['globalatt'][0]['platform_name']= current_expocodes_info['Name'].unique().item()
    all_attributes['globalatt'][0]['id']='GL_PR_BO_'+all_attributes['globalatt'][0]['platform_code']+'-GLODAPv22022.nc'
    all_attributes['globalatt'][0]['wmo_platform_code']=current_expocodes_info['CallSign_WMO'].unique().item()
    all_attributes['globalatt'][0]['ices_platform_code']=current_expocodes_info['ICEScode'].unique().item()

    all_attributes['globalatt'][0]['source_platform_category_code'] = str(current_expocodes_info['PlatformType'].unique().item())
    if all_attributes['globalatt'][0]['source_platform_category_code'] == '31':
        all_attributes['globalatt'][0]['source'] = 'research vessel'
    #    output_file_full_name = output_files_dir + 'VESSEL/' + all_attributes['globalatt'][0]['id']
    elif all_attributes['globalatt'][0]['source_platform_category_code'] == '21':
        all_attributes['globalatt'][0]['source'] = 'propelled manned submersible'
    #    output_file_full_name = output_files_dir + 'ETC/' + all_attributes['globalatt'][0]['id']
    elif all_attributes['globalatt'][0]['source_platform_category_code'] == '62':
        all_attributes['globalatt'][0]['source'] = 'aeroplane'
    #    output_file_full_name = output_files_dir + 'ETC/' + all_attributes['globalatt'][0]['id']

    institution_name=current_expocodes_info['Institution'].unique()
    pi_institution_name=current_expocodes_info['PI_Institution'].unique()
    edmo=current_expocodes_info['EDMO'].unique()
    pi_edmo=current_expocodes_info['PI_EDMO'].unique()
    all_attributes['globalatt'][0]['institution'] =";".join(np.unique([*institution_name, *pi_institution_name]))
    all_attributes['globalatt'][0]['institution_edmo_code'] = ",".join(np.unique([*edmo, *pi_edmo]))

    # Geospatial and time limits
    all_attributes['globalatt'][0]['geospatial_lat_min'] =current_dataframe['G2latitude'].min()
    all_attributes['globalatt'][0]['geospatial_lat_max'] =current_dataframe['G2latitude'].max().__str__()
    all_attributes['globalatt'][0]['geospatial_lon_min'] =current_dataframe['G2longitude'].min().__str__()
    all_attributes['globalatt'][0]['geospatial_lon_max'] =current_dataframe['G2longitude'].max().__str__()
    all_attributes['globalatt'][0]['geospatial_vertical_min'] =current_dataframe['G2depthnominal'].min().__str__()
    all_attributes['globalatt'][0]['geospatial_vertical_max'] =current_dataframe['G2depthnominal'].max().__str__()
    # CORRECT THE DATE FORMATS!!
    #all_attributes['globalatt'][0]['time_coverage_start'] =current_dataframe['G2datetime'].iloc[0].__str__()
    #all_attributes['globalatt'][0]['time_coverage_end'] =current_dataframe['G2datetime'].iloc[-1].__str__()

    # PI attributes
    pi=current_expocodes_info['PI'].unique()
    all_attributes['globalatt'][0]['pi_name'] =";".join(pi)

    all_attributes['globalatt'][0]['date_update'] =datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%SZ")
    all_attributes['globalatt'][0]['history'] = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%SZ")+ ": Creation"

    all_attributes['globalatt'][0]['last_date_observation'] = current_dataframe['G2datetime'].iloc[-1]
    all_attributes['globalatt'][0]['last_latitude_observation'] = current_dataframe['G2latitude'].iloc[-1]
    all_attributes['globalatt'][0]['last_longitude_observation'] = current_dataframe['G2longitude'].iloc[-1]


    return all_attributes['globalatt'][0]

def dimension_attributes_dictionary(dimension_variable_name):
    import json
    with open("CMEMS_INSTAC_metadata.json", "r") as template:
        all_attributes = json.load(template)
    template.close()
    dimension_attributes = all_attributes['dimatt'][0][dimension_variable_name]
    return dimension_attributes

def variable_attributes_dictionary(variable_name):
    import json
    _,variables_dict,_ = common_dictionaries()
    with open("CMEMS_INSTAC_metadata.json", "r") as template:
        all_attributes = json.load(template)
    template.close()

    all_attributes['variable_varatt'][1]["standard_name"] = variables_dict['CF'][variable_name]
    all_attributes['variable_varatt'][1]["units"] = variables_dict['unit'][variable_name]
    all_attributes['variable_varatt'][1]["long_name"] = variables_dict['long'][variable_name]
    all_attributes['variable_varatt'][1]["ancillary_variables"] = variables_dict['SDN'][variable_name] + "_QC"

    all_attributes['QC_varatt'][1]["long_name"] = variables_dict['long'][variable_name] + " quality flag"

    return all_attributes['variable_varatt'], all_attributes['QC_varatt']

