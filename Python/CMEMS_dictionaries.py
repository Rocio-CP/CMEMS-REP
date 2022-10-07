def generate_variables_dictionary(output_files_dir):
    socat_input_vars = ('', '', 'SST [deg.C]', 'sal', '', '', '', '', '', '', '', '', '', '', '', '', '', 'fCO2rec [uatm]')
    socat_input_flag_vars = ('', '', 'SST_flag', 'sal_flag', '', '', '', '', '', '', '', '', '', '', '', '', '', 'fCO2rec_flag')

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
    if 'socat-obs' in output_files_dir:
        variables_dict = {}
        variables_dict['flag'] = dict(zip(socat_input_vars, socat_input_flag_vars))
        variables_dict['SDN'] = dict(zip(socat_input_vars, SDN_var_names))
        variables_dict['unit'] = dict(zip(socat_input_vars, units))
        variables_dict['long'] = dict(zip(socat_input_vars, long_name))
        variables_dict['CF'] = dict(zip(socat_input_vars, CF_standard_name))

    elif 'glodap-obs' in output_files_dir:
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


def variable_attributes_dictionary(variable_name, output_files_dir, platform_type):
    import json
    import numpy as np

    variables_dict = generate_variables_dictionary(output_files_dir)
    with open("CMEMS_INSTAC_metadata.json", "r") as template:
        all_attributes = json.load(template)
    template.close()

    variable_attributes=all_attributes['variable_varatt']
    QC_attributes=all_attributes['QC_varatt']

    variable_attributes[1]["standard_name"] = variables_dict['CF'][variable_name]
    variable_attributes[1]["units"] = variables_dict['unit'][variable_name]
    variable_attributes[1]["long_name"] = variables_dict['long'][variable_name]
    variable_attributes[1]["ancillary_variables"] = variables_dict['SDN'][variable_name] + "_QC"

    if 'glodap-obs' in output_files_dir:
        variable_attributes[1]["sensor_mount"] = "mounted_on_shipborne_profiler"
    elif 'socat-obs' in output_files_dir:
        if platform_type in ['31','32','37']:
            variable_attributes[1]["sensor_mount"] = "mounted_on_shipborne_fixed"
        elif platform_type in ['3B','42']:
            variable_attributes[1]["sensor_mount"] = "mounted_on_glider"
        elif platform_type in ['41']:
            variable_attributes[1]["sensor_mount"] = "mounted_on_surface_buoy"

    #variable_attributes[1]["valid_min"] = np.int32(variable_attributes[1]["valid_min"])
    #variable_attributes[1]["valid_max"] = np.int32(variable_attributes[1]["valid_max"])

    QC_attributes[1]["long_name"] = variables_dict['long'][variable_name] + " quality flag"
    QC_attributes[1]["valid_min"] = np.int8(QC_attributes[1]["valid_min"])
    QC_attributes[1]["valid_max"] = np.int8(QC_attributes[1]["valid_max"])

    return variable_attributes, QC_attributes


def global_attributes_dictionary(current_expocodes_info, current_dataframe, nc_filename):
    import datetime
    import json
    import numpy as np
    import pandas as pd

    with open("CMEMS_INSTAC_metadata.json", "r") as template:
        all_attributes = json.load(template)
    template.close()

# Add eventually: a check of the mandatory value and optional value global attributes
    mandatory_value_att=['platform_code', 'data_mode','title','naming_authority','id','source', 'source_platform_category_code','institution','contact',
                   'geospatial_lat_min','geospatial_lat_max','geospatial_lon_min','geospatial_lon_max','geospatial_vertical_min','geospatial_vertical_max',
                   'time_coverage_start','time_coverage_end','cdm_data_type','data_type',
                   'format_version', 'Conventions', 'netcdf_version', 'references', 'data_assembly_center', 'update_interval','citation', 'distribution_statement'
                   'date_update', 'last_date_observation','last_latitude_observation','last_longitude_observation']
    optional_value_att=['platform_name', 'summary','wmo_platform_code','ices_platform_code','institution_edmo_code','institution_references','site_code','comment',
                        'area','bottom_depth','doi','pi_name','qc_manual','history','wmo_inst_type']
    additional_att=[]

    # Platform attributes
    # Some global attributes may have multiple values (e.g. multiple names, ICES codes, etc...)
    # Platform codes ARE UNIQUE
    # Platform names ALWAYS EXIST
    # ICES codes may not exist
    globattribute_infodf_dict={'platform_code':'PlatformCode',
                               'platform_name':'Name',
                               'ices_platform_code': 'ICEScode'}

    #if current_expocodes_info['Name'].unique().__len__() > 1:
    all_attributes['globalatt'][0]['platform_name'] = ";".join(
                current_expocodes_info['Name'].unique())
    if any(pd.isnull(current_expocodes_info['ICEScode'].unique())):
        all_attributes['globalatt'][0]['ices_platform_code'] = ""
    else:
        all_attributes['globalatt'][0]['ices_platform_code'] = " ".join(
                current_expocodes_info['ICEScode'].unique())

    all_attributes['globalatt'][0]['platform_code'] = current_expocodes_info['PlatformCode'].unique().item()
    all_attributes['globalatt'][0]['wmo_platform_code'] = current_expocodes_info['CallSign_WMO'].unique().item()
    all_attributes['globalatt'][0]['id']= nc_filename[0:-3]

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

    # Geospatial and time limits
    all_attributes['globalatt'][0]['geospatial_lat_min'] = current_dataframe['LATITUDE'].min()
    all_attributes['globalatt'][0]['geospatial_lat_max'] = current_dataframe['LATITUDE'].max()
    all_attributes['globalatt'][0]['geospatial_lon_min'] = current_dataframe['LONGITUDE'].min()
    all_attributes['globalatt'][0]['geospatial_lon_max'] = current_dataframe['LONGITUDE'].max()
    all_attributes['globalatt'][0]['geospatial_vertical_min'] = current_dataframe['DEPH'].min().__str__()
    all_attributes['globalatt'][0]['geospatial_vertical_max'] = current_dataframe['DEPH'].max().__str__()
    all_attributes['globalatt'][0]['time_coverage_start'] =date19502string(current_dataframe['TIME'].iloc[0])
    all_attributes['globalatt'][0]['time_coverage_end'] =date19502string(current_dataframe['TIME'].iloc[-1])

    all_attributes['globalatt'][0]['last_date_observation'] = date19502string(current_dataframe['TIME'].iloc[-1])
    all_attributes['globalatt'][0]['last_latitude_observation'] = current_dataframe['LATITUDE'].iloc[-1]
    all_attributes['globalatt'][0]['last_longitude_observation'] = current_dataframe['LONGITUDE'].iloc[-1]

    # institution and PI attributes
    # There may be multiple PIs associated with single expocodes that would appear repeated if plain unique() was used
    # join and then re-split by ';'
    allpis=';'.join(current_expocodes_info['PI'].unique()).split(';')
    # remove trailing/leading whitespaces and run unique on them
    pis=np.unique([dum.strip() for dum in allpis])
    all_attributes['globalatt'][0]['pi_name'] = ";".join(pis)

    institution_name = current_expocodes_info['Institution'].unique()
    pi_institution_name = current_expocodes_info['PI_Institution'].unique()
    edmo = current_expocodes_info['EDMO'].unique()
    pi_edmo = current_expocodes_info['PI_EDMO'].unique()
    all_attributes['globalatt'][0]['institution'] = ";".join(np.unique([*institution_name, *pi_institution_name]))
    if any(pd.isnull(np.unique([*edmo, *pi_edmo]))):
        all_attributes['globalatt'][0]['institution_edmo_code'] = ""
    else:
        all_attributes['globalatt'][0]['institution_edmo_code'] = " ".join(np.unique([*edmo, *pi_edmo]))

    # History attributes
    all_attributes['globalatt'][0]['date_update'] = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%SZ")
    all_attributes['globalatt'][0]['history'] = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%SZ") + ": Creation"

    # Data product-specific attributes
    if 'socat'.casefold() in nc_filename.casefold():
        all_attributes['globalatt'][0]["title"] = all_attributes['globalatt'][0]["title"] + " - SOCATv2022"
        all_attributes['globalatt'][0]["references"] = all_attributes['globalatt'][0]["references"] + " https://socat.info"
        all_attributes['globalatt'][0]["citation"] = all_attributes['globalatt'][0]["citation"] + "  SOCAT is described in Bakker et al. (2016) https://doi.org/10.5194/essd-8-383-2016; traceable citations are essential for justifying and sustaining the effort."
        all_attributes['globalatt'][0]["doi"] = "https://doi.org/10.25921/1h9f-nb73"

        if all_attributes['globalatt'][0]['source_platform_category_code'] in ['31','32','37','3B','42']:
            all_attributes['globalatt'][0]["cdm_data_type"] = "trajectory"
            all_attributes['globalatt'][0]["data_type"] = "OceanSITES trajectory data"

        elif all_attributes['globalatt'][0]['source_platform_category_code'] in ['41']:
            all_attributes['globalatt'][0]["cdm_data_type"] = "timeSeries"
            all_attributes['globalatt'][0]["data_type"] = "OceanSITES time-series data"

    elif 'glodap'.casefold() in nc_filename.casefold():
        all_attributes['globalatt'][0]["title"] = all_attributes['globalatt'][0]["title"] + " - GLODAPv2.2022"
        all_attributes['globalatt'][0]["references"] = all_attributes['globalatt'][0]["references"] + " https://glodap.info"
        all_attributes['globalatt'][0]["citation"] = all_attributes['globalatt'][0]["citation"] + "  GLODAPv2.2022 is described in Lauvset et al. (2022) DOI MISSING; traceable citations are essential for justifying and sustaining the effort."
        all_attributes['globalatt'][0]["doi"] = "https://doi.org/10.25921/1f4w-0t92"

        all_attributes['globalatt'][0]["cdm_data_type"] = "profile"
        all_attributes['globalatt'][0]["data_type"] = "OceanSITES vertical profile"

    # Add the DOIS of the individual cruises to the "references" global attribute (as doi.org/---)
    doi_list = current_expocodes_info['DOI_link']
    doi_list_nonan = doi_list.loc[doi_list.notna()]
    if not doi_list_nonan.empty:
        indiv_dois=' '.join(doi_list_nonan.unique())
        all_attributes['globalatt'][0]["references"]=all_attributes['globalatt'][0]["references"] + ' ' + indiv_dois

    # Remove NaN attributes; substitute with empty string


    return all_attributes['globalatt'][0]


def date19502string(datenumber):
    import datetime
    base_date=datetime.datetime(1950,1,1,0,0,0)
    datedelta=datetime.timedelta(datenumber)
    datestring=(datedelta + base_date).strftime("%Y-%m-%dT%H:%M:%SZ")
    return(datestring)
