# Some parameters / unchanging values / dictionaries
dimensions=('G2datetime', 'G2depthnominal','G2latitude','G2longitude')
dimension_flags=('G2datetimef','G2depthnominalf','G2positionf','G2positionf')
SDN_dimension_names=('TIME','DEPH','LATITUDE','LONGITUDE')
SDN_dimension_variable_names=('TIME','DEPTH','LATITUDE','LONGITUDE')
dimension_units=("days since 1950-01-01T00:00:00Z","m","degree_north","degree_east")
dimension_long_name=("Time","Depth","Latitude of each location","Longitude of each location")
dimension_CF_standard_name=("time","depth","latitude","longitude")
dimension_axis=("T","Z","Y","X")
dimension_valid_max=(90000.0,12000.0,90.0,180.0)
dimension=valid_min=(-90000.0,-12000.0,-90.0,-180.0)
dimension_dict = {}
dimension_dict['flag'] = dict(zip(dimensions, dimension_flags))
dimension_dict['SDN'] = dict(zip(dimensions, SDN_dimension_names))
dimension_dict['unit'] = dict(zip(dimensions, dimension_units))
dimension_dict['long'] = dict(zip(dimensions, dimension_long_name))
dimension_dict['CF'] = dict(zip(dimensions, dimension_CF_standard_name))



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


