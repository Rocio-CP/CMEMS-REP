clc; clear all; close all

%%% Create CMEMS netcdf files
workdir=['/Users/rpr061/Dropbox/BCDC_Projects/CMEMS_INSTAC/REP_Products/'...
'2018/SOCATv6/workspace/'];
cd(workdir)
inputfile='SOCATv6_tracks_gridded_decades.nc'
outputfile='SOCATv6_tracks_gridded_decades_CMEMS_test2.nc'

% Read info
latin=ncread([workdir,inputfile],'ylat');
lonin=ncread([workdir,inputfile],'xlon');
timeintemp=ncread([workdir,inputfile],'tdecade');

tempin=ncread([workdir,inputfile],'sst_ave_weighted_decade');
salin=ncread([workdir,inputfile],'salinity_ave_weighted_decade');
fco2in=ncread([workdir,inputfile],'fco2_ave_weighted_decade');

% Transform time value from days from 1900 to days from 1950
btw0050 = daysact('01-jan-1900', '01-jan-1950');
timein = double(timeintemp - btw0050);

%% Create new file. Use lower level functions netcdf.*
ncid = netcdf.create([workdir,outputfile], 'NETCDF4');

%% Define dimensions and variables
% Define dimensions
latdim = netcdf.defDim(ncid, 'LATITUDE', 180);
londim = netcdf.defDim(ncid, 'LONGITUDE', 360);
posdim = netcdf.defDim(ncid, 'POSITION', 360);
timedim = netcdf.defDim(ncid, 'TIME', 5);
depthdim = netcdf.defDim(ncid, 'DEPTH', 1);

% Define dimension variables
latvar = netcdf.defVar(ncid, 'LATITUDE', 'float', latdim);
lonvar = netcdf.defVar(ncid, 'LONGITUDE', 'float', londim);
timevar = netcdf.defVar(ncid, 'TIME', 'double', timedim);
depthvar = netcdf.defVar(ncid, 'DEPTH', 'float', depthdim);
posvar = netcdf.defVar(ncid, 'POSITION', 'float', posdim);

% Define other variables
tempvar = netcdf.defVar(ncid, 'TEMP', 'float', [londim, latdim, timedim]);
salvar= netcdf.defVar(ncid, 'PSAL', 'float', [londim, latdim, timedim]);
fco2var= netcdf.defVar(ncid, 'FCO2', 'float', [londim, latdim, timedim]);

% Define _FillValue of variables (does not work after reDef);
netcdf.defVarFill(ncid, tempvar, false, -99999.0);
netcdf.defVarFill(ncid, salvar, false, -99999.0);
netcdf.defVarFill(ncid, fco2var, false, -99999.0);

%% Exit define mode
netcdf.endDef(ncid);

%% Put data
% Put data in dimensions
netcdf.putVar(ncid, latvar, latin);
netcdf.putVar(ncid, lonvar, lonin);
netcdf.putVar(ncid, timevar, timein);
netcdf.putVar(ncid, depthvar, 5.0);

% Put data in variables
netcdf.putVar(ncid, tempvar, tempin);
netcdf.putVar(ncid, salvar, salin);
netcdf.putVar(ncid, fco2var, fco2in);

%% Re-enter define mode (can be done before???)
netcdf.reDef(ncid);

%% Add attributes to dimensions
% Latitude
netcdf.putAtt(ncid, latvar, 'units', 'degrees_north');
netcdf.putAtt(ncid, latvar, 'axis', 'Y');
netcdf.putAtt(ncid, latvar, 'standard_name', 'latitude')
netcdf.putAtt(ncid, latvar, 'long_name', 'Latitude of each location');
netcdf.putAtt(ncid, latvar, 'valid_min', -90.0);
netcdf.putAtt(ncid, latvar, 'valid_max', 90.0);

% Longitude
netcdf.putAtt(ncid, lonvar, 'units', 'degrees_east');
netcdf.putAtt(ncid, lonvar, 'axis', 'X');
netcdf.putAtt(ncid, lonvar, 'standard_name', 'longitude')
netcdf.putAtt(ncid, lonvar, 'long_name', 'Longitude of each location');
netcdf.putAtt(ncid, lonvar, 'valid_min', -180.0);
netcdf.putAtt(ncid, lonvar, 'valid_max', 180.0);

% Time
netcdf.putAtt(ncid, timevar, 'units', 'days since 1950-01-01T00:00:00Z');
netcdf.putAtt(ncid, timevar, 'axis', 'T');
netcdf.putAtt(ncid, timevar, 'standard_name', 'time')
netcdf.putAtt(ncid, timevar, 'long_name', 'Time of measurement');
netcdf.putAtt(ncid, timevar, 'valid_min', 0.0);

% Depth
netcdf.putAtt(ncid, depthvar, 'units', 'meters');
netcdf.putAtt(ncid, depthvar, 'axis', 'Z');
netcdf.putAtt(ncid, depthvar, 'standard_name', 'depth');
netcdf.putAtt(ncid, depthvar, 'long_name', 'Depth of each measurement');
netcdf.putAtt(ncid, depthvar, 'valid_min', 0.0);
netcdf.putAtt(ncid, depthvar, 'valid_max', 12000.0);

% Position

%% Add attributes to variables
% Temperature
netcdf.putAtt(ncid, tempvar, 'units','degree_Celsius');
netcdf.putAtt(ncid, tempvar, 'standard_name','sea_water_temperature');
netcdf.putAtt(ncid, tempvar, 'long_name','sea_water_temperature');
netcdf.putAtt(ncid, tempvar, 'summary',...
    'The weighted cruise means for the months within each decade is averaged.')
netcdf.putAtt(ncid, tempvar, 'valid_min', -100.0);
netcdf.putAtt(ncid, tempvar, 'valid_max', 100.0);
netcdf.putAtt(ncid, tempvar, 'coordinates','longitude latitude time');

% Salinity
netcdf.putAtt(ncid, salvar, 'units','PSU');
netcdf.putAtt(ncid, salvar, 'standard_name','sea_water_salinity');
netcdf.putAtt(ncid, salvar, 'long_name','sea_water_salinity');
netcdf.putAtt(ncid, salvar, 'summary',...
    'The weighted cruise means for the months within each decade is averaged.')
netcdf.putAtt(ncid, salvar, 'valid_min', 0.0);
netcdf.putAtt(ncid, salvar, 'valid_max', 100.0);
netcdf.putAtt(ncid, tempvar, 'coordinates','longitude latitude time');

% fCO2
netcdf.putAtt(ncid, fco2var, 'units','microatmospheres');
netcdf.putAtt(ncid, fco2var, 'standard_name','fugacity_of_CO2');
netcdf.putAtt(ncid, fco2var, 'long_name','fugacity_of_CO2_in_sea_water');
netcdf.putAtt(ncid, fco2var, 'summary',...
    'The weighted cruise means for the months within each decade is averaged.')
netcdf.putAtt(ncid, fco2var, 'valid_min', 0.0);
netcdf.putAtt(ncid, fco2var, 'valid_max', 100000.0);
netcdf.putAtt(ncid, fco2var, 'coordinates','longitude latitude time');


%% Add global attributes
globid = netcdf.getConstant('GLOBAL');

netcdf.putAtt(ncid, globid, 'title', ...
    'SOCAT v6 Decadal 1x1 degree gridded dataset');
netcdf.putAtt(ncid, globid, 'data_mode', 'D');
netcdf.putAtt(ncid, globid, 'format_version', '1.2')
netcdf.putAtt(ncid, globid, 'data_type', 'OceanSITES trajectory data');
netcdf.putAtt(ncid, globid, 'area', 'Global');
netcdf.putAtt(ncid, globid, 'geospatial_lat_min', '-89.5');
netcdf.putAtt(ncid, globid, 'geospatial_lat_max', '89.5');
netcdf.putAtt(ncid, globid, 'geospatial_lat_units', 'degree_north');
netcdf.putAtt(ncid, globid, 'geospatial_lon_min', '-179.5');
netcdf.putAtt(ncid, globid, 'geospatial_lon_max', '179.5');
netcdf.putAtt(ncid, globid, 'geospatial_lon_units', 'degree_east');
netcdf.putAtt(ncid, globid, 'geospatial_vertical_min', '5.0');
netcdf.putAtt(ncid, globid, 'geospatial_vertical_max', '5.0');
netcdf.putAtt(ncid, globid, 'geospatial_vertical_positive', 'down');
netcdf.putAtt(ncid, globid, 'geospatial_vertical_units', 'meter');
netcdf.putAtt(ncid, globid, 'time_coverage_start', '1970-01-01');
netcdf.putAtt(ncid, globid, 'time_coverage_end', '2018-01-07');
netcdf.putAtt(ncid, globid, 'cdm_data_type', 'trajectory');
netcdf.putAtt(ncid, globid, 'keywords_vocabulary','SeaDataNet Parameter Discovery Vocabulary');
netcdf.putAtt(ncid, globid, 'data_assembly_center', 'BERGEN');
netcdf.putAtt(ncid, globid, 'netcdf_version', '4');
netcdf.putAtt(ncid, globid, 'update_interval', 'P1Y');
netcdf.putAtt(ncid, globid, 'format_version', '1.2');
netcdf.putAtt(ncid, globid, 'platform_code','NA');
netcdf.putAtt(ncid, globid, 'site_code','NA');
netcdf.putAtt(ncid, globid, 'date_update','2019-00-00T00:00:00Z');

netcdf.putAtt(ncid, globid, 'summary', ['Surface Ocean Carbon Atlas (SOCAT) Gridded (binned) SOCAT observations, with a spatial grid of ',...
          '1x1 degree and yearly in time. The gridded fields are computed from the monthly 1-degree gridded data, ',...
          'which uses only SOCAT datasets with QC flags of A through D and SOCAT data points flagged with WOCE ',...
          'flag values of 2. This decacal data is computed using data from the start to the end of each decade as ',...
          'described in the summary attribute of each variable. The first decade is 1-Jan-1970 through 31-Dec-1979, ',...
          'the last decade is generally a partial decade. No interpolation was performed. Significant biases are present ',...
          'in these gridded results due to the arbitrar and sparse locations of data values in both space and time']);
netcdf.putAtt(ncid, globid, 'references', ['Bakker et al. (2016). A multi-decade record of high-quality fCO2 data in version 3 of the ',...
          'Surface Ocean CO2 Atlas (SOCAT), Earth Syst. Sci. Data, 8, 383-413, https://doi.org/10.5194/essd-8-383-2016. ',...
          'Bakker et al. (2018): Surface Ocean CO2 Atlas (SOCAT) V6. PANGAEA, https://doi.org/10.1594/PANGAEA.890974'])
netcdf.putAtt(ncid, globid, 'contributor_name', ['Bakker, Dorothee C E; Lauvset, Siv K; Wanninkhof, Rik; Castaño-Primo, Rocío; ',...
          'Currie, Kim I; Jones, Stephen D; Landa, Camilla S; Metzl, Nicolas; Nakaoka, Shin-Ichiro; Nojiri, Yukihiro; ',...
          'Nonaka, Isao; O''Brien, Kevin M; Olsen, Are; Pfeil, Benjamin; Pierrot, Denis; Schuster, Ute; Smith, Karl; ',...
          'Sullivan, Kevin; Sutton, Adrienne; Tilbrook, Bronte; Alin, Simone; Becker, Meike; Benoit-Cattin, Alice; ',...
          'Bott, Randy; Bozec, Yann; Bozzano, Roberto; Burger, Eugene; Burgers, Tonya; Cai, Wei-Jun; Chen, Liqi; ',...
          'Chierici, Melissa; Corredor, Jorge; Cosca, Catherine E; Cross, Jessica; Dandonneau, Yves; De Carlo, Eric Heinen; ',...
          'Dietrich, Colin; Else, Brent; Emerson, Steven R; Farias, Laura; Fransson, Agneta; Garreaud, René D; Gkritzalis, Thanos; ',...
          'Glockzin, Michael; González-Dávila, Melchor; Gregor, Luke; Hartman, Sue E; Hermes, Rudolf; Hoppema, Mario; Howden, Stephan; ',...
          'Hunt, Christopher W; Hydes, David; Ibánhez, J Severino P; Kitidis, Vassilis; Körtzinger, Arne; Kozyr, Alexander; ',...
          'Kuwata, Akira; Lampitt, Richard Stephen; Lefèvre, Nathalie; Lo Monaco, Claire; Maenner, Stacy M; Manke, Ansley; ',...
          'Manzello, Derek P; McGillis, Wade; Mickett, John; Monteiro, Pedro M S; Morell, Julio; Morrison, Ru; Mucci, Alfonso; ',...
          'Munro, David R; Musielewicz, Sylvia; Negri, Ruben M; Newberger, Timothy; Newton, Jan; Noakes, Scott; O''Brien, Chris; ',...
          'Ólafsdóttir, Sólveig Rósa; Ólafsson, Jón; Ono, Tsuneo; Osborne, John; Ouyang, Zhangxian; Padín, Xose Antonio; Papakyriakou, Tim N; ',...
          'Plüddemann, Albert J; Rehder, Gregor; Sabine, Christopher L; Sakurai, Keizo; Salisbury, Joe; Santana-Casiano, Juana Magdalena; ',...
          'Schlitzer, Reiner; Schneider, Bernd; Send, Uwe; Skjelvan, Ingunn; Steinhoff, Tobias; Sulpis, Olivier; Sutherland, Stewart C; ',...
          'Sweeney, Colm; Tadokoro, Kazuaki; Takahashi, Taro; Telszewski, Maciej; Thomas, Helmuth; Tomlinson, Michael; Trull, Tom W; ',...
          'Valdimarsson, Hedinn; van Heuven, Steven; Vandemark, Doug; Wada, Chisato; Wallace, Douglas WR; Watson, Andrew J; Weller, Robert A; Xu, Suqing']);

%% Close nc file
netcdf.close(ncid);