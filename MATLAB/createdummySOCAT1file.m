load SOCATv6.mat
%% DUMMY 1-file SOCAT

latin = latitudedecdegN;
lonin = longitudedecdegE;
timein = datenum(double([yr,mon,day1,hh,mm,ss]));
depthin = 5.0;
posin = zeros(size(latin));

fco2in = fCO2recuatm;
fco2qcin = fCO2rec_flag;
fco2dmin = repmat('d',size(fco2in));

tempin = SSTdegC;
tempqcin = repmat(0, size(tempin));
tempdmin = repmat('d',size(tempin));

salin = sal;
salqcin = repmat(0, size(salin));
saldmin = repmat('d',size(salin));

%
ncid = netcdf.create([workdir,outputfile], 'NETCDF4');

% Define dimensions and variables
% Define dimensions
latdim = netcdf.defDim(ncid, 'LATITUDE', length(latin));
londim = netcdf.defDim(ncid, 'LONGITUDE', length(lonin));
posdim = netcdf.defDim(ncid, 'POSITION', length(lonin));
timedim = netcdf.defDim(ncid, 'TIME', length(timein));
depthdim = netcdf.defDim(ncid, 'DEPTH', 1);

% Define dimension variables
latvar = netcdf.defVar(ncid, 'LATITUDE', 'NC_DOUBLE', latdim);
lonvar = netcdf.defVar(ncid, 'LONGITUDE', 'NC_DOUBLE', londim);
timevar = netcdf.defVar(ncid, 'TIME', 'NC_DOUBLE', timedim);
depthvar = netcdf.defVar(ncid, 'DEPTH', 'NC_DOUBLE', depthdim);
posvar = netcdf.defVar(ncid, 'POSITION', 'NC_DOUBLE', posdim);

% Define other variables
tempvar = netcdf.defVar(ncid, 'TEMP', 'NC_DOUBLE', [depthdim, timedim]);
salvar= netcdf.defVar(ncid, 'PSAL', 'NC_DOUBLE', [depthdim, timedim]);
fco2var= netcdf.defVar(ncid, 'FCO2', 'NC_DOUBLE', [depthdim, timedim]);

tempvarqc = netcdf.defVar(ncid, 'TEMP_QC', 'NC_BYTE', [depthdim, timedim]);
salvarqc = netcdf.defVar(ncid, 'PSAL_QC', 'NC_BYTE', [depthdim, timedim]);
fco2varqc = netcdf.defVar(ncid, 'FCO2_QC', 'NC_BYTE', [depthdim, timedim]);

tempvardm = netcdf.defVar(ncid, 'TEMP_DM', 'NC_CHAR', [depthdim, timedim]);
salvardm = netcdf.defVar(ncid, 'PSAL_DM', 'NC_CHAR', [depthdim, timedim]);
fco2vardm = netcdf.defVar(ncid, 'FCO2_DM', 'NC_CHAR', [depthdim, timedim]);

% Define _FillValue of variables (does not work after reDef);
netcdf.defVarFill(ncid, tempvar, false, -99999.0);
netcdf.defVarFill(ncid, salvar, false, -99999.0);
netcdf.defVarFill(ncid, fco2var, false, -99999.0);

% Exit define mode
netcdf.endDef(ncid);

% Put data
% Put data in dimensions
netcdf.putVar(ncid, latvar, latin);
netcdf.putVar(ncid, lonvar, lonin);
netcdf.putVar(ncid, timevar, timein);
netcdf.putVar(ncid, depthvar, 5.0);

% Put data in variables
netcdf.putVar(ncid, tempvar, tempin);
netcdf.putVar(ncid, salvar, salin);
netcdf.putVar(ncid, fco2var, fco2in);

netcdf.putVar(ncid, tempvarqc, tempqcin);
netcdf.putVar(ncid, salvarqc, salqcin);
netcdf.putVar(ncid, fco2varqc, fco2qcin);

netcdf.putVar(ncid, tempvardm, tempdmin);
netcdf.putVar(ncid, salvardm, saldmin);
netcdf.putVar(ncid, fco2vardm, fco2dmin);
% Re-enter define mode (can be done before???)
netcdf.reDef(ncid);

% Add attributes to dimensions
% Latitude
netcdf.putAtt(ncid, latvar, 'units', 'degree_north');
netcdf.putAtt(ncid, latvar, 'axis', 'Y');
netcdf.putAtt(ncid, latvar, 'standard_name', 'latitude')
netcdf.putAtt(ncid, latvar, 'long_name', 'Latitude of each location');
netcdf.putAtt(ncid, latvar, 'valid_min', -90.0);
netcdf.putAtt(ncid, latvar, 'valid_max', 90.0);

% Longitude
netcdf.putAtt(ncid, lonvar, 'units', 'degree_east');
netcdf.putAtt(ncid, lonvar, 'axis', 'X');
netcdf.putAtt(ncid, lonvar, 'standard_name', 'longitude')
netcdf.putAtt(ncid, lonvar, 'long_name', 'Longitude of each location');
netcdf.putAtt(ncid, lonvar, 'valid_min', -180.0);
netcdf.putAtt(ncid, lonvar, 'valid_max', 180.0);

% Time
netcdf.putAtt(ncid, timevar, 'units', 'days since 1950-01-01T00:00:00Z');
netcdf.putAtt(ncid, timevar, 'axis', 'T');
netcdf.putAtt(ncid, timevar, 'standard_name', 'time')
netcdf.putAtt(ncid, timevar, 'long_name', 'Time');
netcdf.putAtt(ncid, timevar, 'valid_min', 0.0);

% Depth
netcdf.putAtt(ncid, depthvar, 'units', 'm');
netcdf.putAtt(ncid, depthvar, 'axis', 'Z');
netcdf.putAtt(ncid, depthvar, 'standard_name', 'depth');
netcdf.putAtt(ncid, depthvar, 'long_name', 'Depth');
netcdf.putAtt(ncid, depthvar, 'valid_min', 0.0);
netcdf.putAtt(ncid, depthvar, 'valid_max', 12000.0);

% Position

% Add attributes to variables
% Temperature
netcdf.putAtt(ncid, tempvar, 'units','degrees_C');
netcdf.putAtt(ncid, tempvar, 'standard_name','sea_water_temperature');
netcdf.putAtt(ncid, tempvar, 'long_name','sea_water_temperature');
netcdf.putAtt(ncid, tempvar, 'summary', '')
netcdf.putAtt(ncid, tempvar, 'valid_min', -100.0);
netcdf.putAtt(ncid, tempvar, 'valid_max', 100.0);
netcdf.putAtt(ncid, tempvar, 'coordinates','longitude latitude time');

% Salinity
netcdf.putAtt(ncid, salvar, 'units','PSU');
netcdf.putAtt(ncid, salvar, 'standard_name','sea_water_salinity');
netcdf.putAtt(ncid, salvar, 'long_name','sea_water_salinity');
netcdf.putAtt(ncid, salvar, 'summary', '')
netcdf.putAtt(ncid, salvar, 'valid_min', 0.0);
netcdf.putAtt(ncid, salvar, 'valid_max', 100.0);
netcdf.putAtt(ncid, tempvar, 'coordinates','longitude latitude time');

% fCO2
netcdf.putAtt(ncid, fco2var, 'units','microatmospheres');
netcdf.putAtt(ncid, fco2var, 'standard_name','fugacity_of_CO2');
netcdf.putAtt(ncid, fco2var, 'long_name','fugacity_of_CO2_in_sea_water');
netcdf.putAtt(ncid, fco2var, 'summary', '')
netcdf.putAtt(ncid, fco2var, 'valid_min', 0.0);
netcdf.putAtt(ncid, fco2var, 'valid_max', 100000.0);
netcdf.putAtt(ncid, fco2var, 'coordinates','longitude latitude time');

% QC variables


% DM variables


% Add global attributes
globid = netcdf.getConstant('GLOBAL');

netcdf.putAtt(ncid, globid, 'data_type', 'OceanSITES trajectory data');
netcdf.putAtt(ncid, globid, 'Conventions', 'CF-1.6 OceanSITES-Manual-1.2 Copernicus-InSituTAC-SRD-1.3 Copernicus-InSituTAC-ParametersList-3.1.0');
netcdf.putAtt(ncid, globid, 'format_version', '1.2');

netcdf.putAtt(ncid, globid, 'title', 'SOCATv6dummy');
netcdf.putAtt(ncid, globid, 'id', 'NAMEOFFILEWITHOUTEXTENSION');
netcdf.putAtt(ncid, globid, 'data_mode', 'D');
netcdf.putAtt(ncid, globid, 'format_version', '1.2')
netcdf.putAtt(ncid, globid, 'area', 'Global_Ocean');

netcdf.putAtt(ncid, globid, 'source', 'NA');
netcdf.putAtt(ncid, globid, 'source_platform_category_code', 'NA');
netcdf.putAtt(ncid, globid, 'platform_code','NA');
netcdf.putAtt(ncid, globid, 'platform_name','NA');
netcdf.putAtt(ncid, globid, 'site_code','NA');

netcdf.putAtt(ncid, globid, 'citation', 'These data were collected and made freely available by the Copernicus project and the programs that contribute to it.');
netcdf.putAtt(ncid, globid, 'distribution_statement','These data follow Copernicus standards; they are public and free of charge. User assumes all risk for use of data. User must display citation in any publication or product using data. User must contact PI prior to any commercial use of data.');
netcdf.putAtt(ncid, globid, 'naming_authority', 'Copernicus');

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

netcdf.putAtt(ncid, globid, 'cdm_data_type', 'grid');
netcdf.putAtt(ncid, globid, 'keywords_vocabulary','SeaDataNet Parameter Discovery Vocabulary');

netcdf.putAtt(ncid, globid, 'data_assembly_center', 'BERGEN');
netcdf.putAtt(ncid, globid, 'institution', 'University of Bergen Geophysical Institute');
netcdf.putAtt(ncid, globid, 'institution_edmo_code','4595');
netcdf.putAtt(ncid, globid, 'contact', 'post@bcdc.uib.no');
netcdf.putAtt(ncid, globid, 'author', 'cmems-service');

netcdf.putAtt(ncid, globid, 'netcdf_version', '4');
netcdf.putAtt(ncid, globid, 'update_interval', 'yearly');
netcdf.putAtt(ncid, globid, 'date_update', datestr(today, 'yyyy-mm-ddTHH:MM:SSZ'));
netcdf.putAtt(ncid, globid, 'history', [datestr(today, 'yyyy-mm-ddTHH:MM:SSZ'),' : Creation']);

netcdf.putAtt(ncid, globid, 'qc_manual', 'SOCAT Quality Control Cookbook v3: https://www.socat.info/wp-content/uploads/2017/04/2015_SOCAT_QC_Cookbook_v3.pdf');
netcdf.putAtt(ncid, globid, 'quality_control_indicator', '6');
netcdf.putAtt(ncid, globid, 'quality_index', 'A');

netcdf.putAtt(ncid, globid, 'summary', '');
netcdf.putAtt(ncid, globid, 'references', 'http://marine.copernicus.eu https://www.socat.info/');
netcdf.putAtt(ncid, globid, 'comment', ['Further details can be found in Bakker et al. (2016). A multi-decade record of high-quality fCO2 data in version 3 of the ',...
          'Surface Ocean CO2 Atlas (SOCAT), Earth Syst. Sci. Data, 8, 383-413, https://doi.org/10.5194/essd-8-383-2016. ',...
          'and Bakker et al. (2018): Surface Ocean CO2 Atlas (SOCAT) V6. PANGAEA, https://doi.org/10.1594/PANGAEA.890974']);
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

% Close nc file
netcdf.close(ncid);
