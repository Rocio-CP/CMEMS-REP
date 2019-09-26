%% Create QC and DM variables for all variables
datavars={'tempvar','salvar','fco2var'};
allvars={'latvar','lonvar','timevar','depthvar','posvar','tempvar','salvar','fco2var'};

% only fco2 has woce flags, all =2 in SOCAT. fco2 and dimensions -> flag 1; temp, sal -> flag 0 (unknown)
%for ii=1:length(allvars);
%    currentvar=datavars{ii};
depthin=5.0;
depthinqc = 0;
depthindm = 'D';

timeinqc = zeros(size(timein));
possystin=repmat('U',size(lonin));
posinqc=zeros(size(lonin));


fco2inqc = repmat(0, size(fco2in));
fco2inqc(~isnan(fco2in))=1;
fco2indm = repmat('D',size(fco2in));

tempinqc = zeros(size(tempin));
tempindm = repmat('D',size(tempin));

salinqc = zeros(size(salin));
salindm = repmat('D',size(salin));


%end

%% Create new file. Use lower level functions netcdf.*
ncid = netcdf.create([workdir,outputfile], 'NETCDF4');

%% Define dimensions and variables
% Define dimensions
latdim = netcdf.defDim(ncid, 'LATITUDE', length(latin));
londim = netcdf.defDim(ncid, 'LONGITUDE', length(lonin));
posdim = netcdf.defDim(ncid, 'POSITION', length(lonin));
timedim = netcdf.defDim(ncid, 'TIME', netcdf.getConstant('NC_UNLIMITED'));
depthdim = netcdf.defDim(ncid, 'DEPH', 1);

% Define dimension variables
latvar = netcdf.defVar(ncid, 'LATITUDE', 'NC_FLOAT', latdim);
lonvar = netcdf.defVar(ncid, 'LONGITUDE', 'NC_FLOAT', londim);
timevar = netcdf.defVar(ncid, 'TIME', 'NC_DOUBLE', timedim);
timevarqc = netcdf.defVar(ncid, 'TIME_QC', 'NC_BYTE', timedim);

depthvar = netcdf.defVar(ncid, 'DEPH', 'NC_FLOAT', depthdim);
depthvarqc = netcdf.defVar(ncid, 'DEPH_QC', 'NC_BYTE', depthdim);
depthvardm = netcdf.defVar(ncid, 'DEPH_DM', 'NC_CHAR', depthdim);

posvarqc = netcdf.defVar(ncid, 'POSITION_QC', 'NC_BYTE', [posdim]);
possystvar = netcdf.defVar(ncid, 'POSITIONING_SYSTEM', 'NC_CHAR', [posdim]);


% Define other variables
tempvar = netcdf.defVar(ncid, 'TEMP', 'NC_DOUBLE', [londim, latdim, timedim]);
salvar= netcdf.defVar(ncid, 'PSAL', 'NC_DOUBLE', [londim, latdim, timedim]);
fco2var= netcdf.defVar(ncid, 'FCO2', 'NC_DOUBLE', [londim, latdim, timedim]);

tempvarqc = netcdf.defVar(ncid, 'TEMP_QC', 'NC_BYTE', [londim, latdim, timedim]);
salvarqc = netcdf.defVar(ncid, 'PSAL_QC', 'NC_BYTE', [londim, latdim, timedim]);
fco2varqc = netcdf.defVar(ncid, 'FCO2_QC', 'NC_BYTE', [londim, latdim, timedim]);

tempvardm = netcdf.defVar(ncid, 'TEMP_DM', 'NC_CHAR', [londim, latdim, timedim]);
salvardm = netcdf.defVar(ncid, 'PSAL_DM', 'NC_CHAR', [londim, latdim, timedim]);
fco2vardm = netcdf.defVar(ncid, 'FCO2_DM', 'NC_CHAR', [londim, latdim, timedim]);

% Define _FillValue of variables (does not work after reDef);
netcdf.defVarFill(ncid, timevar, false, 999999.0);
netcdf.defVarFill(ncid, latvar, false, 999999.0);
netcdf.defVarFill(ncid, lonvar, false, 999999.0);
netcdf.defVarFill(ncid, depthvar, false, -999999.0);
netcdf.defVarFill(ncid, possystvar, false, '');
netcdf.defVarFill(ncid, timevarqc, false, -128);
netcdf.defVarFill(ncid, posvarqc, false, -128);
netcdf.defVarFill(ncid, depthvarqc, false, -128);
netcdf.defVarFill(ncid, depthvardm, false, '');

netcdf.defVarFill(ncid, tempvar, false, -9999);
netcdf.defVarFill(ncid, salvar, false, -9999);
netcdf.defVarFill(ncid, fco2var, false, -9999);

netcdf.defVarFill(ncid, tempvarqc, false, -128);
netcdf.defVarFill(ncid, salvarqc, false, -128);
netcdf.defVarFill(ncid, fco2varqc, false, -128);

netcdf.defVarFill(ncid, tempvardm, false, '');
netcdf.defVarFill(ncid, salvardm, false, '');
netcdf.defVarFill(ncid, fco2vardm, false, '');

%% Exit define mode
netcdf.endDef(ncid);

%% Put data
% Put data in dimensions
netcdf.putVar(ncid, latvar, latin);
netcdf.putVar(ncid, lonvar, lonin);
netcdf.putVar(ncid, timevar, 0, length(timein), timein);
netcdf.putVar(ncid, timevarqc, timeinqc);
netcdf.putVar(ncid, depthvar, depthin);
netcdf.putVar(ncid, depthvarqc, depthinqc);
netcdf.putVar(ncid, depthvardm, depthindm);
netcdf.putVar(ncid, possystvar, possystin);
netcdf.putVar(ncid, posvarqc, posinqc);

% Put data in variables
netcdf.putVar(ncid, tempvar, tempin);
netcdf.putVar(ncid, salvar, salin);
netcdf.putVar(ncid, fco2var, fco2in);

netcdf.putVar(ncid, tempvarqc, tempinqc);
netcdf.putVar(ncid, salvarqc, salinqc);
netcdf.putVar(ncid, fco2varqc, fco2inqc);

netcdf.putVar(ncid, tempvardm, tempindm);
netcdf.putVar(ncid, salvardm, salindm);
netcdf.putVar(ncid, fco2vardm, fco2indm);

%% Re-enter define mode (can be done before???)
netcdf.reDef(ncid);

%% Add attributes to dimensions
% Add attributes to dimensions
netcdf.putAtt(ncid, latvar, 'units', 'degree_north');
netcdf.putAtt(ncid, latvar, 'reference', 'WGS84');
netcdf.putAtt(ncid, latvar, 'axis', 'Y');
netcdf.putAtt(ncid, latvar, 'standard_name', 'latitude')
netcdf.putAtt(ncid, latvar, 'long_name', 'Latitude of each location');
netcdf.putAtt(ncid, latvar, 'valid_min', -90.0);
netcdf.putAtt(ncid, latvar, 'valid_max', 90.0);

% Longitude
netcdf.putAtt(ncid, lonvar, 'units', 'degree_east');
netcdf.putAtt(ncid, lonvar, 'reference', 'WGS84');
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
netcdf.putAtt(ncid, timevar, 'valid_min', -90000);
netcdf.putAtt(ncid, timevar, 'valid_max', 90000);
netcdf.putAtt(ncid, timevar, 'comment', '');
% Time QC
netcdf.putAtt(ncid, timevarqc, 'long_name', 'quality flag');
netcdf.putAtt(ncid, timevarqc, 'conventions', 'OceanSITES reference table 2');
netcdf.putAtt(ncid, timevarqc, 'valid_min', 0);
netcdf.putAtt(ncid, timevarqc, 'valid_max', 9);
netcdf.putAtt(ncid, timevarqc, 'flag_values', [0  1  2  3  4  5  6  7  8  9]);
netcdf.putAtt(ncid, timevarqc, 'flag_meanings','no_qc_performed good_data probably_good_data bad_data_that_are_potentially_correctable bad_data value_changed not_used nominal_value interpolated_value missing_value');

% Depth
netcdf.putAtt(ncid, depthvar, 'units', 'm');
netcdf.putAtt(ncid, depthvar, 'axis', 'Z');
netcdf.putAtt(ncid, depthvar, 'standard_name', 'depth');
netcdf.putAtt(ncid, depthvar, 'long_name', 'Depth');
netcdf.putAtt(ncid, depthvar, 'valid_min', -12000.0);
netcdf.putAtt(ncid, depthvar, 'valid_max', 12000.0);
netcdf.putAtt(ncid, depthvar, 'positive', 'down');
netcdf.putAtt(ncid, depthvar, 'reference', 'sea_level');
% Depth QC
netcdf.putAtt(ncid, depthvarqc, 'long_name', 'quality flag');
netcdf.putAtt(ncid, depthvarqc, 'conventions', 'OceanSITES reference table 2');
netcdf.putAtt(ncid, depthvarqc, 'valid_min', 0);
netcdf.putAtt(ncid, depthvarqc, 'valid_max', 9);
netcdf.putAtt(ncid, depthvarqc, 'flag_values', [0  1  2  3  4  5  6  7  8  9]);
netcdf.putAtt(ncid, depthvarqc, 'flag_meanings','no_qc_performed good_data probably_good_data bad_data_that_are_potentially_correctable bad_data value_changed not_used nominal_value interpolated_value missing_value');
% Depth DM
netcdf.putAtt(ncid, depthvardm, 'long_name', 'method of data processing');
netcdf.putAtt(ncid, depthvardm, 'conventions', 'OceanSITES reference table 5');
netcdf.putAtt(ncid, depthvardm, 'flag_values', 'R, P, D, M');
netcdf.putAtt(ncid, depthvardm, 'flag_meanings','real-time provisional delayed-mode mixed');

% Positioning system
netcdf.putAtt(ncid, possystvar, 'long_name', 'Positioning system');
netcdf.putAtt(ncid, possystvar, 'flag_values', 'A, G, L, N, U');
netcdf.putAtt(ncid, possystvar, 'flag_meanings', 'Argos, GPS, Loran, Nominal, Unknown');

netcdf.putAtt(ncid, posvarqc, 'long_name', 'quality flag');
netcdf.putAtt(ncid, posvarqc, 'conventions', 'OceanSITES reference table 2');
netcdf.putAtt(ncid, posvarqc, 'valid_min', 0);
netcdf.putAtt(ncid, posvarqc, 'valid_max', 9);
netcdf.putAtt(ncid, posvarqc, 'flag_values', [0  1  2  3  4  5  6  7  8  9]);
netcdf.putAtt(ncid, posvarqc, 'flag_meanings','no_qc_performed good_data probably_good_data bad_data_that_are_potentially_correctable bad_data value_changed not_used nominal_value interpolated_value missing_value');

%% Add attributes to variables
% Temperature
netcdf.putAtt(ncid, tempvar, 'units','degrees_C');
netcdf.putAtt(ncid, tempvar, 'standard_name','sea_water_temperature');
netcdf.putAtt(ncid, tempvar, 'long_name','Sea temperature');
netcdf.putAtt(ncid, tempvar, 'summary', attsummary)
%netcdf.putAtt(ncid, tempvar, 'valid_min', -100.0);
%netcdf.putAtt(ncid, tempvar, 'valid_max', 100.0);

% Salinity
netcdf.putAtt(ncid, salvar, 'units','0.001');
netcdf.putAtt(ncid, salvar, 'standard_name','sea_water_salinity');
netcdf.putAtt(ncid, salvar, 'long_name','Practical salinity');
netcdf.putAtt(ncid, salvar, 'summary', attsummary)
%netcdf.putAtt(ncid, salvar, 'valid_min', 0.0);
%netcdf.putAtt(ncid, salvar, 'valid_max', 100.0);

% fCO2
netcdf.putAtt(ncid, fco2var, 'units','microatmosphere');
netcdf.putAtt(ncid, fco2var, 'standard_name','fugacity_of_CO2');
netcdf.putAtt(ncid, fco2var, 'long_name','fugacity_of_CO2_in_sea_water');
netcdf.putAtt(ncid, fco2var, 'summary', attsummary)
%netcdf.putAtt(ncid, fco2var, 'valid_min', 0.0);
%netcdf.putAtt(ncid, fco2var, 'valid_max', 100000.0);

% QC variables
netcdf.putAtt(ncid, tempvarqc, 'long_name', 'quality flag');
netcdf.putAtt(ncid, tempvarqc, 'conventions', 'OceanSITES reference table 2');
netcdf.putAtt(ncid, tempvarqc, 'valid_min', 0);
netcdf.putAtt(ncid, tempvarqc, 'valid_max', 9);
netcdf.putAtt(ncid, tempvarqc, 'flag_values', [0  1  2  3  4  5  6  7  8  9]);
netcdf.putAtt(ncid, tempvarqc, 'flag_meanings','no_qc_performed good_data probably_good_data bad_data_that_are_potentially_correctable bad_data value_changed not_used nominal_value interpolated_value missing_value');

netcdf.putAtt(ncid, salvarqc, 'long_name', 'quality flag');
netcdf.putAtt(ncid, salvarqc, 'conventions', 'OceanSITES reference table 2');
netcdf.putAtt(ncid, salvarqc, 'valid_min', 0);
netcdf.putAtt(ncid, salvarqc, 'valid_max', 9);
netcdf.putAtt(ncid, salvarqc, 'flag_values', [0  1  2  3  4  5  6  7  8  9]);
netcdf.putAtt(ncid, salvarqc, 'flag_meanings','no_qc_performed good_data probably_good_data bad_data_that_are_potentially_correctable bad_data value_changed not_used nominal_value interpolated_value missing_value');

netcdf.putAtt(ncid, fco2varqc, 'long_name', 'quality flag');
netcdf.putAtt(ncid, fco2varqc, 'conventions', 'OceanSITES reference table 2');
netcdf.putAtt(ncid, fco2varqc, 'valid_min', 0);
netcdf.putAtt(ncid, fco2varqc, 'valid_max', 9);
netcdf.putAtt(ncid, fco2varqc, 'flag_values', [0  1  2  3  4  5  6  7  8  9]);
netcdf.putAtt(ncid, fco2varqc, 'flag_meanings','no_qc_performed good_data probably_good_data bad_data_that_are_potentially_correctable bad_data value_changed not_used nominal_value interpolated_value missing_value');

% DM variables
netcdf.putAtt(ncid, tempvardm, 'long_name', 'method of data processing');
netcdf.putAtt(ncid, tempvardm, 'conventions', 'OceanSITES reference table 5');
netcdf.putAtt(ncid, tempvardm, 'flag_values', 'R, P, D, M');
netcdf.putAtt(ncid, tempvardm, 'flag_meanings','real-time provisional delayed-mode mixed');

netcdf.putAtt(ncid, salvardm, 'long_name', 'method of data processing');
netcdf.putAtt(ncid, salvardm, 'conventions', 'OceanSITES reference table 5');
netcdf.putAtt(ncid, salvardm, 'flag_values', 'R, P, D, M');
netcdf.putAtt(ncid, salvardm, 'flag_meanings','real-time provisional delayed-mode mixed');

netcdf.putAtt(ncid, fco2vardm, 'long_name', 'method of data processing');
netcdf.putAtt(ncid, fco2vardm, 'conventions', 'OceanSITES reference table 5');
netcdf.putAtt(ncid, fco2vardm, 'flag_values', 'R, P, D, M');
netcdf.putAtt(ncid, fco2vardm, 'flag_meanings','real-time provisional delayed-mode mixed');



%% Add global attributes
globid = netcdf.getConstant('GLOBAL');

%netcdf.putAtt(ncid, globid, 'data_type', 'OceanSITES trajectory data');
netcdf.putAtt(ncid, globid, 'Conventions', 'CF-1.6 OceanSITES-Manual-1.2 Copernicus-InSituTAC-SRD-1.3 Copernicus-InSituTAC-ParametersList-3.1.0');
netcdf.putAtt(ncid, globid, 'format_version', '1.2');

netcdf.putAtt(ncid, globid, 'title', filetitle);
netcdf.putAtt(ncid, globid, 'id', outputfile(1:end-3));
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
netcdf.putAtt(ncid, globid, 'author', 'Bjerknes Climate Data Centre');

netcdf.putAtt(ncid, globid, 'netcdf_version', '4');
netcdf.putAtt(ncid, globid, 'update_interval', 'yearly');
netcdf.putAtt(ncid, globid, 'date_update', datestr(today, 'yyyy-mm-ddTHH:MM:SSZ'));
netcdf.putAtt(ncid, globid, 'history', [datestr(today, 'yyyy-mm-ddTHH:MM:SSZ'),' : Creation']);

netcdf.putAtt(ncid, globid, 'qc_manual', 'SOCAT Quality Control Cookbook v3: https://www.socat.info/wp-content/uploads/2017/04/2015_SOCAT_QC_Cookbook_v3.pdf');
netcdf.putAtt(ncid, globid, 'quality_control_indicator', '6');
netcdf.putAtt(ncid, globid, 'quality_index', 'A');

netcdf.putAtt(ncid, globid, 'summary', globalsummary);
netcdf.putAtt(ncid, globid, 'references', 'http://marine.copernicus.eu https://www.socat.info/');
netcdf.putAtt(ncid, globid, 'comment', ['Further details can be found in Bakker et al. (2016). A multi-decade record of high-quality fCO2 data in version 3 of the ',...
    'Surface Ocean CO2 Atlas (SOCAT), Earth Syst. Sci. Data, 8, 383-413, https://doi.org/10.5194/essd-8-383-2016. ',...
    'and Bakker et al. (2018): Surface Ocean CO2 Atlas (SOCAT) V6. PANGAEA, https://doi.org/10.1594/PANGAEA.890974'])
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