%clc; clear all; close all

%%% Create CMEMS netcdf files
workdir=['/Users/rpr061/Dropbox/BCDC_Projects/CMEMS_INSTAC/REP_Products/'...
    '2018/GLODAPv2/GLODAPv2.2016b_MappedClimatologies/'];
cd(workdir)
outputfile=['GLODAPv2_gridded_Copernicus_test.nc'];

%% CHECK TITLES IN THE PIT DOCUMENT!! OR LOICs MAIL. OR PUM?
inputfile='GLODAPv2.2016b.temperature.nc';
% Read info
% Dimension variables are the same in all climatologies
latin=ncread([workdir,inputfile],'lat');
lonin=ncread([workdir,inputfile],'lon');
depthin=ncread([workdir,inputfile],'Depth');
depthinqc=zeros(size(depthin));
depthindm=repmat('D',size(depthin));
timein=0;
timeinqc=0;
posinqc=zeros(size(lonin));
possystin=repmat('U',size(posinqc));


tempin=ncread('./GLODAPv2.2016b.temperature.nc','temperature');
salin=ncread('./GLODAPv2.2016b.salinity.nc','salinity');
oxyin=ncread('./GLODAPv2.2016b.oxygen.nc','oxygen');
nitrain=ncread('./GLODAPv2.2016b.NO3.nc','NO3');
phosin=ncread('./GLODAPv2.2016b.PO4.nc', 'PO4');
siliin=ncread('./GLODAPv2.2016b.silicate.nc', 'silicate');
phin=ncread('./GLODAPv2.2016b.pHtsinsitutp.nc','pHtsinsitutp');
ph25in=ncread('./GLODAPv2.2016b.pHts25p0.nc', 'pHts25p0');
dicin=ncread('./GLODAPv2.2016b.TCO2.nc', 'TCO2');
alkin=ncread('./GLODAPv2.2016b.TAlk.nc', 'TAlk');

vars={'temp','sal','oxy','nitra','phos','sili','ph','ph25','dic','alk'};
osvars={'TEMP','PSAL','DOX2','NTAW','PHOW','SLCW','PHPH','PH25','TICW','ALKW'};
osunits={'degrees_C','0.001','micromol_per_kg',...
    'micromol_per_kg','micromol_per_kg',...
    'micromol_per_kg','','',...
    'micromol_per_kg','micromol_per_kg',};
longname={'Sea temperature','Practical salinity','Dissolved Oxygen',...
    'Nitrate (NO3-N)','Phosphate (PO-P)',...
    'Silicate (SiO4-Si)','Ph','Ph at 25 degrees and 0 dbar',...
    'Dissolved inorganic carbon','Total alkalinity'};
standardname={'sea_water_temperature','sea_water_practical_salinity','moles_of_oxygen_per_unit_mass_in_sea_water',...
    'moles_of_nitrate_per_unit_mass_in_sea_water','moles_of_phosphate_per_unit_mass_in_sea_water',...
    'moles_of_silicate_per_unit_mass_in_sea_water','sea_water_ph_reported_on_total_scale','sea_water_ph_reported_on_total_scale_at_25_degrees_and_0_dbar',...
    'moles_of_dissolved_inorganic_carbon_per_unit_mass_in_sea_water','sea_water_alkalinity_per_unit_mass'};

% QC and DM
for vv=1:length(vars);
     eval([vars{vv},'inqc=repmat(1,size(',vars{vv},'in));']);
    % DM variables
    eval([vars{vv},'indm=repmat(''D'',size(',vars{vv},'in));']);
end

%% 
%% Generate NetCDF
% Create NetCDF file
ncid = netcdf.create([workdir,outputfile], 'NETCDF4');

% Define dimensions and variables
% Define dimensions
latdim = netcdf.defDim(ncid, 'LATITUDE', length(latin));
londim = netcdf.defDim(ncid, 'LONGITUDE', length(lonin));
posdim = netcdf.defDim(ncid, 'POSITION', length(lonin));
timedim = netcdf.defDim(ncid, 'TIME', netcdf.getConstant('NC_UNLIMITED'));
depthdim = netcdf.defDim(ncid, 'DEPH', length(depthin));

% Define dimension-ish variables
latvar = netcdf.defVar(ncid, 'LATITUDE', 'NC_FLOAT', latdim);
lonvar = netcdf.defVar(ncid, 'LONGITUDE', 'NC_FLOAT', londim);
timevar = netcdf.defVar(ncid, 'TIME', 'NC_FLOAT', timedim);
depthvar = netcdf.defVar(ncid, 'DEPH', 'NC_FLOAT', depthdim);
possystvar = netcdf.defVar(ncid, 'POSITIONING_SYSTEM', 'NC_CHAR', posdim);

%posvar = netcdf.defVar(ncid, 'POSITION', 'NC_FLOAT', [timedim]);
% Dimension variables QC DM
depthvarqc = netcdf.defVar(ncid, 'DEPH_QC', 'NC_BYTE', depthdim);
depthvardm = netcdf.defVar(ncid, 'DEPH_DM', 'NC_CHAR', depthdim);
posvarqc = netcdf.defVar(ncid, 'POSITION_QC', 'NC_BYTE', posdim);
timevarqc = netcdf.defVar(ncid, 'TIME_QC', 'NC_BYTE', timedim);


% Define other variables
for vv=1:length(vars);
eval([vars{vv},'var=netcdf.defVar(ncid,''',osvars{vv},''',''NC_FLOAT'',[londim latdim depthdim]);']);
eval([vars{vv},'varqc=netcdf.defVar(ncid,''',osvars{vv},'_QC'',''NC_BYTE'',[londim latdim depthdim]);']);
eval([vars{vv},'vardm=netcdf.defVar(ncid,''',osvars{vv},'_DM'',''NC_CHAR'',[londim latdim depthdim]);']);
end

netcdf.defVarFill(ncid, timevar, false, 999999.0);
netcdf.defVarFill(ncid, latvar, false, 999999.0);
netcdf.defVarFill(ncid, lonvar, false, 999999.0);
netcdf.defVarFill(ncid, depthvar, false, -999999.0);
netcdf.defVarFill(ncid, possystvar, false, '');
netcdf.defVarFill(ncid, timevarqc, false, -128);
netcdf.defVarFill(ncid, posvarqc, false, -128);
netcdf.defVarFill(ncid, depthvarqc, false, -128);
netcdf.defVarFill(ncid, depthvardm, false, '');

% Define _FillValue of variables (does not work after reDef);
for vv=1:length(vars);
eval(['netcdf.defVarFill(ncid,',vars{vv},'var,false,-9999);']);
eval(['netcdf.defVarFill(ncid,',vars{vv},'varqc,false,-128);']);
eval(['netcdf.defVarFill(ncid,',vars{vv},'vardm,false,'''');']);
end

% Exit define mode
netcdf.endDef(ncid);

% Put data
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
for vv=1:length(vars);
eval(['netcdf.putVar(ncid,',vars{vv},'var,',vars{vv},'in);']);
eval(['netcdf.putVar(ncid,',vars{vv},'varqc,',vars{vv},'inqc);']);
eval(['netcdf.putVar(ncid,',vars{vv},'vardm,',vars{vv},'indm);']);
end


% Re-enter define mode (can be done before???)
netcdf.reDef(ncid);

% Add attributes to dimensions
% Latitude
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

%netcdf.putAtt(ncid, posvar, 'units', 'aa');
%netcdf.putAtt(ncid, posvar, 'standard_name', 'Position of each measurement');

% Add attributes to variables
for vv=1:length(vars);
% Variables
eval(['netcdf.putAtt(ncid,',vars{vv},'var,''units'',''',osunits{vv},''');']);
eval(['netcdf.putAtt(ncid,',vars{vv},'var,''standard_name'',''',standardname{vv},''');']);
eval(['netcdf.putAtt(ncid,',vars{vv},'var,''long_name'',''',longname{vv},''');']);
%if vv==2;   
%eval(['netcdf.putAtt(ncid,',vars{vv},'var,''valid_min'',-100.0);']);
%else
%eval(['netcdf.putAtt(ncid,',vars{vv},'var,''valid_min'',0.0);']);
%end
%eval(['netcdf.putAtt(ncid,',vars{vv},'var,''valid_max'',100000.0);']);
%eval(['netcdf.putAtt(ncid,',vars{vv},'var,''coordinates'',''depth time'');']);

% QC
eval(['netcdf.putAtt(ncid,',vars{vv},'varqc,''long_name'',''quality flag'');']);
eval(['netcdf.putAtt(ncid,',vars{vv},'varqc,''conventions'',''OceanSITES reference table 2'');']);
eval(['netcdf.putAtt(ncid,',vars{vv},'varqc,''valid_min'',0);']);
eval(['netcdf.putAtt(ncid,',vars{vv},'varqc,''valid_max'',9);']);
eval(['netcdf.putAtt(ncid,',vars{vv},'varqc,''flag_values'',[0  1  2  3  4  5  6  7  8  9]);']);
eval(['netcdf.putAtt(ncid,',vars{vv},'varqc,''flag_meanings'',''no_qc_performed good_data probably_good_data bad_data_that_are_potentially_correctable bad_data value_changed not_used nominal_value interpolated_value missing_value'');']);

% DM
eval(['netcdf.putAtt(ncid,',vars{vv},'vardm,''long_name'',''method of data processing'');']);
eval(['netcdf.putAtt(ncid,',vars{vv},'vardm,''conventions'',''OceanSITES reference table 5'');']);
eval(['netcdf.putAtt(ncid,',vars{vv},'vardm,''flag_values'',''R, P, D, M'');']);
eval(['netcdf.putAtt(ncid,',vars{vv},'vardm,''flag_meanings'',''real-time provisional delayed-mode mixed'');']);

end

% Add global attributes
globid = netcdf.getConstant('GLOBAL');

%netcdf.putAtt(ncid, globid, 'data_type', 'OceanSITES profile data');
netcdf.putAtt(ncid, globid, 'Conventions', 'CF-1.6 OceanSITES-Manual-1.2 Copernicus-InSituTAC-SRD-1.3 Copernicus-InSituTAC-ParametersList-3.1.0');
netcdf.putAtt(ncid, globid, 'format_version', '1.2');

netcdf.putAtt(ncid, globid, 'title', 'Global Ocean - Global Ocean Data Analysis Project (GLODAP), 1x1 degree Climatology');
netcdf.putAtt(ncid, globid, 'id', outputfile(1:end-3));
netcdf.putAtt(ncid, globid, 'data_mode', 'D');
netcdf.putAtt(ncid, globid, 'format_version', '1.2')
netcdf.putAtt(ncid, globid, 'area', 'Global_Ocean');

netcdf.putAtt(ncid, globid, 'source', 'research vessel');
netcdf.putAtt(ncid, globid, 'source_platform_category_code','31');
netcdf.putAtt(ncid, globid, 'platform_code','NA');
netcdf.putAtt(ncid, globid, 'platform_name','NA');
netcdf.putAtt(ncid, globid, 'site_code','NA');

netcdf.putAtt(ncid, globid, 'citation', 'These data were collected and made freely available by the Copernicus project and the programs that contribute to it.');
netcdf.putAtt(ncid, globid, 'distribution_statement','These data follow Copernicus standards; they are public and free of charge. User assumes all risk for use of data. User must display citation in any publication or product using data. User must contact PI prior to any commercial use of data.');
netcdf.putAtt(ncid, globid, 'naming_authority', 'Copernicus');

netcdf.putAtt(ncid, globid, 'geospatial_lat_min', '-90');
netcdf.putAtt(ncid, globid, 'geospatial_lat_max', '90');
netcdf.putAtt(ncid, globid, 'geospatial_lat_units', 'degree_north');
netcdf.putAtt(ncid, globid, 'geospatial_lon_min', '-180');
netcdf.putAtt(ncid, globid, 'geospatial_lon_max', '180');
netcdf.putAtt(ncid, globid, 'geospatial_lon_units', 'degree_east');
netcdf.putAtt(ncid, globid, 'geospatial_vertical_min', num2str(min(depthin)));
netcdf.putAtt(ncid, globid, 'geospatial_vertical_max', num2str(max(depthin)));
netcdf.putAtt(ncid, globid, 'geospatial_vertical_positive', 'down');
netcdf.putAtt(ncid, globid, 'geospatial_vertical_units', 'meter');
netcdf.putAtt(ncid, globid, 'time_coverage_start', '1972-01-01');
netcdf.putAtt(ncid, globid, 'time_coverage_end', '2013-01-07');
%netcdf.putAtt(ncid, globid, 'last_date_observation', datestr(timeinglodap(end),'yyyy-mm-ddTHH:MM:SSZ'));
%netcdf.putAtt(ncid, globid, 'last_latitude_observation', latin(end));
%netcdf.putAtt(ncid, globid, 'last_longitude_observation', lonin(end));

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

netcdf.putAtt(ncid, globid, 'qc_manual', 'Lauvset et. al (2016), ESSD');
netcdf.putAtt(ncid, globid, 'quality_control_indicator', '6');
netcdf.putAtt(ncid, globid, 'quality_index', 'A');
%    netcdf.putAtt(ncid, globid, 'SOCAT_flag', '');

%    netcdf.putAtt(ncid, globid, 'summary', '');
netcdf.putAtt(ncid, globid, 'references', 'http://marine.copernicus.eu https://www.glodap.info/');
netcdf.putAtt(ncid, globid, 'comment', ['Further details can be found in Lauvset et al. (2016).',...
    'A new global interior ocean mapped climatology: the 1°x1° GLODAP version 2, ',...
    'Earth Syst. Sci. Data, 8, 325–340, 2016, doi:10.5194/essd-8-325-2016']);

% Close nc file
netcdf.close(ncid);
%end







