clc;
clear all;
%clearvars -except SOCATv6
close all
workrootdir=...
    '/Users/rpr061/Documents/localtestarea/CARBON-REP-042020/';

% input folder
outdir=[workrootdir,'SOCAT_REP_OBSERVATIONS/'];%['/Users/rpr061/Documents/localtestarea/CARBON-REP-042020/'];

% Read synthesis, all or both?
whichsocat='all'; % 'all','synthesis','both'

% Read SOCATv2019 synthesis file
if strcmp(whichsocat,'synthesis') | strcmp(whichsocat,'both')
    synthfile=[workrootdir,'SOCATv2019.tsv'];
    zipsynthfile=[synthfile,'.zip'];
    
    if exist([workrootdir,'SOCATV2019_synth.mat'],'file')
        load([workrootdir,'SOCATV2019_synth.mat']);
    else % read file
        if ~exist(synthfile, 'file')
            tmpstr=['unzip ',zipsynthfile,' -d ', workrootdir];
            system(tmpstr);
            clear tmpstr
        end
        SOCATv2019=readSOCATenhancedfile(synthfile);
        save([workrootdir,'SOCATv2019.mat'], 'SOCATv2019');
    end
end

% Read SOCATv2109 enhanced files (A-E, WOCE 2-4)
if strcmp(whichsocat,'all') | strcmp(whichsocat,'both')
    if exist([workrootdir,'SOCATV2019all.mat'],'file')
        load([workrootdir,'SOCATV2019all.mat']);
    else % read files
        allcruisefiles=dir([workrootdir,'SOCATv2019v7All_ABCDE_enhanced_datafiles/**/*.tsv']);
        for acf=1:length(allcruisefiles);
            disp(acf);
            cruisefile=[allcruisefiles(acf).folder,'/',allcruisefiles(acf).name];
            cruisedata=readSOCATenhancedfile(cruisefile);
            if acf==1; SOCATv2019all=cruisedata;
            else SOCATv2019all=CatStructFields(SOCATv2019all,cruisedata,1);
            end
            clear cruisedata
        end
        save([workrootdir,'SOCATv2019all.mat'], 'SOCATv2019all', '-v7.3');
    end
end



%%
SOCATv2019=SOCATv2019all;
% From categorical to cellstr. If done when creating SOCATv2019.mat, the variables is >2GB
SOCATv2019.Expocode=cellstr(SOCATv2019.Expocode);
% In SOCAT tsv NaN are NaN
% CHANGE LONGITUDE TO ±180!!!
SOCATv2019all.longitudedecdegE(SOCATv2019all.longitudedecdegE>180)=...
    SOCATv2019all.longitudedecdegE(SOCATv2019all.longitudedecdegE>180)-360;
unique_expocodes=unique(SOCATv2019.Expocode);
% Read the info sheet
SOCATv2019_info=tdfread([workrootdir,'CMEMS_SOCATv2019.tsv']);
% Convert all char to cellstr
fn=fieldnames(SOCATv2019_info);
for f=1:length(fn)
    if ischar(SOCATv2019_info.(fn{f})); SOCATv2019_info.(fn{f})=cellstr(SOCATv2019_info.(fn{f}));
    elseif isnumeric(SOCATv2019_info.(fn{f})); SOCATv2019_info.(fn{f})=cellstr(num2str(SOCATv2019_info.(fn{f})));
    end
end
nocallsign=cellfun(@isempty,SOCATv2019_info.CallSign_WMO);
allplatformcodes=SOCATv2019_info.CallSign_WMO;
allplatformcodes(nocallsign)=SOCATv2019_info.Name(nocallsign);
unique_platform_code=unique(allplatformcodes);

%% Create NetCDF
cmode = netcdf.getConstant('NETCDF4');
cmode = bitor(cmode,netcdf.getConstant('CLASSIC_MODEL'));

tic

for ec=1:length(unique_platform_code)
    
    currentpc=unique_platform_code{ec};
    filterec=ismember(allplatformcodes,{currentpc});
    filterdata=ismember(SOCATv2019.Expocode,SOCATv2019_info.Expocode(filterec));
    
    % Platform-related info
    platform_code=char(unique(allplatformcodes(filterec)));
    platform_code=platform_code(isstrprop(platform_code, 'alphanum')); % remove not-allowed characters
    platform_name=char(unique(SOCATv2019_info.Name(filterec)));
    wmo_platform_code=char(unique(SOCATv2019_info.CallSign_WMO(filterec)));
    ices_platform_code=char(unique(SOCATv2019_info.ICEScode(filterec))); 
    institution=strjoin(transpose(unique(SOCATv2019_info.Institution(filterec))),'/');
    institution_edmo_code=regexprep(char(join(transpose(unique(SOCATv2019_info.EDMO(filterec))))),' +',' ');
    source_platform_category_code=char(unique(SOCATv2019_info.PlatformType(filterec)));
    switch source_platform_category_code
        case '31'
            source='research vessel';
            data_type='OceanSITES trajectory data';
            cdm_data_type='trajectory';
            outputfile = ['VESSEL/GL_TS_TS_',platform_code,'_SOCATv2019.nc'];
        case '32'
            source='vessel of opportunity';
            data_type='OceanSITES trajectory data';
            cdm_data_type='trajectory';
            outputfile = ['VESSEL/GL_TS_TS_',platform_code,'_SOCATv2019.nc'];
        case '41'
            source='moored surface buoy';
            data_type='OceanSITES time-series data';
            cdm_data_type='station';
            outputfile = ['MOORING/GL_TS_MO_',platform_code,'_SOCATv2019.nc'];
        case '42'
            source='drifting surface float';
            data_type='OceanSITES trajectory data';
            cdm_data_type='trajectory';
            outputfile = ['DRIFTER/GL_TS_DB_',platform_code,'_SOCATv2019.nc'];
    end
    
    % Extract info from SOCAT
    latin = SOCATv2019.latitudedecdegN(filterdata);
    lonin = SOCATv2019.longitudedecdegE(filterdata);
    timeinsocat = datenum(double([SOCATv2019.yr(filterdata),...
        SOCATv2019.mon(filterdata),SOCATv2019.day1(filterdata),...
        SOCATv2019.hh(filterdata),SOCATv2019.mm(filterdata),...
        SOCATv2019.ss(filterdata)]));
    timein = daysdif(datenum('1950-01-01T00:00:00Z','yyyy-mm-ddTHH:MM:SSZ'), ...
        timeinsocat);
    timeinqc=ones(size(timein));
    
    depthin=5.0*ones(size(timein));
    depthinqc=7*ones(size(timein)); % nominal depth
    posinqc=ones(size(timein));
    
    fco2in = round(SOCATv2019.fCO2recuatm(filterdata)*1000.);
    fco2in(isnan(fco2in))=-2147483647;
    % Map QC flags
    fco2inqc1 = SOCATv2019.fCO2rec_flag(filterdata);
    fco2inqc=fco2inqc1;
    fco2inqc(fco2inqc1==2)=1;
    fco2inqc(fco2inqc1==3)=2;
    fco2inqc(fco2inqc1==4)=4;
    fco2inqc(fco2inqc1==9)=9;
    
    tempin = round(SOCATv2019.SSTdegC(filterdata)*1000.);
    tempin(isnan(tempin))=-2147483647;
    tempinqc = zeros(size(tempin));
    
    salin = round(SOCATv2019.sal(filterdata)*1000.);
    salin(isnan(salin))=-2147483647;
    salinqc = zeros(size(salin));
    
    % Create NetCDF file
    ncid = netcdf.create([outdir,outputfile], cmode);
    
    % Define dimensions and variables
    % Define dimensions
    latdim = netcdf.defDim(ncid, 'LATITUDE', length(latin));
    londim = netcdf.defDim(ncid, 'LONGITUDE', length(lonin));
    posdim = netcdf.defDim(ncid, 'POSITION', length(lonin));
    timedim = netcdf.defDim(ncid, 'TIME', length(timein));
    depthdim = netcdf.defDim(ncid, 'DEPTH', 1);
    
    % Define dimension variables
    latvar = netcdf.defVar(ncid, 'LATITUDE', 'NC_FLOAT', latdim);
    lonvar = netcdf.defVar(ncid, 'LONGITUDE', 'NC_FLOAT', londim);
    timevar = netcdf.defVar(ncid, 'TIME', 'NC_DOUBLE', timedim);
    timevarqc = netcdf.defVar(ncid, 'TIME_QC', 'NC_BYTE', timedim);
    depthvar = netcdf.defVar(ncid, 'DEPH', 'NC_FLOAT', [depthdim, timedim]);
    depthvarqc = netcdf.defVar(ncid, 'DEPH_QC', 'NC_BYTE', [depthdim, timedim]);
    posvarqc = netcdf.defVar(ncid, 'POSITION_QC', 'NC_BYTE', posdim);
    
    % Define other variables
    tempvar = netcdf.defVar(ncid, 'TEMP', 'NC_INT', [depthdim, timedim]);
    salvar= netcdf.defVar(ncid, 'PSAL', 'NC_INT', [depthdim, timedim]);
    fco2var= netcdf.defVar(ncid, 'FCO2', 'NC_INT', [depthdim, timedim]);
    
    tempvarqc = netcdf.defVar(ncid, 'TEMP_QC', 'NC_BYTE', [depthdim, timedim]);
    salvarqc = netcdf.defVar(ncid, 'PSAL_QC', 'NC_BYTE', [depthdim, timedim]);
    fco2varqc = netcdf.defVar(ncid, 'FCO2_QC', 'NC_BYTE', [depthdim, timedim]);
    
    % Define _FillValue of variables (does not work after reDef);
    netcdf.defVarFill(ncid, timevar, false, 9.9692099683868690e+36);
    netcdf.defVarFill(ncid, latvar, false, 9.9692099683868690e+36);
    netcdf.defVarFill(ncid, lonvar, false, 9.9692099683868690e+36);
    netcdf.defVarFill(ncid, depthvar, false, 9.9692099683868690e+36);
    netcdf.defVarFill(ncid, timevarqc, false, -127);
    netcdf.defVarFill(ncid, posvarqc, false, -127);
    netcdf.defVarFill(ncid, depthvarqc, false, -127);
    
    netcdf.defVarFill(ncid, tempvar, false, -2147483647);
    netcdf.defVarFill(ncid, salvar, false, -2147483647);
    netcdf.defVarFill(ncid, fco2var, false, -2147483647);
    
    netcdf.defVarFill(ncid, tempvarqc, false, -127);
    netcdf.defVarFill(ncid, salvarqc, false, -127);
    netcdf.defVarFill(ncid, fco2varqc, false, -127);
    
    % Exit define mode
    netcdf.endDef(ncid);
    
    % Put data
    % Put data in dimensions
    netcdf.putVar(ncid, latvar, latin);
    netcdf.putVar(ncid, lonvar, lonin);
    netcdf.putVar(ncid, timevar, timein);
    netcdf.putVar(ncid, timevarqc, timeinqc);
    netcdf.putVar(ncid, depthvar, depthin);
    netcdf.putVar(ncid, depthvarqc, depthinqc);
    netcdf.putVar(ncid, posvarqc, posinqc);
    
    % Put data in variables
    netcdf.putVar(ncid, tempvar, tempin);
    netcdf.putVar(ncid, salvar, salin);
    netcdf.putVar(ncid, fco2var, fco2in);
    
    netcdf.putVar(ncid, tempvarqc, tempinqc);
    netcdf.putVar(ncid, salvarqc, salinqc);
    netcdf.putVar(ncid, fco2varqc, fco2inqc);
    
    % Re-enter define mode (can be done before???)
    netcdf.reDef(ncid);
    
    % Add attributes to dimensions
    netcdf.putAtt(ncid, latvar, 'long_name', 'Latitude of each location');
    netcdf.putAtt(ncid, latvar, 'standard_name', 'latitude');
    netcdf.putAtt(ncid, latvar, 'units', 'degree_north');
    netcdf.putAtt(ncid, latvar, 'valid_min', -90.0);
    netcdf.putAtt(ncid, latvar, 'valid_max', 90.0);
    netcdf.putAtt(ncid, latvar, 'QC_indicator', 1);
    netcdf.putAtt(ncid, latvar, 'QC_procedure', 1);
    netcdf.putAtt(ncid, latvar, 'axis', 'Y');
    
    % Longitude
    netcdf.putAtt(ncid, lonvar, 'long_name', 'Longitude of each location');
    netcdf.putAtt(ncid, lonvar, 'standard_name', 'longitude');
    netcdf.putAtt(ncid, lonvar, 'units', 'degree_east');
    netcdf.putAtt(ncid, lonvar, 'valid_min', -180.0);
    netcdf.putAtt(ncid, lonvar, 'valid_max', 180.0);
    netcdf.putAtt(ncid, lonvar, 'QC_indicator', 1);
    netcdf.putAtt(ncid, lonvar, 'QC_procedure', 1);
    netcdf.putAtt(ncid, lonvar, 'axis', 'X');
    
    % Time
    netcdf.putAtt(ncid, timevar, 'long_name', 'Time');
    netcdf.putAtt(ncid, timevar, 'standard_name', 'time');
    netcdf.putAtt(ncid, timevar, 'units', 'days since 1950-01-01T00:00:00Z');
    netcdf.putAtt(ncid, timevar, 'valid_min', -90000.);
    netcdf.putAtt(ncid, timevar, 'valid_max', 90000.);
    netcdf.putAtt(ncid, timevar, 'QC_indicator', 1);
    netcdf.putAtt(ncid, timevar, 'QC_procedure', 1);
    netcdf.putAtt(ncid, timevar, 'axis', 'T');
    
    % Depth
    netcdf.putAtt(ncid, depthvar, 'long_name', 'Depth');
    netcdf.putAtt(ncid, depthvar, 'standard_name', 'depth');
    netcdf.putAtt(ncid, depthvar, 'valid_min', -12000.0);
    netcdf.putAtt(ncid, depthvar, 'valid_max', 12000.0);
    netcdf.putAtt(ncid, depthvar, 'units', 'm');
    netcdf.putAtt(ncid, depthvar, 'positive', 'down');
    netcdf.putAtt(ncid, depthvar, 'reference', 'sea_level');
    netcdf.putAtt(ncid, depthvar, 'QC_indicator', 1);
    netcdf.putAtt(ncid, depthvar, 'axis', 'Z');
    
    % Temperature
    netcdf.putAtt(ncid, tempvar, 'standard_name','sea_water_temperature');
    netcdf.putAtt(ncid, tempvar, 'units','degrees_C');
    netcdf.putAtt(ncid, tempvar, 'add_offset',0);
    netcdf.putAtt(ncid, tempvar, 'long_name','Sea temperature');
    netcdf.putAtt(ncid, tempvar, 'scale_factor',0.001);
    netcdf.putAtt(ncid, tempvar, 'valid_min', -100.0);
    netcdf.putAtt(ncid, tempvar, 'valid_max', 100.0);
    netcdf.putAtt(ncid, tempvar, 'coordinates', 'TIME LATITUDE LONGITUDE DEPH')
    
    % Salinity
    netcdf.putAtt(ncid, salvar, 'standard_name','sea_water_practical_salinity');
    netcdf.putAtt(ncid, salvar, 'units','0.001');
    netcdf.putAtt(ncid, salvar, 'add_offset',0);
    netcdf.putAtt(ncid, salvar, 'long_name','Practical salinity');
    netcdf.putAtt(ncid, salvar, 'scale_factor',0.001);
    netcdf.putAtt(ncid, salvar, 'valid_min', 0.0);
    netcdf.putAtt(ncid, salvar, 'valid_max', 100.0);
    netcdf.putAtt(ncid, salvar, 'coordinates', 'TIME LATITUDE LONGITUDE DEPH')
    
    % fCO2
    netcdf.putAtt(ncid, fco2var, 'standard_name','fugacity_of_CO2');
    netcdf.putAtt(ncid, fco2var, 'units','µatm');
    netcdf.putAtt(ncid, fco2var, 'add_offset',0);
    netcdf.putAtt(ncid, fco2var, 'long_name','fugacity_of_CO2_in_sea_water');
    netcdf.putAtt(ncid, fco2var, 'scale_factor',0.001);
    netcdf.putAtt(ncid, fco2var, 'valid_min', 0.0);
    netcdf.putAtt(ncid, fco2var, 'valid_max', 100000.0);
    netcdf.putAtt(ncid, fco2var, 'accuracy', '±20 µatm');
    netcdf.putAtt(ncid, fco2var, 'coordinates', 'TIME LATITUDE LONGITUDE DEPH')
    
    % QC variables
    attributes_qc(ncid,timevarqc)
    attributes_qc(ncid,depthvarqc)
    attributes_qc(ncid,posvarqc)
    attributes_qc(ncid,tempvarqc)
    attributes_qc(ncid,salvarqc)
    attributes_qc(ncid,fco2varqc)
    
    % Add global attributes
    globid = netcdf.getConstant('GLOBAL');
    
    netcdf.putAtt(ncid, globid, 'data_type', data_type);
    netcdf.putAtt(ncid, globid, 'platform_code', platform_code);
    netcdf.putAtt(ncid, globid, 'platform_name', platform_name);
    netcdf.putAtt(ncid, globid, 'data_mode', 'D');
    netcdf.putAtt(ncid, globid, 'title', 'Global Ocean - In Situ reprocessed carbon observations - SOCATv2019');
    % netcdf.putAtt(ncid, globid, 'summary', '');
    netcdf.putAtt(ncid, globid, 'naming_authority', 'Copernicus Marine in situ');
    netcdf.putAtt(ncid, globid, 'id', outputfile(1:end-3));
    % netcdf.putAtt(ncid, globid, 'wmo_platform_code', '');
    netcdf.putAtt(ncid, globid, 'source', source);
    netcdf.putAtt(ncid, globid, 'source_platform_category_code', source_platform_category_code);
    netcdf.putAtt(ncid, globid, 'institution_edmo_code', institution_edmo_code);
    netcdf.putAtt(ncid, globid, 'institution', institution);
    
    netcdf.putAtt(ncid, globid, 'geospatial_lat_min', num2str(min(latin)));
    netcdf.putAtt(ncid, globid, 'geospatial_lat_max', num2str(max(latin)));
    netcdf.putAtt(ncid, globid, 'geospatial_lon_min', num2str(min(lonin)));
    netcdf.putAtt(ncid, globid, 'geospatial_lon_max', num2str(max(lonin)));
    netcdf.putAtt(ncid, globid, 'geospatial_vertical_min', '5.0');
    netcdf.putAtt(ncid, globid, 'geospatial_vertical_max', '5.0');
    netcdf.putAtt(ncid, globid, 'time_coverage_start', datestr(min(timeinsocat),'yyyy-mm-ddTHH:MM:SSZ'));
    netcdf.putAtt(ncid, globid, 'time_coverage_end', datestr(max(timeinsocat),'yyyy-mm-ddTHH:MM:SSZ'));
    netcdf.putAtt(ncid, globid, 'cdm_data_type', cdm_data_type);
    %netcdf.putAtt(ncid, globid, 'summary', '');
    
    netcdf.putAtt(ncid, globid, 'format_version', '1.4');
    netcdf.putAtt(ncid, globid, 'Conventions', 'CF-1.6 Copernicus-InSituTAC-Manual-1.0 Copernicus-InSituTAC-SRD-1.4 Copernicus-InSituTAC-ParametersList-3.1.0');
    netcdf.putAtt(ncid, globid, 'netcdf_version', 'netCDF-4 classic model');
    
    netcdf.putAtt(ncid, globid, 'references', 'http://marine.copernicus.eu https://www.socat.info');
    netcdf.putAtt(ncid, globid, 'data_assembly_center', 'BERGEN');
    netcdf.putAtt(ncid, globid, 'update_interval', 'yearly');
    netcdf.putAtt(ncid, globid, 'citation', ['These data were collected and made freely available by the Copernicus project and the programs that contribute to it. ',...
        'Cite as Bakker et al. (2016), Bakker et al. (2019); see dois']);
    netcdf.putAtt(ncid, globid, 'doi', '0.5194/essd-8-383-2016 10.25921/cpbz-qa92');
    netcdf.putAtt(ncid, globid, 'date_update', datestr(now, 'yyyy-mm-ddTHH:MM:SSZ'));
    netcdf.putAtt(ncid, globid, 'history', [datestr(now, 'yyyy-mm-ddTHH:MM:SSZ'),' : Creation']);
    netcdf.putAtt(ncid, globid, 'last_date_observation',datestr(timeinsocat(end),'yyyy-mm-ddTHH:MM:SSZ'));
    netcdf.putAtt(ncid, globid, 'last_latitude_observation',latin(end));
    netcdf.putAtt(ncid, globid, 'last_longitude_observation',lonin(end));
    netcdf.putAtt(ncid, globid, 'distribution_statement','These data follow Copernicus standards; they are public and free of charge. User assumes all risk for use of data. User must display citation in any publication or product using data. User must contact PI prior to any commercial use of data.');
    
    % Close nc file
    netcdf.close(ncid);
    %end
    
    
    clearvars '-except' workrootdir outdir SOCATv2019 SOCATv2019_info ...
        cmode allplatformcodes unique_platform_code ec
    
end

toc

system(['cd /Users/rpr061/Dropbox/BCDC_Projects/CMEMS_INSTAC/REP_Products/current_FormatChecker/; for f in ', outdir,'VESSEL/*.nc; do ./control.csh $f >> ',outdir,'VESSEL/formatcheckoutSOCAT; done'])
system(['cd /Users/rpr061/Dropbox/BCDC_Projects/CMEMS_INSTAC/REP_Products/current_FormatChecker/; for f in ', outdir,'MOORING/*.nc; do ./control.csh $f >> ',outdir,'MOORING/formatcheckoutSOCAT; done'])
system(['cd /Users/rpr061/Dropbox/BCDC_Projects/CMEMS_INSTAC/REP_Products/current_FormatChecker/; for f in ', outdir,'DRIFTER/*.nc; do ./control.csh $f >> ',outdir,'DRIFTER/formatcheckoutSOCAT; done'])

%% ======================================================================
%%==================== FUNCTIONS =========================================
function [outstructure]=readSOCATenhancedfile(infile)

% Get # of header lines to ignore
tmpstr=['grep -n -m 1 ''Expocode''$''\t''''version''$''\t''''SOCAT_DOI'' ', infile,' | cut -d: -f1'];
[status,nheaderlines]=system(tmpstr);
if status~=0; error('Cant grep SOCAT synthesis file'); end
nheaderlines=str2num(nheaderlines);

% Check column headers
fileID = fopen(infile,'r');
colheaders = textscan(fileID, '%s', 1, 'Delimiter', '\n', 'HeaderLines', nheaderlines-1);
fclose(fileID);


if ~strcmp(char(colheaders{1}), sprintf(['Expocode\tversion\tSOCAT_DOI\tQC_Flag\t',...
        'yr\tmon\tday\thh\tmm\tss\tlongitude [dec.deg.E]\tlatitude [dec.deg.N]\t',...
        'sample_depth [m]\tsal\tSST [deg.C]\tTequ [deg.C]\t',...
        'PPPP [hPa]\tPequ [hPa]\tWOA_SSS\tNCEP_SLP [hPa]\tETOPO2_depth [m]\tdist_to_land [km]\t',...
        'GVCO2 [umol/mol]\tfCO2rec [uatm]\tfCO2rec_src\tfCO2rec_flag'])) && ... % header of synthesis file
        ~strcmp(char(colheaders{1}), sprintf(['Expocode\tversion\tSOCAT_DOI\tQC_Flag\t',...
        'yr\tmon\tday\thh\tmm\tss\tlongitude [dec.deg.E]\tlatitude [dec.deg.N]\t',...
        'sample_depth [m]\tsal\tSST [deg.C]\tTequ [deg.C]\t',...
        'PPPP [hPa]\tPequ [hPa]\tWOA_SSS\tNCEP_SLP [hPa]\tETOPO2_depth [m]\tdist_to_land [km]\t',...
        'GVCO2 [umol/mol]\txCO2water_equ_dry [umol/mol]\txCO2water_SST_dry [umol/mol]\t',...
        'pCO2water_equ_wet [uatm]\tpCO2water_SST_wet [uatm]\tfCO2water_equ_wet [uatm]\tfCO2water_SST_wet [uatm]\t',...
        'fCO2rec [uatm]\tfCO2rec_src\tfCO2rec_flag'])) % header of SOCAT-enhanced files
    error('Revise the columns headers and formatSpec');
end

% Read columns of data as text
% For more information, see the TEXTSCAN documentation.

if length(strsplit(char(colheaders{1}), '\t'))<=27 %synthesis file
    formatSpec = '%s%s%s%s%d%d%d%d%d%d%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%d%[^\n\r]';
    ifCO2rec=24;
    ifCO2rec_flag=26;
else % enhanced files
    formatSpec = '%s%s%s%s%d%d%d%d%d%d%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%d%d%[^\n\r]';
    ifCO2rec=30;
    ifCO2rec_flag=32;
end
fileID = fopen(infile,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '\t', 'TextType', 'string', 'HeaderLines' ,nheaderlines, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);

outstructure.Expocode = categorical(dataArray{:, 1});
%outstructure.version1 = (dataArray{:, 2});
%outstructure.SOCAT_DOI = categorical(dataArray{:, 3});
outstructure.QC_Flag = categorical(dataArray{:, 4});
outstructure.yr = (dataArray{:, 5});
outstructure.mon = (dataArray{:, 6});
outstructure.day1 = (dataArray{:, 7});
outstructure.hh = (dataArray{:, 8});
outstructure.mm = (dataArray{:, 9});
outstructure.ss = (dataArray{:, 10});
outstructure.longitudedecdegE = (dataArray{:, 11});
outstructure.latitudedecdegN = (dataArray{:, 12});
%outstructure.sample_depthm = (dataArray{:, 13});
outstructure.sal = (dataArray{:, 14});
outstructure.SSTdegC = (dataArray{:, 15});
%outstructure.TequdegC = (dataArray{:, 16});
%outstructure.PPPPhPa = (dataArray{:, 17});
%outstructure.PequhPa = (dataArray{:, 18});
%outstructure.WOA_SSS = (dataArray{:, 19});
%outstructure.NCEP_SLPhPa = (dataArray{:, 20});
%outstructure.ETOPO2_depthm = (dataArray{:, 21});
%outstructure.dist_to_landkm = (dataArray{:, 22});
outstructure.GVCO2umolmol = (dataArray{:, 23});
outstructure.fCO2recuatm = (dataArray{:, ifCO2rec});
%outstructure.fCO2rec_src = (dataArray{:, 25});
outstructure.fCO2rec_flag = (dataArray{:, ifCO2rec_flag});
end

function S = CatStructFields(S, T, dim)
fields = fieldnames(S);
for k = 1:numel(fields)
    aField     = fields{k}; % EDIT: changed to {}
    S.(aField) = cat(dim, S.(aField), T.(aField));
end; end

function attributes_qc(ncid, inputvar)
netcdf.putAtt(ncid, inputvar, 'long_name', 'quality flag');
netcdf.putAtt(ncid, inputvar, 'conventions', 'Copernicus Marine in situ reference table 2');
netcdf.putAtt(ncid, inputvar, 'valid_min', int8(0));
netcdf.putAtt(ncid, inputvar, 'valid_max', int8(9));
netcdf.putAtt(ncid, inputvar, 'flag_values',  [int8(0)  int8(1)  int8(2)  int8(3)  int8(4)  int8(5)  int8(6)  int8(7)  int8(8)  int8(9)]);
netcdf.putAtt(ncid, inputvar, 'flag_meanings','no_qc_performed good_data probably_good_data bad_data_that_are_potentially_correctable bad_data value_changed not_used nominal_value interpolated_value missing_value');
end
