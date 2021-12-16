clc;
%clear all;
clearvars -except SOCAT
close all
workrootdir=...
    '/Users/rpr061/Documents/localtestarea/CARBON-REP-122021/';

% input folder
outdir=[workrootdir,'SOCAT_REP_OBSERVATIONS/'];%['/Users/rpr061/Documents/localtestarea/CARBON-REP-042020/'];
mkdir(outdir, 'VESSEL')
mkdir(outdir, 'MOORING')
mkdir(outdir, 'DRIFTER')
mkdir(outdir, 'AUV')


%%
% Read synthesis, all or both?
whichsocat='synthesis'; % 'all','synthesis','both'
socatversion='SOCATv2021';
%%
% Read SOCATv2019 synthesis file
%if strcmp(whichsocat,'synthesis') | strcmp(whichsocat,'both')
    synthfile=[workrootdir,[socatversion,'.tsv']];
    zipsynthfile=[synthfile,'.zip'];
    synthfileE=[workrootdir,[socatversion,'_FlagE.tsv']];
    
   if exist('SOCAT','var')
       disp('SOCAT data is already a variable')
    
   %elseif exist([workrootdir,[socatversion,'_synthAE.mat']],'file')
   %     disp('loading the matfile')
   %     load([workrootdir,[socatversion,'_synthAE.mat']]);
    else % read file
        if ~exist(synthfile, 'file')
            tmpstr=['unzip ',zipsynthfile,' -d ', workrootdir];
            system(tmpstr);
            clear tmpstr
        end
        SOCATAD=readSOCATenhancedfile(synthfile);
        SOCATE=readSOCATenhancedfile(synthfileE);
        fields = fieldnames(SOCATAD);
for k = 1:numel(fields)
  aField     = fields{k}; % EDIT: changed to {}
  SOCAT.(aField) = cat(1, SOCATAD.(aField), SOCATE.(aField));
end
        save([workrootdir,[socatversion,'_synthAE.mat']], 'SOCAT', '-v7.3');
    end
%end

%%
%SOCAT=SOCATall;
% From categorical to cellstr. If done when creating SOCAT.mat, the variables is >2GB
SOCAT.Expocode=cellstr(SOCAT.Expocode);
% In SOCAT tsv NaN are NaN
% CHANGE LONGITUDE TO ±180!!!
SOCAT.longitudedecdegE(SOCAT.longitudedecdegE>180)=...
    SOCAT.longitudedecdegE(SOCAT.longitudedecdegE>180)-360;
unique_expocodes=unique(SOCAT.Expocode);
%%
% Read the info sheet
%SOCAT_info=tdfread([workrootdir,'CMEMS_SOCAT.tsv']);
%SOCAT_info=tdfread(['/Users/rpr061/Downloads/SOCATv2020CMEMS.tsv']); %
%TDFread does not do well with UTF8
file=fopen([workrootdir,'SOCATv2021CMEMS.tsv'],'r','n','UTF-8');
header=fgets(file);
header=strsplit(header,'\t');
header{end}=header{end}(1:end-2);
data=textscan(file, [repmat('%s',1,11),'%s\r'], 'Delimiter','\t');
fclose(file)
% Create structure
for ff=1:length(header);
    SOCAT_info.(header{ff})=data{ff}; end

% Convert all char to cellstr
% fn=fieldnames(SOCAT_info);
% for f=1:length(fn)
%     if ischar(SOCAT_info.(fn{f})); SOCAT_info.(fn{f})=cellstr(SOCAT_info.(fn{f}));
%     elseif isnumeric(SOCAT_info.(fn{f})); SOCAT_info.(fn{f})=cellstr(num2str(SOCAT_info.(fn{f})));
%     end
% end
nocallsign=cellfun(@isempty,SOCAT_info.CallSign_WMO);
allplatformcodes=SOCAT_info.CallSign_WMO;
allplatformcodes(nocallsign)=SOCAT_info.Name(nocallsign);
unique_platform_code=unique(allplatformcodes);


%% Create NetCDF
cmode = netcdf.getConstant('NETCDF4');
cmode = bitor(cmode,netcdf.getConstant('CLASSIC_MODEL'));

tic

for ec=1:2%length(unique_platform_code) 
  
    currentpc=unique_platform_code{ec};
    filterec=ismember(allplatformcodes,{currentpc});
    filterdata=ismember(SOCAT.Expocode,SOCAT_info.Expocode(filterec));
    
    % Platform-related info
    platform_code=char(unique(allplatformcodes(filterec)));
    platform_code=platform_code(isstrprop(platform_code, 'alphanum')); % remove not-allowed characters
    platform_name=char(unique(SOCAT_info.Name(filterec)));
    wmo_platform_code=char(unique(SOCAT_info.CallSign_WMO(filterec)));
    if isempty(wmo_platform_code); wmo_platform_code=' '; end
    ices_platform_code=char(unique(SOCAT_info.ICEScode(filterec)));
    if isempty(ices_platform_code); ices_platform_code=' '; end

    % List of institutions
    institution=strjoin(transpose(unique(SOCAT_info.Institution(filterec))),'/');
    institution_edmo_code=regexprep(char(join(transpose(unique(SOCAT_info.EDMO(filterec))))),' +',' ');
        if strcmp(institution_edmo_code, " NaN"); institution_edmo_code=' '; end
    pi_institution=strjoin(transpose(unique(SOCAT_info.PI_Institution(filterec))),'/');
    %pi_institution_edmo_code=regexprep(char(join(transpose(unique(SOCAT_info.PI_EDMO(filterec))))),' +',' ');
     pi_institution_edmo_code=regexprep(char(join(transpose(unique(SOCAT_info.PI_EDMO(filterec))))),' +',' ');
   %allinstitution=[institution,';',pi_institution];
   if ~strcmp(institution_edmo_code,pi_institution_edmo_code)
           allinstitution=[institution,' ',pi_institution];
    allinstitution_edmo_code=[institution_edmo_code,' ',pi_institution_edmo_code];
   else allinstitution_edmo_code=[institution_edmo_code];
           allinstitution=[institution];
   end
   if strcmp(allinstitution_edmo_code(1),' '); allinstitution_edmo_code=allinstitution_edmo_code(2:end); end 
   
   
    % List of PIs
    allpis=unique(SOCAT_info.PI(filterec));
    for aa=1:length(allpis)
    indpi{aa}=strsplit(allpis{aa},';'); end
    pis=strjoin(unique([indpi{:}]),';'); 
    %pis=strjoin(transpose(unique(SOCAT_info.PI(filterec))),';');
    
    source_platform_category_code=char(unique(SOCAT_info.PlatformType(filterec)));
    %if ~strcmp(source_platform_category_code,'3B'); continue; end

    switch source_platform_category_code
        case '31'
            source='research vessel';
            data_type='OceanSITES trajectory data';
            cdm_data_type='trajectory';
            outputfile = ['VESSEL/GL_TS_CO_',platform_code,'-',socatversion,'.nc'];
                                    sensormount='mounted_on_shipborne_fixed';

        case '32'
            source='vessel of opportunity';
            data_type='OceanSITES trajectory data';
            cdm_data_type='trajectory';
            outputfile = ['VESSEL/GL_TS_CO_',platform_code,'-',socatversion,'.nc'];
                        sensormount='mounted_on_shipborne_fixed';


        case '3B'
            source='autonomous surface water vehicle';
            data_type='OceanSITES trajectory data';
            cdm_data_type='trajectory';
            if contains(platform_name, 'SailDrone');
            outputfile = ['AUV/GL_TS_SD_',platform_code,'-',socatversion,'.nc'];
            elseif contains(platform_name, 'WaveGlider')
            outputfile = ['AUV/GL_TS_GL_',platform_code,'-',socatversion,'.nc'];
            end
                        sensormount='mounted_on_glider';                        
                        
        case '41'
            source='moored surface buoy';
            data_type='OceanSITES time-series data';
            cdm_data_type='timeSeries';
            outputfile = ['MOORING/GL_TS_MO_',platform_code,'-',socatversion,'.nc'];
                        sensormount='mounted_on_surface_buoy';

        case '42'
            source='drifting surface float';
            data_type='OceanSITES trajectory data';
            cdm_data_type='trajectory';
            outputfile = ['DRIFTER/GL_TS_DB_',platform_code,'-',socatversion,'.nc'];
            sensormount='mounted_on_glider';
    end
    splitoutfile=split(outputfile,'/');
    identifier=splitoutfile{2}(1:end-3);
    
    % Extract info from SOCAT
    latin = SOCAT.latitudedecdegN(filterdata);
    lonin = SOCAT.longitudedecdegE(filterdata);
    timeinsocat = datenum(double([SOCAT.yr(filterdata),...
        SOCAT.mon(filterdata),SOCAT.day1(filterdata),...
        SOCAT.hh(filterdata),SOCAT.mm(filterdata),...
        SOCAT.ss(filterdata)]));
    timein = timeinsocat - datenum('1950-01-01T00:00:00Z','yyyy-mm-ddTHH:MM:SSZ');
    timeinqc=ones(size(timein));
    
    depthin=5.0*ones(size(timein));
    depthinqc=7*ones(size(timein)); % nominal depth
    posinqc=ones(size(timein));
    
    fco2in = round(SOCAT.fCO2recuatm(filterdata)*1000.);
    fco2in(isnan(fco2in))=-2147483647;
    % Map QC flags
    fco2inqc1 = SOCAT.fCO2rec_flag(filterdata);
    fco2inqc=fco2inqc1;
    fco2inqc(fco2inqc1==2)=1;
    fco2inqc(fco2inqc1==3)=2;
    fco2inqc(fco2inqc1==4)=4;
    fco2inqc(fco2inqc1==9)=9;
    
    tempin = round(SOCAT.SSTdegC(filterdata)*1000.);
    tempin(isnan(tempin))=-2147483647;
    tempinqc = zeros(size(tempin));
    
    salin = round(SOCAT.sal(filterdata)*1000.);
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
    netcdf.putAtt(ncid, latvar, 'QC_indicator', int8(1));
    netcdf.putAtt(ncid, latvar, 'uncertainty', '');
    netcdf.putAtt(ncid, latvar, 'comment', '');
    netcdf.putAtt(ncid, latvar, 'axis', 'Y');
    netcdf.putAtt(ncid, latvar, 'ancillary_variables', 'POSITION_QC'); 
    netcdf.delAtt(ncid, latvar, '_FillValue'); 

    
    % Longitude
    netcdf.putAtt(ncid, lonvar, 'long_name', 'Longitude of each location');
    netcdf.putAtt(ncid, lonvar, 'standard_name', 'longitude');
    netcdf.putAtt(ncid, lonvar, 'units', 'degree_east');
    netcdf.putAtt(ncid, lonvar, 'valid_min', -180.0);
    netcdf.putAtt(ncid, lonvar, 'valid_max', 180.0);
    netcdf.putAtt(ncid, lonvar, 'QC_indicator', int8(1));
    netcdf.putAtt(ncid, lonvar, 'uncertainty', '');
    netcdf.putAtt(ncid, lonvar, 'comment', '');
    netcdf.putAtt(ncid, lonvar, 'axis', 'X');
    netcdf.putAtt(ncid, lonvar, 'ancillary_variables', 'POSITION_QC');
        netcdf.delAtt(ncid, lonvar, '_FillValue'); 

    % Time
    netcdf.putAtt(ncid, timevar, 'long_name', 'Time');
    netcdf.putAtt(ncid, timevar, 'standard_name', 'time');
    netcdf.putAtt(ncid, timevar, 'units', 'days since 1950-01-01T00:00:00Z');
    netcdf.putAtt(ncid, timevar, 'valid_min', -90000.);
    netcdf.putAtt(ncid, timevar, 'valid_max', 90000.);
    netcdf.putAtt(ncid, timevar, 'QC_indicator', int8(1));
    netcdf.putAtt(ncid, timevar, 'uncertainty', '');
    netcdf.putAtt(ncid, timevar, 'comment', '');
    netcdf.putAtt(ncid, timevar, 'axis', 'T');
    netcdf.putAtt(ncid, timevar, 'ancillary_variables', 'TIME_QC');
    netcdf.putAtt(ncid, timevar, 'calendar', 'standard');
    netcdf.delAtt(ncid, timevar, '_FillValue'); 

    
    % Depth
    netcdf.putAtt(ncid, depthvar, 'long_name', 'Depth');
    netcdf.putAtt(ncid, depthvar, 'standard_name', 'depth');
    netcdf.putAtt(ncid, depthvar, 'units', 'm');
    netcdf.putAtt(ncid, depthvar, 'positive', 'down');
    netcdf.putAtt(ncid, depthvar, 'valid_min', -12000.0);
    netcdf.putAtt(ncid, depthvar, 'valid_max', 12000.0);
    netcdf.putAtt(ncid, depthvar, 'uncertainty', '');
    netcdf.putAtt(ncid, depthvar, 'comment', '');
    netcdf.putAtt(ncid, depthvar, 'axis', 'Z');
    netcdf.putAtt(ncid, depthvar, 'reference', 'sea_level');
    netcdf.putAtt(ncid, depthvar, 'data_mode', 'D');
    netcdf.putAtt(ncid, depthvar, 'ancillary_variables', 'DEPH_QC');

    
    % Temperature
    netcdf.putAtt(ncid, tempvar, 'standard_name','sea_water_temperature');
    netcdf.putAtt(ncid, tempvar, 'units','degrees_C');
    netcdf.putAtt(ncid, tempvar, 'add_offset',0);
    netcdf.putAtt(ncid, tempvar, 'scale_factor',0.001);
    netcdf.putAtt(ncid, tempvar, 'long_name','Sea temperature');
    netcdf.putAtt(ncid, tempvar, 'valid_min', -100.0);
    netcdf.putAtt(ncid, tempvar, 'valid_max', 100.0);
    %netcdf.putAtt(ncid, tempvar, 'comment', '');
    %netcdf.putAtt(ncid, tempvar, 'uncertainty', '');
    %netcdf.putAtt(ncid, tempvar, 'accuracy', '');
    %netcdf.putAtt(ncid, tempvar, 'precision', '');
    %netcdf.putAtt(ncid, tempvar, 'resolution', '');
    %netcdf.putAtt(ncid, tempvar, 'cell_methods', '');
    netcdf.putAtt(ncid, tempvar, 'coordinates', 'TIME LATITUDE LONGITUDE DEPH')
    %netcdf.putAtt(ncid, tempvar, 'type_of_analysis', '');
    %netcdf.putAtt(ncid, tempvar, 'sensor_depth', '');
    netcdf.putAtt(ncid, tempvar, 'sensor_mount', sensormount);
    %netcdf.putAtt(ncid, tempvar, 'sensor_orientation', '');
    netcdf.putAtt(ncid, tempvar, 'data_mode', 'D');
    netcdf.putAtt(ncid, tempvar, 'ancillary_variables', 'TEMP_QC');

    
    % Salinity
    netcdf.putAtt(ncid, salvar, 'standard_name','sea_water_practical_salinity');
    netcdf.putAtt(ncid, salvar, 'units','0.001');
    netcdf.putAtt(ncid, salvar, 'add_offset',0);
    netcdf.putAtt(ncid, salvar, 'scale_factor',0.001);
    netcdf.putAtt(ncid, salvar, 'long_name','Practical salinity');
    netcdf.putAtt(ncid, salvar, 'valid_min', 0.0);
    netcdf.putAtt(ncid, salvar, 'valid_max', 100.0);
    netcdf.putAtt(ncid, salvar, 'coordinates', 'TIME LATITUDE LONGITUDE DEPH')
    %netcdf.putAtt(ncid, salvar, 'comment', '');
    %netcdf.putAtt(ncid, salvar, 'uncertainty', '');
    %netcdf.putAtt(ncid, salvar, 'accuracy', '');
    %netcdf.putAtt(ncid, salvar, 'precision', '');
    %netcdf.putAtt(ncid, salvar, 'resolution', '');
    %netcdf.putAtt(ncid, salvar, 'cell_methods', '');
    netcdf.putAtt(ncid, salvar, 'coordinates', 'TIME LATITUDE LONGITUDE DEPH')
    %netcdf.putAtt(ncid, salvar, 'type_of_analysis', '');
    %netcdf.putAtt(ncid, salvar, 'sensor_depthh', '');
    %netcdf.putAtt(ncid, salvar, 'sensor_mount', sensormount);
    %netcdf.putAtt(ncid, salvar, 'sensor_orientation', '');
    netcdf.putAtt(ncid, salvar, 'data_mode', 'D');
    netcdf.putAtt(ncid, salvar, 'ancillary_variables', 'PSAL_QC');
    
    % fCO2
    netcdf.putAtt(ncid, fco2var, 'standard_name','fugacity_of_carbon_dioxide_in_sea_water');
    netcdf.putAtt(ncid, fco2var, 'units','µatm');
    netcdf.putAtt(ncid, fco2var, 'add_offset',0);
    netcdf.putAtt(ncid, fco2var, 'scale_factor',0.001);
    netcdf.putAtt(ncid, fco2var, 'long_name','CO2 fugacity');
    netcdf.putAtt(ncid, fco2var, 'valid_min', 0.0);
    netcdf.putAtt(ncid, fco2var, 'valid_max', 100000.0);
    %netcdf.putAtt(ncid, fco2var, 'comment', '');
    %netcdf.putAtt(ncid, fco2var, 'uncertainty', '');
    netcdf.putAtt(ncid, fco2var, 'accuracy', '±20 µatm');
    %netcdf.putAtt(ncid, fco2var, 'precision', '');
    %netcdf.putAtt(ncid, fco2var, 'resolution', '');
    %netcdf.putAtt(ncid, fco2var, 'cell_methods', '');
    netcdf.putAtt(ncid, fco2var, 'coordinates', 'TIME LATITUDE LONGITUDE DEPH')
    %netcdf.putAtt(ncid, fco2var, 'type_of_analysis', '');
    %netcdf.putAtt(ncid, fco2var, 'sensor_depth', '');
    netcdf.putAtt(ncid, fco2var, 'sensor_mount', sensormount);
    %netcdf.putAtt(ncid, fco2var, 'sensor_orientation', '');
    netcdf.putAtt(ncid, fco2var, 'data_mode', 'D');
    netcdf.putAtt(ncid, fco2var, 'ancillary_variables', 'FCO2_QC');    
    
    % QC variables
    attributes_qc(ncid,timevarqc); netcdf.putAtt(ncid, timevarqc, 'long_name', 'Time quality flag');
    attributes_qc(ncid,depthvarqc); netcdf.putAtt(ncid, depthvarqc, 'long_name', 'Depth quality flag');
    attributes_qc(ncid,posvarqc); netcdf.putAtt(ncid, posvarqc, 'long_name', 'Position quality flag');
    attributes_qc(ncid,tempvarqc); netcdf.putAtt(ncid, tempvarqc, 'long_name', 'Sea temperature quality flag');
    attributes_qc(ncid,salvarqc); netcdf.putAtt(ncid, salvarqc, 'long_name', 'Practical salinity quality flag');
    attributes_qc(ncid,fco2varqc); netcdf.putAtt(ncid, fco2varqc, 'long_name', 'CO2 fugacity quality flag');
    
    % Add global attributes
    globid = netcdf.getConstant('GLOBAL');
    
    netcdf.putAtt(ncid, globid, 'platform_code', platform_code);
    netcdf.putAtt(ncid, globid, 'platform_name', platform_name);
    netcdf.putAtt(ncid, globid, 'data_mode', 'D');
    netcdf.putAtt(ncid, globid, 'title', 'Global Ocean - In Situ reprocessed carbon observations - SOCATv2021');
    netcdf.putAtt(ncid, globid, 'summary', ['Observations of temperature, salinity and',...
        ' fCO2 from the platform ',platform_name,' included in SOCATv2021. ',...
        'The Surface Ocean CO2 Atlas (SOCAT)',...
        ' is a synthesis activity for quality-controlled, surface ocean fCO2',...
        ' (fugacity of carbon dioxide) observations by the international marine',...
        ' carbon research community (>100 contributors). SOCAT enables quantification',...
        'of the ocean carbon sink and ocean acidification and evaluation of ocean',...
        ' biogeochemical models. SOCAT is a core Global Ocean Observing System',...
        ' data product for biogeochemistry endorsed by the Global Ocean Observing System GOOS.']);
    netcdf.putAtt(ncid, globid, 'naming_authority', 'Copernicus Marine In Situ');
    netcdf.putAtt(ncid, globid, 'id', identifier);
    
    netcdf.putAtt(ncid, globid, 'wmo_platform_code', wmo_platform_code);
    netcdf.putAtt(ncid, globid, 'ices_platform_code', ices_platform_code);
    
    %     if (strcmp(source_platform_category_code,'41') | strcmp(source_platform_category_code,'42')) & ~isempty(wmo_platform_code)
%         netcdf.putAtt(ncid, globid, 'wmo_platform_code', wmo_platform_code);
%         netcdf.putAtt(ncid, globid, 'ices_platform_code', ' ');
%     elseif (strcmp(source_platform_category_code,'31') | strcmp(source_platform_category_code,'32')) & ~isempty(ices_platform_code)
%         netcdf.putAtt(ncid, globid, 'ices_platform_code', ices_platform_code);
%         netcdf.putAtt(ncid, globid, 'wmo_platform_code', ' ');
%     end
    netcdf.putAtt(ncid, globid, 'source', source);
    netcdf.putAtt(ncid, globid, 'source_platform_category_code', source_platform_category_code);
    netcdf.putAtt(ncid, globid, 'institution_edmo_code', allinstitution_edmo_code);
    netcdf.putAtt(ncid, globid, 'institution', allinstitution);
    netcdf.putAtt(ncid, globid, 'institution_references', '');
    netcdf.putAtt(ncid, globid, 'site_code', ''); %OceanSites code; e.g. STATION-M
    netcdf.putAtt(ncid, globid, 'comment', ''); 
    netcdf.putAtt(ncid, globid, 'contact', 'data.manager@bcdc.no cmems-service@ifremer.fr'); 

    netcdf.putAtt(ncid, globid, 'area', 'Global');     
    netcdf.putAtt(ncid, globid, 'geospatial_lat_min', num2str(min(latin)));
    netcdf.putAtt(ncid, globid, 'geospatial_lat_max', num2str(max(latin)));
    netcdf.putAtt(ncid, globid, 'geospatial_lon_min', num2str(min(lonin)));
    netcdf.putAtt(ncid, globid, 'geospatial_lon_max', num2str(max(lonin)));
    netcdf.putAtt(ncid, globid, 'geospatial_vertical_min', '5.0');
    netcdf.putAtt(ncid, globid, 'geospatial_vertical_max', '5.0');
    netcdf.putAtt(ncid, globid, 'time_coverage_start', datestr(min(timeinsocat),'yyyy-mm-ddTHH:MM:SSZ'));
    netcdf.putAtt(ncid, globid, 'time_coverage_end', datestr(max(timeinsocat),'yyyy-mm-ddTHH:MM:SSZ'));
    netcdf.putAtt(ncid, globid, 'cdm_data_type', cdm_data_type);
    netcdf.putAtt(ncid, globid, 'data_type', data_type);
    netcdf.putAtt(ncid, globid, 'bottom_depth', ' ');

    netcdf.putAtt(ncid, globid, 'format_version', '1.4');
    netcdf.putAtt(ncid, globid, 'Conventions', 'CF-1.6 Copernicus-InSituTAC-FormatManual-1.42 Copernicus-InSituTAC-SRD-1.5 Copernicus-InSituTAC-ParametersList-3.2.0');
    netcdf.putAtt(ncid, globid, 'netcdf_version', 'netCDF-4 classic model');
    netcdf.putAtt(ncid, globid, 'references', 'http://marine.copernicus.eu http://www.marineinsitu.eu https://www.socat.info');
    netcdf.putAtt(ncid, globid, 'data_assembly_center', 'University of Bergen');
    netcdf.putAtt(ncid, globid, 'update_interval', 'P1Y');
    
    netcdf.putAtt(ncid, globid, 'citation', ['These data were collected and made freely available by the Copernicus project and the programs that contribute to it. ',...
        'The Surface Ocean CO2 Atlas (SOCAT) is an international effort, endorsed by the ',...
        'International Ocean Carbon Coordination Project (IOCCP), the Surface Ocean ',...
        'Lower Atmosphere Study (SOLAS) and the Integrated Marine Biosphere Research ',...
        '(IMBeR) program, to deliver a uniformly quality-controlled surface ocean CO2 ',...
        'database. The many researchers and funding agencies responsible for the collection ',...
        'of data and quality control are thanked for their contributions to SOCAT. ',...
        'Cite as Bakker et al. (2016), Bakker et al. (2021).']);
    netcdf.putAtt(ncid, globid, 'distribution_statement',['These data follow Copernicus standards; they are public ',...
        'and free of charge. User assumes all risk for use of data. User must display citation in any publication or ',...
        'product using data. User must contact PI prior to any commercial use of data.']);
    
    netcdf.putAtt(ncid, globid, 'doi', 'https://doi.org/10.5194/essd-8-383-2016 https://doi.org/10.25921/yg69-jd96');
    netcdf.putAtt(ncid, globid, 'pi_name', pis);
    netcdf.putAtt(ncid, globid, 'qc_manual', '');
    
    netcdf.putAtt(ncid, globid, 'date_update', datestr(now, 'yyyy-mm-ddTHH:MM:SSZ'));
    netcdf.putAtt(ncid, globid, 'history', [datestr(now, 'yyyy-mm-ddTHH:MM:SSZ'),' : Creation']);
    netcdf.putAtt(ncid, globid, 'last_date_observation',datestr(timeinsocat(end),'yyyy-mm-ddTHH:MM:SSZ'));
    netcdf.putAtt(ncid, globid, 'last_latitude_observation',latin(end));
    netcdf.putAtt(ncid, globid, 'last_longitude_observation',lonin(end));
    netcdf.putAtt(ncid, globid, 'wmo_inst_type', '');

    
%    netcdf.putAtt(ncid, globid, 'PI_institution_edmo_code', pi_institution_edmo_code);
%    netcdf.putAtt(ncid, globid, 'PI_institution', pi_institution);
        
    
    
    % Close nc file
    netcdf.close(ncid);    
    
    clearvars '-except' workrootdir outdir SOCAT SOCAT_info socatversion ...
        cmode allplatformcodes unique_platform_code ec
    
end

toc
%%
%system(['cd /Users/rpr061/Dropbox/BCDC_projects/CMEMS_INSTAC/Releases/current_FormatChecker; for f in ', outdir,'VESSEL/*.nc; do ./control.csh $f >> ',outdir,'VESSEL_formatcheckoutSOCAT; done'])
%system(['cd /Users/rpr061/Dropbox/BCDC_Projects/CMEMS_INSTAC/Releases/current_FormatChecker/; for f in ', outdir,'MOORING/*.nc; do ./control.csh $f >> ',outdir,'MOORING_formatcheckoutSOCAT; done'])
%system(['cd /Users/rpr061/Dropbox/BCDC_Projects/CMEMS_INSTAC/Releases/current_FormatChecker/; for f in ', outdir,'DRIFTER/*.nc; do ./control.csh $f >> ',outdir,'DRIFTER_formatcheckoutSOCAT; done'])
%system(['cd /Users/rpr061/Dropbox/BCDC_Projects/CMEMS_INSTAC/Releases/current_FormatChecker/; for f in ', outdir,'AUV/*.nc; do ./control.csh $f >> ',outdir,'AUV_formatcheckoutSOCAT; done'])

%% ======================================================================
%%==================== FUNCTIONS =========================================
function [outstructure]=readSOCATenhancedfile(infile)

% Get # of header lines to ignore
if str2num(infile(end-7:end-4))<2020;
tmpstr=['grep -n -m 1 ''Expocode''$''\t''''version''$''\t''''SOCAT_DOI'' ', infile,' | cut -d: -f1'];
else
    tmpstr=['grep -n -m 1 ''Expocode''$''\t''''version''$''\t''''Source_DOI'' ', infile,' | cut -d: -f1'];
end
[status,nheaderlines]=system(tmpstr);
if status~=0; error('Cant grep SOCAT synthesis file'); end
if isempty(nheaderlines); error('Number of header lines not valid'); end
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
        'fCO2rec [uatm]\tfCO2rec_src\tfCO2rec_flag'])) && ... % header of SOCAT-enhanced files 
         ~strcmp(char(colheaders{1}{1}), sprintf(['Expocode\tversion\tSource_DOI\tQC_Flag\t',...
        'yr\tmon\tday\thh\tmm\tss\tlongitude [dec.deg.E]\tlatitude [dec.deg.N]\t',...
        'sample_depth [m]\tsal\tSST [deg.C]\tTequ [deg.C]\t',...
        'PPPP [hPa]\tPequ [hPa]\tWOA_SSS\tNCEP_SLP [hPa]\tETOPO2_depth [m]\tdist_to_land [km]\t',...
        'GVCO2 [umol/mol]\txCO2water_equ_dry [umol/mol]\txCO2water_SST_dry [umol/mol]\t',...
        'pCO2water_equ_wet [uatm]\tpCO2water_SST_wet [uatm]\tfCO2water_equ_wet [uatm]\tfCO2water_SST_wet [uatm]\t',...
        'fCO2rec [uatm]\tfCO2rec_src\tfCO2rec_flag'])) % header of SOCAT synthesis files from v2020
    error('Revise the columns headers and formatSpec');
end

% Read columns of data as text
% For more information, see the TEXTSCAN documentation.

if length(strsplit(char(colheaders{1}), '\t'))<=27 %synthesis file
    formatSpec = '%s%s%s%s%d%d%d%d%d%d%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%d%[^\n\r]';
    ifCO2rec=24;
    ifCO2rec_flag=26;
else % enhanced files and v2020 synthesis
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
%netcdf.putAtt(ncid, inputvar, 'long_name', 'quality flag');
netcdf.putAtt(ncid, inputvar, 'conventions', 'Copernicus Marine In Situ reference table 2');
netcdf.putAtt(ncid, inputvar, 'valid_min', int8(0));
netcdf.putAtt(ncid, inputvar, 'valid_max', int8(9));
netcdf.putAtt(ncid, inputvar, 'flag_values',  [int8(0)  int8(1)  int8(2)  int8(3)  int8(4)  int8(5)  int8(6)  int8(7)  int8(8)  int8(9)]);
netcdf.putAtt(ncid, inputvar, 'flag_meanings','no_qc_performed good_data probably_good_data bad_data_that_are_potentially_correctable bad_data value_changed value_below_detection nominal_value interpolated_value missing_value');
end
