clc; clear all; close all
workrootdir=...
    '/Users/rpr061/Dropbox/BCDC_Projects/CMEMS_INSTAC/REP_Products/';

% input folder
load([workrootdir,'04_2018/GLODAPv2/workspace/GLODAPv2.mat']);
load([workrootdir,'04_2018/ReferenceData/PlatformCodes.mat'], 'GLODAPcruises');
%outdir=[workrootdir,'2018/GLODAPv2/workspace/GLODAPsplit_test/'];
outdir=['/Users/rpr061/Documents/localtestarea/CMEMS/'...
    'GLODAP_REP_OBSERVATIONS/VESSEL/'];

GLODAPcruises.Cruise=str2double(GLODAPcruises.Cruise);
% Change all -9999 with NaN
fn=fieldnames(GLODAPv2);
for nfn=1:length(fn)
    GLODAPv2.(fn{nfn})(GLODAPv2.(fn{nfn})==-9999)=NaN; end

% Map flag variables to OceanSITES: NaN -> 0; 0->8; 1->9; 2->1; 3->2; 4->4;
% 5->9; 6->8; 7->0; 8->0; 9->9
flaggedvars={'salinityf','oxygenf','nitratef',...
    'nitritef','phosphatef','silicatef','phtsinsitutpf','phts25p0f','tco2f',...
    'talkf','docf','donf','tdnf','chlaf'};

for nfv=1:length(flaggedvars)
    %  GLODAPv2.(flaggedvars{nfv})(isnan(GLODAPv2.(flaggedvars{nfv})))=0;
    GLODAPv2.(flaggedvars{nfv})(GLODAPv2.(flaggedvars{nfv})==0)=8;
    GLODAPv2.(flaggedvars{nfv})(GLODAPv2.(flaggedvars{nfv})==1)=9;
    GLODAPv2.(flaggedvars{nfv})(GLODAPv2.(flaggedvars{nfv})==2)=1;
    GLODAPv2.(flaggedvars{nfv})(GLODAPv2.(flaggedvars{nfv})==3)=2;
    GLODAPv2.(flaggedvars{nfv})(GLODAPv2.(flaggedvars{nfv})==4)=4;
    GLODAPv2.(flaggedvars{nfv})(GLODAPv2.(flaggedvars{nfv})==5)=9;
    GLODAPv2.(flaggedvars{nfv})(GLODAPv2.(flaggedvars{nfv})==6)=8;
    GLODAPv2.(flaggedvars{nfv})(GLODAPv2.(flaggedvars{nfv})==7)=0;
    GLODAPv2.(flaggedvars{nfv})(GLODAPv2.(flaggedvars{nfv})==8)=0;
    GLODAPv2.(flaggedvars{nfv})(GLODAPv2.(flaggedvars{nfv})==9)=9;
    
    GLODAPv2.(flaggedvars{nfv})=int8(GLODAPv2.(flaggedvars{nfv}));
end

GLODAPv2.bottomdepthf=zeros(size(GLODAPv2.salinityf));
GLODAPv2.temperaturef=ones(size(GLODAPv2.salinityf));
GLODAPv2.doc1f=GLODAPv2.docf;

unique_cruise=unique(GLODAPv2.cruise);

%% Variables to be included in Copernicus

invars={'bottomdepth','temperature','salinity','oxygen','nitrate',...
    'nitrite','phosphate','silicate','phtsinsitutp','phts25p0','tco2',...
    'talk','doc1','don','tdn','chla'};
vars={'bath','temp','sal','oxy','nitra',...
    'nitri','phos','sili','ph','ph25','dic',...
    'alk','doc','don','tdn','chla'};
osvars={'BATH','TEMP','PSAL','DOX2','NTAW',...
    'NTIW','PHOW','SLCW','PHPH','PH25','TICW',...
    'ALKW','CORG','NODW','NT1D','CPHL'};
osunits={'meters','degrees_C','0.001','µmol kg-1','µmol kg-1',...
    'µmol kg-1','µmol kg-1','µmol kg-1','1','1','µmol kg-1',...
    'µmol kg-1','µmol kg-1','µmol kg-1','µmol kg-1','mg m-3'};
oslongname={'Bottom depth','Sea temperature','Practical salinity',...
    'Dissolved oxygen','Nitrate (NO3-N)',...
    'Nitrite (NO2-N)','Phosphate (PO-P)','Silicate (SiO4-Si)',...
    'Ph','Ph at 25 degrees and 0 dbar','Dissolved inorganic carbon',...
    'Total alkalinity','Dissolved organic carbon','Dissolved organic nitrogen',...
    'Total dissolved nitrogen','Chlorophyll-a'};
osstandardname={'bottom_depth','sea_water_temperature',...
    'sea_water_practical_salinity','moles_of_oxygen_per_unit_mass_in_sea_water',...
    'moles_of_nitrate_per_unit_mass_in_sea_water','moles_of_nitrite_per_unit_mass_in_sea_water',...
    'moles_of_phosphate_per_unit_mass_in_sea_water','moles_of_silicate_per_unit_mass_in_sea_water',...
    'sea_water_ph_reported_on_total_scale','sea_water_ph_reported_on_total_scale_at_25_degrees_and_0_dbar',...
    'moles_of_dissolved_inorganic_carbon_per_unit_mass_in_sea_water','sea_water_alkalinity_per_unit_mass',...
    'moles_of_dissolved_organic_carbon_per_unit_mass_in_sea_water','moles_of_dissolved_organic_nitrogen_per_unit_mass_in_sea_water',...
    'moles_of_dissolved_total_nitrogen_per_unit_mass_in_sea_water','mass_concentration_of_chlorophyll_a_in_sea_water'};

cmode = netcdf.getConstant('NETCDF4');
cmode = bitor(cmode,netcdf.getConstant('CLASSIC_MODEL'));




%%
% Either fix the issue of missing cruises or add the multiple ship cruises
for ec=[unique_cruise]'
    currentcruise=ec;
    currenticescode=GLODAPcruises.ICEScode{GLODAPcruises.Cruise==currentcruise};
    
    if strcmp(currenticescode,'1891') | strcmp(currenticescode,'ZZIC')
        platformcategory='62';
        source='aeroplane';
    elseif strcmp(currenticescode,'33PY')
        platformcategory='21';
        source='propelled manned submersible';
    else
        platformcategory= '31';
        source='research vessel';
        
    end
    ostype='OceanSITES vertical profile';
    cdmtype='profile';
    platformcode=GLODAPcruises.Callsign{GLODAPcruises.Cruise==currentcruise};
    platformname=GLODAPcruises.Name{GLODAPcruises.Cruise==currentcruise}
    platformIMOnumber=GLODAPcruises.IMO{GLODAPcruises.Cruise==currentcruise};
    
    if isempty(platformcode)
        platformcode=platformname(isstrprop(platformname, 'alphanum'));
    end
    outputfile = ['GL_PR_BO_',platformcode,'.nc'];
    if exist([outdir,outputfile], 'file')
        continue
    else
        whichcruises=GLODAPcruises.Cruise(...
            contains(GLODAPcruises.ICEScode,currenticescode));
        cruisefilter=ismember(GLODAPv2.cruise,whichcruises);
    end
    
    %% Extract the most frequent EDMO code (un poco churro; mejorar)
    [unique_EDMO, ~, EDMO_map]=unique(GLODAPcruises.EDMO(whichcruises));
    platformEDMO=unique_EDMO{mode(EDMO_map)}
    platforminstitutes=GLODAPcruises.Institute(whichcruises);
    platformInstitute=platforminstitutes{mode(EDMO_map)}
    
    %% Input variables
    depthinglodap=GLODAPv2.depth(cruisefilter);
    [depthin,inddepthin]=unique(depthinglodap);
    
    timeinglodapMATLAB = datenum(([GLODAPv2.year1(cruisefilter),...
        GLODAPv2.month1(cruisefilter),GLODAPv2.day1(cruisefilter),...
        GLODAPv2.hour1(cruisefilter),GLODAPv2.minute1(cruisefilter),...
        zeros(size(GLODAPv2.minute1(cruisefilter)))]));
    timeinglodap = daysdif(datenum('1950-01-01T00:00:00Z','yyyy-mm-ddTHH:MM:SSZ'), ...
        timeinglodapMATLAB);
    
    [timein,indtimein]=unique(timeinglodap);
    
    depthinmatrix=nan(length(depthin),length(timein));
    depthinqcmatrix=nan(length(depthin),length(timein));
    %depthindmmatrix=nan(length(depthin),length(timein));
    
    
    for oo=1:length(depthinglodap)
        depthind=find(depthin==depthinglodap(oo));
        timeind=find(timein==timeinglodap(oo));
        
        depthinmatrix(depthind,timeind)=nanmean(depthinglodap(oo));
        depthinqcmatrix(depthind,timeind)=1;
        depthindmmatrix(depthind,timeind)='D';
    end
    
    latinglodap = GLODAPv2.latitude(cruisefilter);
    latin=latinglodap(indtimein);
    loninglodap=GLODAPv2.longitude(cruisefilter);
    lonin=loninglodap(indtimein);
    
    for vv=1:length(vars)
        % Variables
        eval([vars{vv},'inglodap=round(GLODAPv2.',invars{vv},'(cruisefilter)*1000.);']);
        eval([vars{vv},'finglodap=GLODAPv2.',invars{vv},'f(cruisefilter);']);
        %       eval([vars{vv},'inglodap(',vars{vv},'inglodap==-9999)=NaN;']);
        
        % Create depthxtime variables
        eval([vars{vv},'in=nan(length(depthin),length(',...
            'timein));']);
        eval([vars{vv},'inqc=nan(length(depthin),length(',...
            'timein));']);
        
        for oo=1:length(depthinglodap)
            depthind=find(depthin==depthinglodap(oo));
            timeind=find(timein==timeinglodap(oo));
            
            eval([vars{vv},'in(',num2str(depthind),...
                ',',num2str(timeind),')=round(nanmean([',...
                vars{vv},'in(',num2str(depthind),...
                ',',num2str(timeind),'),',...
                vars{vv},'inglodap(oo)]));']);
            
            % If value is average of >1 bottle, flag is 8 (interpolated)
            if  eval(['isnan(',vars{vv},'inqc(',num2str(depthind),...
                    ',',num2str(timeind),'));'])
                
                eval([vars{vv},'inqc(',num2str(depthind),...
                    ',',num2str(timeind),')=',vars{vv},'finglodap(oo);']);
            else
                eval([vars{vv},'inqc(',num2str(depthind),...
                    ',',num2str(timeind),')=8;']);
            end
            
            % eval([vars{vv},'inqc(',num2str(depthind),...
            %     ',',num2str(timeind),')=0;']);
            
            eval([vars{vv},'indm(',num2str(depthind),...
                ',',num2str(timeind),')=''D'';']);
        end
        eval([vars{vv},'in(isnan(',vars{vv},'in))=-9999;']);
        eval([vars{vv},'inqc(isnan(',vars{vv},'inqc))=-9999;']);
        eval([vars{vv},'inqc=int8(',vars{vv},'inqc);']);
        
        
    end
    
    %    for vv=1:length(vars)
    %        eval([vars{vv},'in=round(',vars{vv},'in*1000.);']);
    %    end
    
    depthinqc=ones(size(depthin));
    depthindm=repmat('D',size(depthin));
    timeinqc=ones(size(timein));
    
    possystin=repmat('U', size(timein));
    posinqc=ones(size(timein));
    
    %% Generate NetCDF
    % Create NetCDF file
    ncid = netcdf.create([outdir,outputfile], cmode);
    
    % Define dimensions and variables
    % Define dimensions
    latdim = netcdf.defDim(ncid, 'LATITUDE', length(latin));
    londim = netcdf.defDim(ncid, 'LONGITUDE', length(lonin));
    posdim = netcdf.defDim(ncid, 'POSITION', length(lonin));
    timedim = netcdf.defDim(ncid, 'TIME', netcdf.getConstant('NC_UNLIMITED'));
    depthdim = netcdf.defDim(ncid, 'DEPTH', length(depthin));
    
    % Define dimension-ish variables
    latvar = netcdf.defVar(ncid, 'LATITUDE', 'NC_FLOAT', latdim);
    lonvar = netcdf.defVar(ncid, 'LONGITUDE', 'NC_FLOAT', londim);
    timevar = netcdf.defVar(ncid, 'TIME', 'NC_DOUBLE', timedim);
    depthvar = netcdf.defVar(ncid, 'DEPH', 'NC_FLOAT', [depthdim,timedim]);
    possystvar = netcdf.defVar(ncid, 'POSITIONING_SYSTEM', 'NC_CHAR', posdim);
    
    %posvar = netcdf.defVar(ncid, 'POSITION', 'NC_FLOAT', [timedim]);
    % Dimension variables QC DM
    depthvarqc = netcdf.defVar(ncid, 'DEPH_QC', 'NC_BYTE', [depthdim,timedim]);
    depthvardm = netcdf.defVar(ncid, 'DEPH_DM', 'NC_CHAR', [depthdim,timedim]);
    posvarqc = netcdf.defVar(ncid, 'POSITION_QC', 'NC_BYTE', posdim);
    timevarqc = netcdf.defVar(ncid, 'TIME_QC', 'NC_BYTE', timedim);
    
    % Define other variables
    for vv=1:length(vars)
        eval([vars{vv},'var=netcdf.defVar(ncid,''',osvars{vv},...
            ''',''NC_INT'',[depthdim,timedim]);']);
        eval([vars{vv},'varqc=netcdf.defVar(ncid,''',osvars{vv},...
            '_QC'',''NC_BYTE'',[depthdim,timedim]);']);
        eval([vars{vv},'vardm=netcdf.defVar(ncid,''',osvars{vv},...
            '_DM'',''NC_CHAR'',[depthdim,timedim]);']);
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
    for vv=1:length(vars)
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
    netcdf.putVar(ncid, depthvar, depthinmatrix);
    netcdf.putVar(ncid, depthvarqc, depthinqcmatrix);
    netcdf.putVar(ncid, depthvardm, depthindmmatrix);
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
    for vv=1:length(vars)
        % Variables
        eval(['netcdf.putAtt(ncid,',vars{vv},'var,''units'',''',osunits{vv},''');']);
        eval(['netcdf.putAtt(ncid,',vars{vv},'var,''standard_name'',''',osstandardname{vv},''');']);
        eval(['netcdf.putAtt(ncid,',vars{vv},'var,''long_name'',''',oslongname{vv},''');']);
        eval(['netcdf.putAtt(ncid,',vars{vv},'var,''scale_factor'',0.001);']);
        eval(['netcdf.putAtt(ncid,',vars{vv},'var,''add_offset'',0);']);
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
    
    netcdf.putAtt(ncid, globid, 'data_type', ostype);
    netcdf.putAtt(ncid, globid, 'Conventions', 'CF-1.6 OceanSITES-Manual-1.2 Copernicus-InSituTAC-SRD-1.4 Copernicus-InSituTAC-ParametersList-3.1.0');
    netcdf.putAtt(ncid, globid, 'format_version', '1.2');
    
    netcdf.putAtt(ncid, globid, 'title', 'Global Ocean Data Analysis Project for Carbon (GLODAP)');
    netcdf.putAtt(ncid, globid, 'id', outputfile(1:end-3));
    netcdf.putAtt(ncid, globid, 'data_mode', 'D');
    netcdf.putAtt(ncid, globid, 'format_version', '1.2')
    netcdf.putAtt(ncid, globid, 'area', 'Global_Ocean');
    
    netcdf.putAtt(ncid, globid, 'source', source);
    netcdf.putAtt(ncid, globid, 'source_platform_category_code',platformcategory);
    netcdf.putAtt(ncid, globid, 'platform_code',platformcode);
    netcdf.putAtt(ncid, globid, 'platform_name',platformname);
    netcdf.putAtt(ncid, globid, 'ices_nodc_platform_code',currenticescode);
    netcdf.putAtt(ncid, globid, 'site_code','');
    
    if  strcmp(platformcategory,'31')
        netcdf.putAtt(ncid, globid, 'imo_platform_number',platformIMOnumber);
    end
    
    netcdf.putAtt(ncid, globid, 'citation', ['These data were collected and made freely available by the Copernicus project and the programs that contribute to it. ',...
        'Cite as Olsen, A.; Key, R. M.; van Heuven, S.; Lauvset, S. K.; Velo,  A.; Lin, X.; Schirnick, C.; ',...
        'Kozyr, A.; Tanhua, T.; Hoppema, M.; Jutterström, S.; Steinfeldt, R.; Jeansson, E.; Ishii, M.; ',...
        'Pérez, F. F.; and T. Suzuki, T. (2016).',...
        'The Global Ocean Data Analysis Project version 2 (GLODAPv2) – an internally consistent data product for the world ocean, ',...
        'Earth Syst. Sci. Data, 8, 297-323, https://doi.org/10.5194/essd-8-297-2016']);
    netcdf.putAtt(ncid, globid, 'distribution_statement','These data follow Copernicus standards; they are public and free of charge. User assumes all risk for use of data. User must display citation in any publication or product using data. User must contact PI prior to any commercial use of data.');
    netcdf.putAtt(ncid, globid, 'naming_authority', 'Copernicus');
    
    netcdf.putAtt(ncid, globid, 'geospatial_lat_min', num2str(min(latin)));
    netcdf.putAtt(ncid, globid, 'geospatial_lat_max', num2str(max(latin)));
    netcdf.putAtt(ncid, globid, 'geospatial_lat_units', 'degree_north');
    netcdf.putAtt(ncid, globid, 'geospatial_lon_min', num2str(min(lonin)));
    netcdf.putAtt(ncid, globid, 'geospatial_lon_max', num2str(max(lonin)));
    netcdf.putAtt(ncid, globid, 'geospatial_lon_units', 'degree_east');
    netcdf.putAtt(ncid, globid, 'geospatial_vertical_min', num2str(min(depthin)));
    netcdf.putAtt(ncid, globid, 'geospatial_vertical_max', num2str(max(depthin)));
    netcdf.putAtt(ncid, globid, 'geospatial_vertical_positive', 'down');
    netcdf.putAtt(ncid, globid, 'geospatial_vertical_units', 'meter');
    netcdf.putAtt(ncid, globid, 'time_coverage_start', datestr(min(timeinglodap),'yyyy-mm-ddTHH:MM:SSZ'));
    netcdf.putAtt(ncid, globid, 'time_coverage_end', datestr(max(timeinglodap),'yyyy-mm-ddTHH:MM:SSZ'));
    netcdf.putAtt(ncid, globid, 'last_date_observation', datestr(timeinglodap(end),'yyyy-mm-ddTHH:MM:SSZ'));
    netcdf.putAtt(ncid, globid, 'last_latitude_observation', num2str(latin(end)));
    netcdf.putAtt(ncid, globid, 'last_longitude_observation', num2str(lonin(end)));
    
    netcdf.putAtt(ncid, globid, 'cdm_data_type', cdmtype);
    netcdf.putAtt(ncid, globid, 'keywords_vocabulary','SeaDataNet Parameter Discovery Vocabulary');
    
    netcdf.putAtt(ncid, globid, 'data_assembly_center', 'BERGEN');
    netcdf.putAtt(ncid, globid, 'institution', platformInstitute);
    netcdf.putAtt(ncid, globid, 'institution_edmo_code', platformEDMO);
    netcdf.putAtt(ncid, globid, 'contact', 'post@bcdc.uib.no');
    netcdf.putAtt(ncid, globid, 'author', 'Bjerknes Climate Data Centre');
    
    netcdf.putAtt(ncid, globid, 'netcdf_version', '4');
    netcdf.putAtt(ncid, globid, 'update_interval', 'yearly');
    netcdf.putAtt(ncid, globid, 'date_update', datestr(today, 'yyyy-mm-ddTHH:MM:SSZ'));
    netcdf.putAtt(ncid, globid, 'history', [datestr(today, 'yyyy-mm-ddTHH:MM:SSZ'),' : Creation']);
    
    netcdf.putAtt(ncid, globid, 'qc_manual', 'Olsen et. al (2016), ESSD');
    netcdf.putAtt(ncid, globid, 'quality_control_indicator', '6');
    netcdf.putAtt(ncid, globid, 'quality_index', 'A');
    
    netcdf.putAtt(ncid, globid, 'references', 'http://marine.copernicus.eu https://www.glodap.info');
    if ec==720 | ec==721
        netcdf.putAtt(ncid, globid, 'cruise_comment', 'Data gathered by RV Bjarni Saemundsson and RV Arni Fridrikkson (Callsigns TFNA/TFJA; IMO numbers 9192404/6710358)');
    elseif ec==722
        netcdf.putAtt(ncid, globid, 'cruise_comment', 'Data gathered by RV Belgica, RV Charles Darwin and RV Meteor');
    end
    
    % Close nc file
    netcdf.close(ncid);
    
    
    clearvars '-except' cmode ec workrootdir outdir GLODAPcruises GLODAPv2 ...
        unique_cruise invars vars osvars osunits oslongname osstandardname
    
end

system(['cd /Users/rpr061/Dropbox/BCDC_Projects/CMEMS_INSTAC/REP_Products/current_FormatChecker/; for f in ', outdir,'*.nc; do ./control.csh $f >> ',outdir,'formatcheckoutGLODAP; done'])






