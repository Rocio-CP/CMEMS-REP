clc;
clear all;
%clearvars -except SOCATv6
close all
workrootdir=...
    '/Users/rpr061/Documents/localtestarea/CARBON-REP-042020/';

% input folder
outdir=workrootdir;%['/Users/rpr061/Documents/localtestarea/CARBON-REP-042020/'];


%% Read SOCATv2019 synthesis file
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

%% Read SOCATv2109 enhanced files (A-E, WOCE 2-4)
if exist([workrootdir,'SOCATV2019all.mat'],'file')
    load([workrootdir,'SOCATV2019all.mat']);
else % read files
    allcruisefiles=dir([workrootdir,'SOCATv2019v7All_ABCDE_enhanced_datafiles/**/*.tsv']);
    for acf=1:length(allcruisefiles);
        cruisefile=[allcruisefiles(acf).folder,'/',allcruisefiles(acf).name];
        cruisedata=readSOCATenhancedfile(cruisefile);
        if acf==1; SOCATv2019all=cruisedata;
        else SOCATv2019all=CatStructFields(SOCATv2019all,cruisedata,1);
        end
        clear cruisedata
    end
    save([workrootdir,'SOCATv2019all.mat'], 'SOCATv2019all', '-v7.3');
end


%%


% From categorical to cellstr
% If done when creating SOCATv2019.mat, the variables is >2GB
SOCATv2019.Expocode=cellstr(SOCATv2019.Expocode);
SOCATv2019.ICEScode=cellstr(SOCATv2019.ICEScode);
SOCATv2019.PlatformName=cellstr(SOCATv2019.PlatformName);

% Remap some of the ICES codes; some ships have 2 when it should be 1
SOCATv2019.ICEScode(contains(SOCATv2019.ICEScode,'06MT'))={'06M3'};
SOCATv2019.ICEScode(contains(SOCATv2019.ICEScode,'06P1'))={'06PO'};
SOCATv2019.ICEScode(contains(SOCATv2019.ICEScode,'35MV'))={'35MF'};

% In SOCAT tsv NaN are NaN
unique_expocodes=unique(SOCATv2019.Expocode);

cmode = netcdf.getConstant('NETCDF4');
cmode = bitor(cmode,netcdf.getConstant('CLASSIC_MODEL'));
%%
for ec=1:length(unique_expocodes)
    
    currentexpocode=unique_expocodes{ec};
    currenticescode=currentexpocode(1:4);
    
    if strcmp(currenticescode,'06MT'); currenticescode='06M3';
    elseif  strcmp(currenticescode,'06P1'); currenticescode='06PO';
    elseif  strcmp(currenticescode,'35MV'); currenticescode='35MF';
    end
    
    %Which platform it is & set corresponding attributes
    % Research Vessel
    if any(ismember(ResearchVessels.SOCAT_code,...
            currenticescode))
        
        platformcategory='31';
        source='research vessel';
        platformcode=ResearchVessels.Callsign{...
            ismember(ResearchVessels.SOCAT_code,...
            currenticescode)}
        platformname=ResearchVessels.Name{...
            ismember(ResearchVessels.SOCAT_code,...
            currenticescode)}
        platformIMOnumber=ResearchVessels.IMO{...
            ismember(ResearchVessels.SOCAT_code,...
            currenticescode)};
        platformEDMOnumber=ResearchVessels.EDMO{...
            ismember(ResearchVessels.SOCAT_code,...
            currenticescode)}
        platformInstitute=ResearchVessels.Institute{...
            ismember(ResearchVessels.SOCAT_code,...
            currenticescode)}
        cdmtype='trajectory';
        ostype='OceanSITES trajectory data';
        
        if isempty(platformcode)
            platformcode=platformname(isstrprop(platformname, 'alphanum'));
        end
        outputfile = ['VESSEL/GL_TS_TS_',platformcode,'_SOCATv2019.nc'];
        if exist([outdir,outputfile], 'file')
            continue
        else
            expocodefilt=contains(SOCATv2019.ICEScode,currenticescode);
        end
        
        % Vessels of Opportunity
    elseif any(ismember(ShipsofOpportunity.SOCAT_code,...
            currenticescode))
        
        platformcategory='32';
        source='vessel of opportunity';
        platformcode=ShipsofOpportunity.Callsign{...
            ismember(ShipsofOpportunity.SOCAT_code,...
            currenticescode)}
        platformname=ShipsofOpportunity.Name{...
            ismember(ShipsofOpportunity.SOCAT_code,...
            currenticescode)}
        platformIMOnumber=ShipsofOpportunity.IMO{...
            ismember(ShipsofOpportunity.SOCAT_code,...
            currenticescode)};
        platformEDMOnumber=ShipsofOpportunity.EDMO{...
            ismember(ShipsofOpportunity.SOCAT_code,...
            currenticescode)}
        platformInstitute=ShipsofOpportunity.Institute{...
            ismember(ShipsofOpportunity.SOCAT_code,...
            currenticescode)}
        cdmtype='trajectory';
        ostype='OceanSITES trajectory data';
        
        if isempty(platformcode)
            platformcode=platformname(isstrprop(platformname, 'alphanum'));
        end
        
        outputfile = ['VESSEL/GL_TS_TS_',platformcode,'_SOCATv2019.nc'];
        if exist([outdir,outputfile], 'file')
            continue
        else
            expocodefilt=contains(SOCATv2019.ICEScode,currenticescode);
        end
        
        % Moorings
    elseif any(ismember(FixedMoorings.SOCAT_code,...
            currentexpocode))
        platformcategory='41';
        source='moored surface buoy';
        platformcode=FixedMoorings.WMO{...
            ismember(FixedMoorings.SOCAT_code,...
            currentexpocode)};
        platformname=FixedMoorings.Name{...
            ismember(FixedMoorings.SOCAT_code,...
            currentexpocode)}
        platformOScode=FixedMoorings.OceanSITES{...
            ismember(FixedMoorings.SOCAT_code,...
            currentexpocode)};
        platformEDMOnumber=FixedMoorings.EDMO{...
            ismember(FixedMoorings.SOCAT_code,...
            currentexpocode)}
        platformInstitute=FixedMoorings.Institute{...
            ismember(FixedMoorings.SOCAT_code,...
            currentexpocode)}
        cdmtype='station';
        ostype='OceanSITES time-series data';
        
        if isempty(platformcode)
            platformcode=platformname(isstrprop(platformname, 'alphanum'));
        end
        outputfile = ['MOORING/GL_TS_MO_',platformcode,'_SOCATv2019.nc'];
        if exist([outdir,outputfile], 'file')
            continue
        else
            expocodefilt=contains(SOCATv2019.PlatformName,platformname);
        end
        
        % Drifting buoys
    elseif any(ismember(DriftingBuoys.SOCAT_code,...
            currentexpocode))
        platformcategory='42';
        source='drifting surface float';
        platformcode=DriftingBuoys.WMO{...
            ismember(DriftingBuoys.SOCAT_code,...
            currenticescode)};
        platformname=DriftingBuoys.Name{...
            ismember(DriftingBuoys.SOCAT_code,...
            currenticescode)}
        platformOScode='';
        platformEDMOnumber=DriftingBuoys.EDMO{...
            ismember(DriftingBuoys.SOCAT_code,...
            currenticescode)}
        platformInstitute=DriftingBuoys.Institute{...
            ismember(DriftingBuoys.SOCAT_code,...
            currenticescode)}
        cdmtype='trajectory';
        ostype='OceanSITES trajectory data';
        
        if isempty(platformcode)
            platformcode=platformname(...
                isstrprop(platformname, 'alphanum')); end
        outputfile = ['GL_TS_DB_',platformcode,'_SOCATv2019.nc'];
        if exist([outdir,outputfile], 'file')
            continue
        else
            expocodefilt=contains(SOCATv2019.Expocode,currentexpocode);
        end
        
    else
        disp(['Can''t find the platform ', currentexpocode]);
    end
    
    % Extract info from SOCAT
    latin = SOCATv2019.latitudedecdegN(expocodefilt);
    lonin = SOCATv2019.longitudedecdegE(expocodefilt);
    
    timeinsocat = datenum(double([SOCATv2019.yr(expocodefilt),...
        SOCATv2019.mon(expocodefilt),SOCATv2019.day1(expocodefilt),...
        SOCATv2019.hh(expocodefilt),SOCATv2019.mm(expocodefilt),...
        SOCATv2019.ss(expocodefilt)]));
    timein = daysdif(datenum('1950-01-01T00:00:00Z','yyyy-mm-ddTHH:MM:SSZ'), ...
        timeinsocat);
    timeinqc=ones(size(timein));
    
    depthin=5.0*ones(size(timein));
    depthinqc=7*ones(size(timein)); % nominal depth
    depthindm=repmat('D',size(timein));
    
    possystin=repmat('U', size(timein));
    posinqc=ones(size(timein));
    
    fco2in = round(SOCATv2019.fCO2recuatm(expocodefilt)*1000.);
    fco2in(isnan(fco2in))=-9999;
    fco2inqc1 = SOCATv2019.fCO2rec_flag(expocodefilt);
    fco2inqc=fco2inqc1;
    fco2inqc(fco2inqc1==2)=1;
    fco2indm = repmat('D',size(fco2in));
    
    tempin = round(SOCATv2019.SSTdegC(expocodefilt)*1000.);
    tempin(isnan(tempin))=-9999;
    tempinqc = zeros(size(tempin));
    tempindm = repmat('D',size(tempin));
    
    salin = round(SOCATv2019.sal(expocodefilt)*1000.);
    salin(isnan(salin))=-9999;
    salinqc = zeros(size(salin));
    salindm = repmat('D',size(salin));
    
    % Create NetCDF file
    ncid = netcdf.create([outdir,outputfile], cmode);
    
    % Define dimensions and variables
    % Define dimensions
    latdim = netcdf.defDim(ncid, 'LATITUDE', length(latin));
    londim = netcdf.defDim(ncid, 'LONGITUDE', length(lonin));
    posdim = netcdf.defDim(ncid, 'POSITION', length(lonin));
    timedim = netcdf.defDim(ncid, 'TIME', netcdf.getConstant('NC_UNLIMITED'));
    depthdim = netcdf.defDim(ncid, 'DEPTH', 1);
    
    % Define dimension variables
    latvar = netcdf.defVar(ncid, 'LATITUDE', 'NC_FLOAT', latdim);
    lonvar = netcdf.defVar(ncid, 'LONGITUDE', 'NC_FLOAT', londim);
    timevar = netcdf.defVar(ncid, 'TIME', 'NC_DOUBLE', timedim);
    timevarqc = netcdf.defVar(ncid, 'TIME_QC', 'NC_BYTE', timedim);
    
    depthvar = netcdf.defVar(ncid, 'DEPH', 'NC_FLOAT', [depthdim, timedim]);
    depthvarqc = netcdf.defVar(ncid, 'DEPH_QC', 'NC_BYTE', [depthdim, timedim]);
    depthvardm = netcdf.defVar(ncid, 'DEPH_DM', 'NC_CHAR', [depthdim, timedim]);
    
    posvarqc = netcdf.defVar(ncid, 'POSITION_QC', 'NC_BYTE', posdim);
    possystvar = netcdf.defVar(ncid, 'POSITIONING_SYSTEM', 'NC_CHAR', posdim);
    
    % Define other variables
    tempvar = netcdf.defVar(ncid, 'TEMP', 'NC_INT', [depthdim, timedim]);
    salvar= netcdf.defVar(ncid, 'PSAL', 'NC_INT', [depthdim, timedim]);
    fco2var= netcdf.defVar(ncid, 'FCO2', 'NC_INT', [depthdim, timedim]);
    
    tempvarqc = netcdf.defVar(ncid, 'TEMP_QC', 'NC_BYTE', [depthdim, timedim]);
    salvarqc = netcdf.defVar(ncid, 'PSAL_QC', 'NC_BYTE', [depthdim, timedim]);
    fco2varqc = netcdf.defVar(ncid, 'FCO2_QC', 'NC_BYTE', [depthdim, timedim]);
    
    tempvardm = netcdf.defVar(ncid, 'TEMP_DM', 'NC_CHAR', [depthdim, timedim]);
    salvardm = netcdf.defVar(ncid, 'PSAL_DM', 'NC_CHAR', [depthdim, timedim]);
    fco2vardm = netcdf.defVar(ncid, 'FCO2_DM', 'NC_CHAR', [depthdim, timedim]);
    
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
    netcdf.putVar(ncid, tempvar, tempin);
    netcdf.putVar(ncid, salvar, salin);
    netcdf.putVar(ncid, fco2var, fco2in);
    
    netcdf.putVar(ncid, tempvarqc, tempinqc);
    netcdf.putVar(ncid, salvarqc, salinqc);
    netcdf.putVar(ncid, fco2varqc, fco2inqc);
    
    netcdf.putVar(ncid, tempvardm, tempindm);
    netcdf.putVar(ncid, salvardm, salindm);
    netcdf.putVar(ncid, fco2vardm, fco2indm);
    
    % Re-enter define mode (can be done before???)
    netcdf.reDef(ncid);
    
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
    netcdf.putAtt(ncid, depthvarqc, 'flag_values', int8([0  1  2  3  4  5  6  7  8  9]));
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
    
    % Add attributes to variables
    % Temperature
    netcdf.putAtt(ncid, tempvar, 'units','degrees_C');
    netcdf.putAtt(ncid, tempvar, 'standard_name','sea_water_temperature');
    netcdf.putAtt(ncid, tempvar, 'long_name','Sea temperature');
    netcdf.putAtt(ncid, tempvar, 'scale_factor',0.001);
    netcdf.putAtt(ncid, tempvar, 'add_offset',0);
    %  netcdf.putAtt(ncid, tempvar, 'summary', '')
    % netcdf.putAtt(ncid, tempvar, 'valid_min', -100.0);
    % netcdf.putAtt(ncid, tempvar, 'valid_max', 100.0);
    
    % Salinity
    netcdf.putAtt(ncid, salvar, 'units','0.001');
    netcdf.putAtt(ncid, salvar, 'standard_name','sea_water_practical_salinity');
    netcdf.putAtt(ncid, salvar, 'long_name','Practical salinity');
    netcdf.putAtt(ncid, salvar, 'scale_factor',0.001);
    netcdf.putAtt(ncid, salvar, 'add_offset',0);    %  netcdf.putAtt(ncid, salvar, 'summary', '')
    %   netcdf.putAtt(ncid, salvar, 'valid_min', 0.0);
    %   netcdf.putAtt(ncid, salvar, 'valid_max', 100.0);
    
    % fCO2
    netcdf.putAtt(ncid, fco2var, 'units','µatm');
    netcdf.putAtt(ncid, fco2var, 'standard_name','fugacity_of_CO2');
    netcdf.putAtt(ncid, fco2var, 'long_name','fugacity_of_CO2_in_sea_water');
    netcdf.putAtt(ncid, fco2var, 'scale_factor',0.001);
    netcdf.putAtt(ncid, fco2var, 'add_offset',0); % netcdf.putAtt(ncid, fco2var, 'summary', '')
    %  netcdf.putAtt(ncid, fco2var, 'valid_min', 0.0);
    %  netcdf.putAtt(ncid, fco2var, 'valid_max', 100000.0);
    
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
    netcdf.putAtt(ncid, depthvardm, 'long_name', 'method of data processing');
    netcdf.putAtt(ncid, depthvardm, 'conventions', 'OceanSITES reference table 5');
    netcdf.putAtt(ncid, depthvardm, 'flag_values', 'R, P, D, M');
    netcdf.putAtt(ncid, depthvardm, 'flag_meanings','real-time provisional delayed-mode mixed');
    
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
    
    
    % Add global attributes
    globid = netcdf.getConstant('GLOBAL');
    
    netcdf.putAtt(ncid, globid, 'data_type', ostype);
    netcdf.putAtt(ncid, globid, 'Conventions', 'CF-1.6 OceanSITES-Manual-1.2 Copernicus-InSituTAC-SRD-1.4 Copernicus-InSituTAC-ParametersList-3.1.0');
    netcdf.putAtt(ncid, globid, 'format_version', '1.2');
    
    netcdf.putAtt(ncid, globid, 'title', 'Global Ocean - Surface Ocean CO2 Atlas (SOCAT)');
    netcdf.putAtt(ncid, globid, 'id', outputfile(1:end-3));
    netcdf.putAtt(ncid, globid, 'data_mode', 'D');
    netcdf.putAtt(ncid, globid, 'area', 'Global_Ocean');
    
    netcdf.putAtt(ncid, globid, 'source', source);
    netcdf.putAtt(ncid, globid, 'source_platform_category_code',platformcategory);
    netcdf.putAtt(ncid, globid, 'platform_code',platformcode);
    netcdf.putAtt(ncid, globid, 'platform_name',platformname);
    netcdf.putAtt(ncid, globid, 'site_code','');
    % Platform type-specific attributes
    if strcmp(platformcategory,'31') | strcmp(platformcategory,'32')
        netcdf.putAtt(ncid, globid, 'imo_platform_number',platformIMOnumber);
        netcdf.putAtt(ncid, globid, 'ices_nodc_platform_code',currenticescode);
    elseif    strcmp(platformcode,'41')
        netcdf.putAtt(ncid, globid, 'wmo_platform_code',platformcode);
        netcdf.putAtt(ncid, globid, 'site_code',platformOScode);
    elseif   strcmp(platformcode,'43')
        netcdf.putAtt(ncid, globid, 'wmo_platform_code',platformcode);
    end
    
    netcdf.putAtt(ncid, globid, 'citation', ['These data were collected and made freely available by the Copernicus project and the programs that contribute to it. ',...
        'Cite as Bakker et al. (2016) and Bakker et al. (2018)']);
    netcdf.putAtt(ncid, globid, 'distribution_statement','These data follow Copernicus standards; they are public and free of charge. User assumes all risk for use of data. User must display citation in any publication or product using data. User must contact PI prior to any commercial use of data.');
    netcdf.putAtt(ncid, globid, 'naming_authority', 'Copernicus');
    
    netcdf.putAtt(ncid, globid, 'geospatial_lat_min', num2str(min(latin)));
    netcdf.putAtt(ncid, globid, 'geospatial_lat_max', num2str(max(latin)));
    netcdf.putAtt(ncid, globid, 'geospatial_lat_units', 'degree_north');
    netcdf.putAtt(ncid, globid, 'geospatial_lon_min', num2str(min(lonin)));
    netcdf.putAtt(ncid, globid, 'geospatial_lon_max', num2str(max(lonin)));
    netcdf.putAtt(ncid, globid, 'geospatial_lon_units', 'degree_east');
    netcdf.putAtt(ncid, globid, 'geospatial_vertical_min', '5.0');
    netcdf.putAtt(ncid, globid, 'geospatial_vertical_max', '5.0');
    netcdf.putAtt(ncid, globid, 'geospatial_vertical_positive', 'down');
    netcdf.putAtt(ncid, globid, 'geospatial_vertical_units', 'meter');
    netcdf.putAtt(ncid, globid, 'time_coverage_start', datestr(min(timeinsocat),'yyyy-mm-ddTHH:MM:SSZ'));
    netcdf.putAtt(ncid, globid, 'time_coverage_end', datestr(max(timeinsocat),'yyyy-mm-ddTHH:MM:SSZ'));
    netcdf.putAtt(ncid, globid, 'last_date_observation', datestr(timeinsocat(end),'yyyy-mm-ddTHH:MM:SSZ'));
    netcdf.putAtt(ncid, globid, 'last_latitude_observation', num2str(latin(end)));
    netcdf.putAtt(ncid, globid, 'last_longitude_observation', num2str(lonin(end)));
    
    netcdf.putAtt(ncid, globid, 'cdm_data_type', cdmtype);
    netcdf.putAtt(ncid, globid, 'keywords_vocabulary','SeaDataNet Parameter Discovery Vocabulary');
    
    netcdf.putAtt(ncid, globid, 'data_assembly_center', 'BERGEN');
    netcdf.putAtt(ncid, globid, 'institution', platformInstitute);
    netcdf.putAtt(ncid, globid, 'institution_edmo_code', platformEDMOnumber);
    netcdf.putAtt(ncid, globid, 'contact', 'post@bcdc.uib.no');
    netcdf.putAtt(ncid, globid, 'author', 'Bjerknes Climate Data Centre');
    
    netcdf.putAtt(ncid, globid, 'netcdf_version', '4');
    netcdf.putAtt(ncid, globid, 'update_interval', 'yearly');
    netcdf.putAtt(ncid, globid, 'date_update', datestr(today, 'yyyy-mm-ddTHH:MM:SSZ'));
    netcdf.putAtt(ncid, globid, 'history', [datestr(today, 'yyyy-mm-ddTHH:MM:SSZ'),' : Creation']);
    
    netcdf.putAtt(ncid, globid, 'qc_manual', 'SOCAT Quality Control Cookbook v3: https://www.socat.info/wp-content/uploads/2017/04/2015_SOCAT_QC_Cookbook_v3.pdf');
    netcdf.putAtt(ncid, globid, 'quality_control_indicator', '6');
    netcdf.putAtt(ncid, globid, 'quality_index', 'A');
    %    netcdf.putAtt(ncid, globid, 'SOCAT_flag', '');
    
    %    netcdf.putAtt(ncid, globid, 'summary', '');
    netcdf.putAtt(ncid, globid, 'references', 'http://marine.copernicus.eu https://www.socat.info');
    netcdf.putAtt(ncid, globid, 'comment', ['Further details can be found in Bakker et al. (2016). A multi-decade record of high-quality fCO2 data in version 3 of the ',...
        'Surface Ocean CO2 Atlas (SOCAT), Earth Syst. Sci. Data, 8, 383-413, https://doi.org/10.5194/essd-8-383-2016. ',...
        'and Bakker et al. (2018): Surface Ocean CO2 Atlas (SOCAT) V6. PANGAEA, https://doi.org/10.1594/PANGAEA.890974']);
    netcdf.putAtt(ncid, globid, 'SOCAT_contributor_name', ['Bakker, Dorothee C E; Lauvset, Siv K; Wanninkhof, Rik; Castaño-Primo, Rocío; ',...
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
    %end
    
    
    clearvars '-except' cmode ec workrootdir outdir SOCATv2019 ...
        unique_expocodes ResearchVessels ShipsofOpportunity ...
        FixedMoorings DriftingBuoys
    
end
system(['cd /Users/rpr061/Dropbox/BCDC_Projects/CMEMS_INSTAC/REP_Products/current_FormatChecker/; for f in ', outdir,'VESSEL/*.nc; do ./control.csh $f >> ',outdir,'VESSEL/formatcheckoutSOCAT; done'])
system(['cd /Users/rpr061/Dropbox/BCDC_Projects/CMEMS_INSTAC/REP_Products/current_FormatChecker/; for f in ', outdir,'MOORING/*.nc; do ./control.csh $f >> ',outdir,'MOORING/formatcheckoutSOCAT; done'])

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
        'PPPP [hPa]	Pequ [hPa]\tWOA_SSS\tNCEP_SLP [hPa]\tETOPO2_depth [m]\tdist_to_land [km]\t',...
        'GVCO2 [umol/mol]\tfCO2rec [uatm]\tfCO2rec_src\tfCO2rec_flag'])) && ... % header of synthesis file
        ~strcmp(char(colheaders{1}), sprintf(['Expocode\tversion\tSOCAT_DOI\tQC_Flag\t',...
        'yr\tmon\tday\thh\tmm\tss\tlongitude [dec.deg.E]\tlatitude [dec.deg.N]\t',...
        'sample_depth [m]\tsal\tSST [deg.C]\tTequ [deg.C]\t',...
        'PPPP [hPa]	Pequ [hPa]\tWOA_SSS\tNCEP_SLP [hPa]\tETOPO2_depth [m]\tdist_to_land [km]\t',...
        'GVCO2 [umol/mol]\txCO2water_equ_dry [umol/mol]\txCO2water_SST_dry [umol/mol]\t',...
        'pCO2water_equ_wet [uatm]\tpCO2water_SST_wet [uatm]\tfCO2water_equ_wet [uatm]\tfCO2water_SST_wet [uatm]\t',...
        'fCO2rec [uatm]\tfCO2rec_src\tfCO2rec_flag'])) % header of SOCAT-enhanced files
    error('Revise the columns headers and formatSpec');
end

% Read columns of data as text
% For more information, see the TEXTSCAN documentation.
if length(colheaders)<=27
    formatSpec = '%s%s%s%s%d%d%d%d%d%d%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%d%[^\n\r]';
    ifCO2rec=24;
    ifCO2rec_flag=26;
else
    formatSpec = '%s%s%s%s%d%d%d%d%d%d%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%d%[^\n\r]';
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
