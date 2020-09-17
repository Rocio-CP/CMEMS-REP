% SOCAT gridded file. According to Corentin, keep it as-is, with added
% Copernicus-specific global attributes
% Split by variable???
clc; clear all; close all
workrootdir=...
    '/Users/rpr061/Documents/localtestarea/CARBON-REP-112020/';
indir=[workrootdir,'original_files/'];
%outdir=[workrootdir,'2018/SOCATv6/workspace/'];
outdir=[workrootdir,'SOCAT_REP_GRIDDED_FIELDS/'];

gridfiles=dir([indir,'*gridded*.nc']);

%
% Load existing netcdfs
for f=1:length(gridfiles)
    
    system(['cp ',indir,gridfiles(f).name,' ',outdir,gridfiles(f).name(1:end-3),'_CMEMS_OG.nc'  ]);
    system(['chflags nouchg ',outdir,gridfiles(f).name(1:end-3),'_CMEMS_OG.nc']);
    
%        cmode = netcdf.getConstant('NETCDF4');
%cmode = bitor(cmode,netcdf.getConstant('CLASSIC_MODEL'));
%  
%     ncid=netcdf.open([outdir,gridfiles(f).name(1:end-3),'_CMEMS.nc'], ...
%        cmode);
%    
%        netcdf.close(ncid);

   
    ncid=netcdf.open([outdir,gridfiles(f).name(1:end-3),'_CMEMS_OG.nc'], ...
        'NC_WRITE');
    
    netcdf.reDef(ncid);
    globid=netcdf.getConstant('GLOBAL');


    
    netcdf.putAtt(ncid, globid, 'data_type', 'Copernicus Marine gridded data');
    netcdf.putAtt(ncid, globid, 'data_mode', 'D');
    netcdf.putAtt(ncid, globid, 'title', 'Global Ocean - Gridded In Situ reprocessed carbon observations - SOCATv2020');
    netcdf.putAtt(ncid, globid, 'references', 'http://marine.copernicus.eu, http://www.socat.info/');
    netcdf.putAtt(ncid, globid, 'naming_authority', 'Copernicus Marine in situ');
    netcdf.putAtt(ncid, globid, 'id', [gridfiles(f).name(1:end-3),'_CMEMS']);
    netcdf.putAtt(ncid, globid, 'institution', 'Pacific Marine Environmental Laboratory (PMEL), National Oceanic and Atmospheric Administration (NOAA)');
    netcdf.putAtt(ncid, globid, 'institution_edmo_code', '1440');
    
    netcdf.putAtt(ncid, globid, 'area', 'Global Ocean');
    netcdf.putAtt(ncid, globid, 'geospatial_lat_min', '-90');
    netcdf.putAtt(ncid, globid, 'geospatial_lat_max', '90');
    netcdf.putAtt(ncid, globid, 'geospatial_lon_min', '-180');
    netcdf.putAtt(ncid, globid, 'geospatial_lon_max', '180');
    netcdf.putAtt(ncid, globid, 'geospatial_vertical_min', '5');
    netcdf.putAtt(ncid, globid, 'geospatial_vertical_max', '5');
    netcdf.putAtt(ncid, globid, 'time_coverage_start', '1970-01-01T00:00:00Z');
    netcdf.putAtt(ncid, globid, 'time_coverage_end', '2019-02-25T07:57:00Z');
    netcdf.putAtt(ncid, globid, 'longitude_resolution', '1.0');
    netcdf.putAtt(ncid, globid, 'latitude_resolution', '1.0');
    if contains(gridfiles(f).name,'coast')
        netcdf.putAtt(ncid, globid, 'longitude_resolution', '0.25');
        netcdf.putAtt(ncid, globid, 'latitude_resolution', '0.25');
    end
    netcdf.putAtt(ncid, globid, 'cdm_data_type', 'grid');
    
    netcdf.putAtt(ncid, globid, 'format_version', '1.4');
    netcdf.putAtt(ncid, globid, 'Conventions', 'CF-1.6 Copernicus-InSituTAC-FormatManual-1.4 Copernicus-InSituTAC-SRD-1.41 Copernicus-InSituTAC-ParametersList-3.2.0');
    netcdf.putAtt(ncid, globid, 'netcdf_version', 'netCDF-4');
    
    netcdf.putAtt(ncid, globid, 'references', 'http://marine.copernicus.eu https://www.socat.info');
    netcdf.putAtt(ncid, globid, 'data_assembly_center', 'BERGEN');
    netcdf.putAtt(ncid, globid, 'update_interval', 'yearly');
    netcdf.putAtt(ncid, globid, 'citation', ['These data were collected and made freely available by the Copernicus project and the programs that contribute to it. ',...
        'The Surface Ocean CO2 Atlas (SOCAT) is an international effort, endorsed by the ',...
        'International Ocean Carbon Coordination Project (IOCCP), the Surface Ocean ',...
        'Lower Atmosphere Study (SOLAS) and the Integrated Marine Biosphere Research ',...
        '(IMBeR) program, to deliver a uniformly quality-controlled surface ocean CO2 ',...
        'database. The many researchers and funding agencies responsible for the collection ',...
        'of data and quality control are thanked for their contributions to SOCAT. ',...
        'Cite as Bakker et al. (2016), Bakker et al. (2020).']);
    netcdf.putAtt(ncid, globid, 'distribution_statement',['These data follow Copernicus standards; they are public ',...
        'and free of charge. User assumes all risk for use of data. User must display citation in any publication or ',...
        'product using data. User must contact PI prior to any commercial use of data.']);
    
    netcdf.putAtt(ncid, globid, 'doi', 'https://doi.org/10.5194/essd-8-383-2016 https://doi.org/10.25921/4xkx-ss49');

    netcdf.putAtt(ncid, globid, 'date_update', datestr(today, 'yyyy-mm-ddTHH:MM:SSZ'));
    netcdf.putAtt(ncid, globid, 'history', [datestr(now, 'yyyy-mm-ddTHH:MM:SSZ'),' : Creation']);
    
    netcdf.putAtt(ncid, globid, 'contact', 'post@bcdc.uib.no');
    netcdf.putAtt(ncid, globid, 'author', 'SOCAT and Copernicus data provider');
    netcdf.putAtt(ncid, globid, 'data_assembly_center', 'BERGEN');
    netcdf.putAtt(ncid, globid, 'update_interval', 'yearly');
    netcdf.putAtt(ncid, globid, 'SOCAT_history', 'PyFerret V7.52 (optimized)  5-Jun-20');
    netcdf.putAtt(ncid, globid, 'SOCAT_Notes', 'SOCAT gridded v2020 05-June-2020');
    netcdf.putAtt(ncid, globid, 'summary', 'Surface Ocean Carbon Atlas (SOCAT) Gridded (binned) SOCAT observations, with a spatial grid of 1x1 degree and yearly in time. The gridded fields are computed from the monthly 1-degree gridded data, which uses only SOCAT datasets with QC flags of A through D and SOCAT data points flagged with WOCE flag values of 2. This yearly data is computed using data from the start to the end of each year as described in the summary attribute of each variable.');
    netcdf.putAtt(ncid, globid, 'caution', 'NO INTERPOLATION WAS PERFORMED. SIGNIFICANT BIASES ARE PRESENT IN THESE GRIDDED RESULTS DUE TO THE ARBITRARY AND SPARSE LOCATIONS OF DATA VALUES IN BOTH SPACE AND TIME.');
    netcdf.putAtt(ncid, globid, 'distribution_statement','These data follow Copernicus standards; they are public and free of charge. User assumes all risk for use of data. User must display citation in any publication or product using data. User must contact PI prior to any commercial use of data.');
    
    netcdf.close(ncid);
    
    S=ncinfo([outdir,gridfiles(f).name(1:end-3),'_CMEMS_OG.nc']);
    file_fmt=S.Format
    S.Format='netcdf4';
    ncwriteschema([outdir,gridfiles(f).name(1:end-3),'_CMEMS.nc'],S);
    
end

%system(['cd /Users/rpr061/Dropbox/BCDC_Projects/CMEMS_INSTAC/Releases/current_FormatChecker/; for f in ', outdir,'*.nc; do ./control.csh $f >> ',outdir,'formatcheckoutSOCATgrid; done'])
