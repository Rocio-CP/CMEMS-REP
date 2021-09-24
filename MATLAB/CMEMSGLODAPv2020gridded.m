    
clc; clear all; close all
%workrootdir=...
%    '/Users/rpr061/Dropbox/BCDC_Projects/CMEMS_INSTAC/REP_Products/';
indir=['/Users/rpr061/Documents/localtestarea/CARBON-REP-112020/glodadgriddednewformat/GLODAPv2.2016b_MappedClimatologies/'];
%outdir=[workrootdir,'2018/GLODAPv2/workspace/'];
outdir=['/Users/rpr061/Documents/localtestarea/CARBON-REP-112020/glodadgriddednewformat/'];
filestoadd={'temperature','salinity','oxygen','NO3','PO4','silicate',...
    'pHtsinsitutp','pHts25p0','TCO2','Talk'};

for v=1:length(filestoadd)

   system(['cp ',indir,'GLODAPv2.2016b.',filestoadd{v},'.nc ',outdir,'GLODAPv2.2016b.',filestoadd{v},'_CMEMS.nc']);
   system(['chflags nouchg ',outdir,'GLODAPv2.2016b.',filestoadd{v},'_CMEMS.nc']);


   
%% Make CF-compliant
ncwriteatt([outdir,'GLODAPv2.2016b.',filestoadd{v},'_CMEMS.nc'],'lon','standard_name','longitude');
ncwriteatt([outdir,'GLODAPv2.2016b.',filestoadd{v},'_CMEMS.nc'],'lon','units','degree_east');
ncwriteatt([outdir,'GLODAPv2.2016b.',filestoadd{v},'_CMEMS.nc'],'lat','standard_name','latitude');
ncwriteatt([outdir,'GLODAPv2.2016b.',filestoadd{v},'_CMEMS.nc'],'lat','units','degree_north');
ncwriteatt([outdir,'GLODAPv2.2016b.',filestoadd{v},'_CMEMS.nc'],'Depth','standard_name','depth');

% Fix units:
fileinfo=ncinfo([outdir,'GLODAPv2.2016b.',filestoadd{v},'_CMEMS.nc']);
allvars={fileinfo.Variables.Name};
for av=1:length(allvars);
    try ncreadatt([outdir,'GLODAPv2.2016b.',filestoadd{v},'_CMEMS.nc'],allvars{av},'units')
    catch ME
        continue
    end
    
    if isempty(ncreadatt([outdir,'GLODAPv2.2016b.',filestoadd{v},'_CMEMS.nc'],allvars{av},'units'))
continue
    elseif strcmp(ncreadatt([outdir,'GLODAPv2.2016b.',filestoadd{v},'_CMEMS.nc'],allvars{av},'units'),...
            'degrees celcius');
        ncwriteatt([outdir,'GLODAPv2.2016b.',filestoadd{v},'_CMEMS.nc'],allvars{av},'units','degree_C');
        
    elseif strcmp(ncreadatt([outdir,'GLODAPv2.2016b.',filestoadd{v},'_CMEMS.nc'],allvars{av},'units'),...
            'micro-mol kg-1');
        ncwriteatt([outdir,'GLODAPv2.2016b.',filestoadd{v},'_CMEMS.nc'],allvars{av},'units','umol kg-1');
    end
end


    %% Make Copernicus-compliants (global attributes
       ncid=netcdf.open([outdir,'GLODAPv2.2016b.',filestoadd{v},'_CMEMS.nc'], ...
        'NC_WRITE');
netcdf.reDef(ncid);
globid=netcdf.getConstant('GLOBAL');

    % Remove some global attributes (??

    netcdf.delAtt(ncid, globid,'Description');
    netcdf.delAtt(ncid, globid,'Created');
    netcdf.delAtt(ncid, globid,'Institution name');
    netcdf.delAtt(ncid, globid,'Contact information');
    netcdf.delAtt(ncid, globid,'Citation');
    netcdf.delAtt(ncid, globid,'Comment');

% INSTAC attributes    
    netcdf.putAtt(ncid, globid, 'data_type', 'Copernicus Marine gridded data');
    netcdf.putAtt(ncid, globid, 'data_mode', 'D');
    netcdf.putAtt(ncid, globid, 'title', 'Global Ocean - Gridded In Situ reprocessed carbon observations - GLODAPv2.2016b');
    netcdf.putAtt(ncid, globid, 'references', 'http://marine.copernicus.eu, http://www.glodap.info/');
    netcdf.putAtt(ncid, globid, 'naming_authority', 'Copernicus Marine in situ');
    netcdf.putAtt(ncid, globid, 'id', ['GLODAPv2.2016b.',filestoadd{v},'_CMEMS.nc']);
    netcdf.putAtt(ncid, globid, 'institution', 'Bjerknes Centre for Climate Research');
    netcdf.putAtt(ncid, globid, 'institution_edmo_code', '1411');
    
    netcdf.putAtt(ncid, globid, 'area', 'Global Ocean');
    netcdf.putAtt(ncid, globid, 'geospatial_lat_min', '-90');
    netcdf.putAtt(ncid, globid, 'geospatial_lat_max', '90');
    netcdf.putAtt(ncid, globid, 'geospatial_lon_min', '-180');
    netcdf.putAtt(ncid, globid, 'geospatial_lon_max', '180');
    netcdf.putAtt(ncid, globid, 'geospatial_vertical_min', '5');
    netcdf.putAtt(ncid, globid, 'geospatial_vertical_max', '5');
    netcdf.putAtt(ncid, globid, 'time_coverage_start', '1972-01-01T00:00:00Z');
    netcdf.putAtt(ncid, globid, 'time_coverage_end', '2013-12-31T23:59:59Z');
    netcdf.putAtt(ncid, globid, 'longitude_resolution', '1.0');
    netcdf.putAtt(ncid, globid, 'latitude_resolution', '1.0');

    netcdf.putAtt(ncid, globid, 'cdm_data_type', 'grid');
    
    netcdf.putAtt(ncid, globid, 'format_version', '1.4');
    netcdf.putAtt(ncid, globid, 'Conventions', 'CF-1.6 Copernicus-InSituTAC-FormatManual-1.4 Copernicus-InSituTAC-SRD-1.41 Copernicus-InSituTAC-ParametersList-3.2.0');
    netcdf.putAtt(ncid, globid, 'netcdf_version', 'netCDF-4');
    
    netcdf.putAtt(ncid, globid, 'references', 'http://marine.copernicus.eu https://www.glodap.info');
    netcdf.putAtt(ncid, globid, 'citation', ['These data were collected and made freely available by the Copernicus project and the programs that contribute to it. ',...
        'Cite as: Lauvset, Siv K.; Key, Robert M.; Olsen, Are; van Heuven, Steven; ',...
    'Velo, Anton; Lin, Xiaohua; Schirnick, Carsten; Kozyr, Alex; ',...
    'Tanhua, Toste; Hoppema, Mario; Jutterstrom, Sara; Steinfeldt, Reiner; ',...
    'Jeansson, Emil; Ishii, Masao; Perez, Fiz F.; Suzuki, Toru; ',...
    'and Watelet, Sylvain: A new global interior ocean mapped climatology: ',...
    'the 1°x1° GLODAP version 2, Earth Syst. Sci. Data, 8, 325–340, 2016, ',...
    'doi:10.5194/essd-8-325-2016.']);
    netcdf.putAtt(ncid, globid, 'distribution_statement',['These data follow Copernicus standards; they are public ',...
        'and free of charge. User assumes all risk for use of data. User must display citation in any publication or ',...
        'product using data. User must contact PI prior to any commercial use of data.']);
    
    netcdf.putAtt(ncid, globid, 'doi', 'https://doi.org/10.5194/10.5194/essd-2015-43');

    netcdf.putAtt(ncid, globid, 'date_update', datestr(today, 'yyyy-mm-ddTHH:MM:SSZ'));
    netcdf.putAtt(ncid, globid, 'history', ['2016-05-12T18:41:01: Creation by Siv K. Lauvset \n',...
    '2019-05-24T11:23:50Z: Created CMEMS INSTAC file \n',...
    datestr(now, 'yyyy-mm-ddTHH:MM:SSZ'),' : Updated to CMEMS INSTAC format']);
    
    netcdf.putAtt(ncid, globid, 'contact', 'post@bcdc.uib.no');
    netcdf.putAtt(ncid, globid, 'author', 'GLODAP and Copernicus data provider');
    netcdf.putAtt(ncid, globid, 'data_assembly_center', 'BERGEN');
    netcdf.putAtt(ncid, globid, 'update_interval', ' ');

    netcdf.putAtt(ncid, globid, 'summary', '1 X 1 global mapped field from the GLODAPv2 data product. Mapping is performed using the DIVA software (Troupin et al., 2012). Error fields are calculated using the clever poor mans error calculation method in DIVA (Beckers et al., 2014). The error fields represent the mapping error only, and does not include measurement or calculation uncertainties in the input data.');

            netcdf.close(ncid);

end