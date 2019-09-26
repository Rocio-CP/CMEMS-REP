globid = netcdf.getConstant('GLOBAL');
 
if datatype == 'observations'
netcdf.putAtt(ncid, globid, 'data_type', ['OceanSITES ',OStype,' data']);
end
netcdf.putAtt(ncid, globid, 'Conventions', 'CF-1.6 OceanSITES-Manual-1.2 Copernicus-InSituTAC-SRD-1.3 Copernicus-InSituTAC-ParametersList-3.1.0');
netcdf.putAtt(ncid, globid, 'format_version', '1.2');
 
netcdf.putAtt(ncid, globid, 'title', filetitle);
netcdf.putAtt(ncid, globid, 'id', outputfile(1:end-3));
netcdf.putAtt(ncid, globid, 'data_mode', 'D');
netcdf.putAtt(ncid, globid, 'format_version', '1.2')
netcdf.putAtt(ncid, globid, 'area', 'Global_Ocean');
 
netcdf.putAtt(ncid, globid, 'source', source);
netcdf.putAtt(ncid, globid, 'source_platform_category_code',sourceplatformcategorycode);
netcdf.putAtt(ncid, globid, 'platform_code',platformcode);
netcdf.putAtt(ncid, globid, 'platform_name',platformname);
netcdf.putAtt(ncid, globid, 'site_code', sitecode);
 
netcdf.putAtt(ncid, globid, 'citation', 'These data were collected and made freely available by the Copernicus project and the programs that contribute to it.');
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
netcdf.putAtt(ncid, globid, 'time_coverage_start', datestr(min(timeinog),'yyyy-mm-ddTHH:MM:SSZ'));
netcdf.putAtt(ncid, globid, 'time_coverage_end', datestr(max(timeinog),'yyyy-mm-ddTHH:MM:SSZ'));
netcdf.putAtt(ncid, globid, 'last_date_observation', datestr(timeinog(end),'yyyy-mm-ddTHH:MM:SSZ'));
netcdf.putAtt(ncid, globid, 'last_latitude_observation', latin(end));
netcdf.putAtt(ncid, globid, 'last_longitude_observation', lonin(end));
 
netcdf.putAtt(ncid, globid, 'cdm_data_type', cdmtype);
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
 
netcdf.putAtt(ncid, globid, 'qc_manual', qcmanual);
netcdf.putAtt(ncid, globid, 'quality_control_indicator', '6');
netcdf.putAtt(ncid, globid, 'quality_index', 'A');
 
netcdf.putAtt(ncid, globid, 'references', ['http://marine.copernicus.eu ',refwebsite]);
netcdf.putAtt(ncid, globid, 'comment', globcomment);
 

