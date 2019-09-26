% GLODAP gridded file. According to Corentin, keep it as-is, with added
% Copernicus-specific global attributes
clc; clear all; close all
workrootdir=...
    '/Users/rpr061/Dropbox/BCDC_Projects/CMEMS_INSTAC/REP_Products/';
indir=[workrootdir,'04_2018/GLODAPv2/original/GLODAPv2.2016b_MappedClimatologies/'];
%outdir=[workrootdir,'2018/GLODAPv2/workspace/'];
outdir=['/Users/rpr061/Documents/localtestarea/CMEMS/'...
    'GLODAP_REP_GRIDDED_FIELDS/CLIMATOLOGY/'];
filestoadd={'temperature','salinity','oxygen','NO3','PO4','silicate',...
    'pHtsinsitutp','pHts25p0','TCO2','Talk'};

%variables (not dimension ones)


% Load existing netcdfs
for v=1:length(filestoadd)

   system(['cp ',indir,'GLODAPv2.2016b.',filestoadd{v},'.nc ',outdir]);
   system(['chflags nouchg ',outdir,'GLODAPv2.2016b.',filestoadd{v},'.nc']);

   varschangedim={filestoadd{v}, [filestoadd{v},'_error'],[filestoadd{v},'_relerr'],...
       'Input_mean','Input_std','Input_N'};
   filestoadd{v}
   S=ncinfo([indir,'GLODAPv2.2016b.',filestoadd{v},'.nc']);
   
   
   
   for vcd=1:length(varschangedim)
       varind=find(strcmp(varschangedim(vcd),{S.Variables.Name}));
       if strcmp(filestoadd{v},'Talk')
            varind=find(strcmp('TAlk',{S.Variables.Name}));
       end
       varschangedim{vcd}
       S.Variables(varind).Dimensions(4).Name='time';
       S.Variables(varind).Dimensions(4).Length=1;
       S.Variables(varind).Dimensions(4).Unlimited=false;
   end
   
   % Set all fillvalues; otherwise it'll create its own with 9.9e36
   
   for allv=1:length(S.Variables)
       S.Variables(allv).FillValue=-999;
   end
   
   % Global attributes

  S.Attributes(1).Name='data_type'; S.Attributes(1).Value='Copernicus Marine gridded data';
  S.Attributes(2).Name='format_version'; S.Attributes(2).Value='1.0';
  S.Attributes(3).Name='date_update'; S.Attributes(3).Value=datestr(today, 'yyyy-mm-ddTHH:MM:SSZ');
  S.Attributes(4).Name='institution'; S.Attributes(4).Value='Bjerknes Centre for Climate Research';
  S.Attributes(5).Name='institution_edmo_code'; S.Attributes(5).Value='1411';
  S.Attributes(6).Name='data_mode'; S.Attributes(6).Value='D';
  S.Attributes(7).Name='references'; S.Attributes(7).Value='http://marine.copernicus.eu, https://www.glodap.info/';
  S.Attributes(8).Name='Conventions'; S.Attributes(8).Value='CMEMS-INS-PUM-013-049-050';
  S.Attributes(9).Name='history'; S.Attributes(9).Value=[datestr(today, 'yyyy-mm-ddTHH:MM:SSZ'),' : Creation'];
  S.Attributes(10).Name='naming_authority'; S.Attributes(10).Value='Copernicus Marine';
  S.Attributes(11).Name='id'; S.Attributes(11).Value=['GLODAPv2.2016b.',filestoadd{v},'_CMEMS'];
  S.Attributes(12).Name='cdm_data_type'; S.Attributes(12).Value='grid';
  S.Attributes(13).Name='area'; S.Attributes(13).Value='Global Ocean';
  S.Attributes(14).Name='geospatial_lat_min'; S.Attributes(14).Value='-90';
  S.Attributes(15).Name='geospatial_lat_max'; S.Attributes(15).Value='90';
  S.Attributes(16).Name='geospatial_lon_min'; S.Attributes(16).Value='-180';
  S.Attributes(17).Name='geospatial_lon_max'; S.Attributes(17).Value='180';
  S.Attributes(18).Name='geospatial_vertical_min'; S.Attributes(18).Value='0';
  S.Attributes(19).Name='geospatial_vertical_max'; S.Attributes(19).Value='5500';
  S.Attributes(20).Name='time_coverage_start'; S.Attributes(20).Value='1972-01-01T00:00:00Z';
  S.Attributes(21).Name='time_coverage_end'; S.Attributes(21).Value='2013-12-31T23:59:59Z';
  S.Attributes(22).Name='longitude_resolution'; S.Attributes(22).Value='1';
  S.Attributes(23).Name='latitude_resolution'; S.Attributes(23).Value='1';
  S.Attributes(24).Name='contact'; S.Attributes(24).Value='post@bcdc.uib.no';
  S.Attributes(25).Name='author'; S.Attributes(25).Value='GLODAP and Copernicus data provider';
  S.Attributes(26).Name='data_assembly_center'; S.Attributes(26).Value='BERGEN';
  S.Attributes(27).Name='distribution_statement'; S.Attributes(27).Value='These data follow Copernicus standards; they are public and free of charge. User assumes all risk for use of data. User must display citation in any publication or product using data. User must contact PI prior to any commercial use of data.';
  S.Attributes(28).Name='citation'; S.Attributes(28).Value=['These data were collected and made freely available by the Copernicus project and the programs that contribute to it. ',...
      'Cite as: Lauvset, Siv K.; Key, Robert M.; Olsen, Are; van Heuven, Steven; ',...
    'Velo, Anton; Lin, Xiaohua; Schirnick, Carsten; Kozyr, Alex; ',...
    'Tanhua, Toste; Hoppema, Mario; Jutterstrom, Sara; Steinfeldt, Reiner; ',...
    'Jeansson, Emil; Ishii, Masao; Perez, Fiz F.; Suzuki, Toru; ',...
    'and Watelet, Sylvain: A new global interior ocean mapped climatology: ',...
    'the 1°x1° GLODAP version 2, Earth Syst. Sci. Data, 8, 325–340, 2016, ',...
    'doi:10.5194/essd-8-325-2016'];
  S.Attributes(29).Name='update_interval'; S.Attributes(29).Value='yearly';
  S.Attributes(30).Name='qc_manual'; S.Attributes(30).Value='QUID_CMEMS-IN_GLO_CARBON_REP_OBSERVATIONS_013_050';
  S.Attributes(31).Name='description'; S.Attributes(31).Value='1 X 1 global mapped field of temperature from the GLODAPv2 data product. Mapping is performed using the DIVA software (Troupin et al., 2012). Error fields are calculated using the clever poor mans error calculation method in DIVA (Beckers et al., 2014). The error fields represent the mapping error only, and does not include measurement or calculation uncertainties in the input data.';
  S.Attributes(32).Name='created'; S.Attributes(32).Value='Created by Siv K. Lauvset on 12-May-2016 18:47:05';
  S.Attributes(33).Name='title'; S.Attributes(33).Value='Global Ocean Data Analysis Project for Carbon (GLODAP)';

  
  % end of defining attributes / dimensions
  
  % Create a new netcdf from S "template"
  ncwriteschema([outdir,'GLODAPv2.2016b.',filestoadd{v},'_CMEMS.nc'],S);
  
  
  % Add values
  for nv = 1:length(S.Variables)
      ncwrite([outdir,'GLODAPv2.2016b.',filestoadd{v},'_CMEMS.nc'],...
          S.Variables(nv).Name, ...
          ncread([indir,'GLODAPv2.2016b.',filestoadd{v},'.nc'],S.Variables(nv).Name));
      
  end
     system(['rm ',outdir,'GLODAPv2.2016b.',filestoadd{v},'.nc ']);
  
end
%%

system(['cd /Users/rpr061/Dropbox/BCDC_Projects/CMEMS_INSTAC/REP_Products/current_FormatChecker/; for f in ', outdir,'*.nc; do ./control.csh $f >> ',outdir,'formatcheckoutGLODAPgrid; done'])

