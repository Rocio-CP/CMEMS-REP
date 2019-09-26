clc; clear all; close all

%%% Create CMEMS netcdf files
workdir=['/Users/rpr061/Dropbox/BCDC_Projects/CMEMS_INSTAC/REP_Products/'...
    '2018/SOCATv6/workspace/'];
cd(workdir)
typegrid = 'coastal';
outputfile=['SOCATv6_gridded_',typegrid,'_Copernicus_test.nc'];

%% CHECK TITLES IN THE PIT DOCUMENT!! OR LOICs MAIL. OR PUM?
% Different attributes depending on climatology
switch typegrid
    case 'decade'
        inputfile = 'SOCATv6_tracks_gridded_decades.nc';
        timevarname = 'tdecade';
        varappend = '_decade';
        filetitle = 'Global Ocean - Surface Ocean CO₂ Atlas (SOCAT) v6 Decadal 1x1 degree gridded dataset';
        attsummary = 'The weighted cruise means for the months within each decade is averaged.'
        globalsummary = ['Surface Ocean Carbon Atlas (SOCAT) Gridded (binned) SOCAT observations, with a spatial grid of ',...
            '1x1 degree and yearly in time. The gridded fields are computed from the monthly 1-degree gridded data, ',...
            'which uses only SOCAT datasets with QC flags of A through D and SOCAT data points flagged with WOCE ',...
            'flag values of 2. This decadal data is computed using data from the start to the end of each decade as ',...
            'described in the summary attribute of each variable. The first decade is 1-Jan-1970 through 31-Dec-1979, ',...
            'the last decade is generally a partial decade. No interpolation was performed. Significant biases are present ',...
            'in these gridded results due to the arbitrary and sparse locations of data values in both space and time'];
        
    case 'year'
        inputfile = 'SOCATv6_tracks_gridded_yearly.nc';
        timevarname = 'tyear';
        varappend = '_year'
        filetitle = 'Global Ocean - Surface Ocean CO₂ Atlas (SOCAT) v6 Yearly 1x1 degree gridded data';
        attsummary = 'The weighted cruise means for the months within the year is averaged.';
        globalsummary = ['Surface Ocean Carbon Atlas (SOCAT) Gridded (binned) SOCAT observations, with a spatial grid of ',...
            '1x1 degree and yearly in time. The gridded fields are computed from the monthly 1-degree gridded data, ',...
            'which uses only SOCAT datasets with QC flags of A through D and SOCAT data points flagged with WOCE ',...
            'flag values of 2. This yearly data is computed using data from the start to the end of each year as ',...
            'described in the summary attribute of each variable. No interpolation was performed. Significant biases are present ',...
            'in these gridded results due to the arbitrary and sparse locations of data values in both space and time'];
        
    case 'month'
        inputfile = 'SOCATv6_tracks_gridded_monthly.nc';
        timevarname = 'tmnth';
        varappend = '';
        filetitle = 'Global Ocean - Surface Ocean CO₂ Atlas (SOCAT) v6 Monthly 1x1 degree gridded data';
        attsummary = ['Mean computed by calculating the arithmetic mean value for each cruise passing through the cell ',...
            'and then averaging these datasets.'];
        globalsummary = ['Surface Ocean Carbon Atlas (SOCAT) Gridded (binned) SOCAT observations, with a spatial ',...
            'grid of 1x1 degree and monthly in time. The gridded fields are computed using only SOCAT ',...
            'datasets with QC flags of A through D and SOCAT data points flagged with WOCE flag values of 2. ',...
            'No interpolation was performed. Significant biases are present in these gridded results ',...
            'due to the arbitrary and sparse locations of data values in both space and time'];
        
    case 'coastal'
        inputfile = 'SOCATv6_qrtrdeg_gridded_coast_monthly.nc';
        timevarname = 'tmnth';
        varappend = '';
        filetitle = 'Global Ocean - Surface Ocean CO₂ Atlas (SOCAT) v6 Coastal Monthly quarter-degree gridded data';
        attsummary = ['Mean computed by calculating the arithmetic mean value for each cruise passing through the cell ',...
            'and then averaging these cruises.'];
        globalsummary = ['Surface Ocean Carbon Atlas (SOCAT) Gridded (binned) SOCAT observations, masked for coastal data only, ',...
            'with a spatial grid of quarter-degree and monthly in time. The gridded fields are computed using only ',...
            'SOCAT cruises with QC flags of A through D and SOCAT data points flagged with WOCE flag values of 2. ',...
            'The grid is monthly in time and quarter-degree in longitude and latitude. A coastal mask is applied to ',...
            'match the 400-Km coastal region of the SOCAT database. No interpolation was performed. Significant biases are present ',...
            'in these gridded results due to the arbitrary and sparse locations of data values in both space and time'];
        
end

% Read info
latin=ncread([workdir,inputfile],'ylat');
lonin=ncread([workdir,inputfile],'xlon');
timeintemp=ncread([workdir,inputfile], timevarname);

if strcmp (typegrid, 'coastal');
    tempin=ncread([workdir,inputfile], ['coast_sst_ave_weighted', varappend]);
    salin=ncread([workdir,inputfile], ['coast_salinity_ave_weighted', varappend]);
    fco2in=ncread([workdir,inputfile], ['coast_fco2_ave_weighted', varappend]);
else
    tempin=ncread([workdir,inputfile], ['sst_ave_weighted', varappend]);
    salin=ncread([workdir,inputfile], ['salinity_ave_weighted', varappend]);
    fco2in=ncread([workdir,inputfile], ['fco2_ave_weighted', varappend]);
end

% Transform time value from days from 1900 to days from 1950
if strcmp (typegrid, 'decade');
    btw0050 = daysact('01-jan-1900', '01-jan-1950');
else
    btw0050 = daysact('01-jan-1970', '01-jan-1950');
end
timein = double(timeintemp - btw0050);
lengthtime = length(timein);

run CMEMS_SOCATv6_attributes.m