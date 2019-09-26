
# Clear plots
if(!is.null(dev.list())) dev.off()
# Clean workspace
rm(list=ls())

library(ncdf4)

ncpath <- "/Users/rpr061/Dropbox/BCDC_projects/CMEMS_INSTAC/REP_Products/2018/SOCATv6/workspace/"
ncname <- "SOCATv6_tracks_gridded_decades_CMEMS_test"
ncfname <- paste(ncpath, ncname, ".nc", sep="")
setwd(ncpath)
# If output file exist, remove
if (file.exists(ncfname))
  file.remove(ncfname)

### ------ 
# LOAD AND READ ORIGINAL NETCDF FILE
ncinput <- nc_open("SOCATv6_tracks_gridded_decades.nc")

# Extract info
timein1900 <- ncvar_get(ncinput, "tdecade")
timeout1950 <-as.numeric(as.Date(timein1900,origin="1900-01-01") - as.Date("1950-01-01"))
latin <- ncvar_get(ncinput, "ylat")
lonin <- ncvar_get(ncinput, "xlon")

tempin <- ncvar_get(ncinput, "sst_ave_weighted_decade")
salin <- ncvar_get(ncinput, "salinity_ave_weighted_decade")
fco2in <- ncvar_get(ncinput, "fco2_ave_weighted_decade")

nc_close(ncinput)

### ----
# DEFINE AND CREATE CMEMS NETCDF
# Define dimensions
timedim <- ncdim_def(name="TIME",units="days since 1950-01-01T00:00:00Z",  vals=timeout1950, unlim=TRUE)
depthdim <- ncdim_def(name="DEPTH", create_dimvar=F, units='', vals=as.integer(1))#units="meters", vals=c(5.0)) # for SOCAT; more for GLODAP; different for mapped product
latdim <- ncdim_def(name="LATITUDE", units="degrees_north", vals=latin)
#latdim <- ncdim_def(name="LATITUDE", create_dimvar=F, units='', vals=as.integer(c(1:length(latin))))# units="degrees_north", vals=latin)
londim <- ncdim_def(name="LONGITUDE", units="degrees_east", vals=lonin)
#posdim <- ncdim_def(name="POSITION", units="")

#latvar <- ncvar_def(name="LATITUDE", units="degrees_north", dim=list(), prec="float" )
depthvar <- ncvar_def(name="DEPTH", units="meters", dim=list(), prec="float")

# Define variables
tempvar <- ncvar_def(name="TEMP", units="degree_Celsius", dim=list(londim,latdim, timedim), missval=-99999.0, prec="float")
salvar <- ncvar_def(name="PSAL", units="PSU", dim=list(londim,latdim, timedim), missval=-99999.0, prec="float")
pco2var <- ncvar_def(name="FCO2", units="microatmospheres", dim=list(londim,latdim,timedim), missval=-99999.0, prec="float")


# Create NetCDF files
vars <- list(depthvar, tempvar, salvar, pco2var)
ncout <- nc_create(ncfname, vars, force_v4=TRUE)

# Put variables values
#ncvar_put(ncout, latvar, latin)
ncvar_put(ncout, depthvar, c(5.0))


ncvar_put(ncout, tempvar, tempin)
ncvar_put(ncout, salvar, salin)
ncvar_put(ncout, pco2var, fco2in)

# Add attributes to dimensions
ncatt_put(ncout, "TIME", "axis", "T")
ncatt_put(ncout, "TIME", "standard_name", "time")
ncatt_put(ncout, "TIME", "long_name", "Time of measurement")
ncatt_put(ncout, "TIME", "valid_min", 0.0)

ncatt_put(ncout, "LATITUDE", "axis", "Y")
ncatt_put(ncout, "LATITUDE", "standard_name", "latitude")
ncatt_put(ncout, "LATITUDE", "long_name", "Latitude of each location")
ncatt_put(ncout, "LATITUDE", "valid_min", -90)
ncatt_put(ncout, "LATITUDE", "valid_max", 90)

ncatt_put(ncout, "LONGITUDE", "axis", "X")
ncatt_put(ncout, "LONGITUDE", "standard_name", "longitude")
ncatt_put(ncout, "LONGITUDE", "long_name", "Longitude of each location")
ncatt_put(ncout, "LONGITUDE", "valid_min", -180.)
ncatt_put(ncout, "LONGITUDE", "valid_max", 180.)

ncatt_put(ncout, "DEPTH", "axis", "Z")
ncatt_put(ncout, "DEPTH", "standard_name", "depth")
ncatt_put(ncout, "DEPTH", "long_name", "Depth of each measurement")
ncatt_put(ncout, "DEPTH", "valid_min", 0.0)
ncatt_put(ncout, "DEPTH", "valid_max", 12000.0)

# ADD ATTRIBUTES TO VARIABLES
ncatt_put(ncout, "TEMP", "standard_name","sea_water_temperature")
ncatt_put(ncout, "TEMP", "long_name","sea_water_temperature")
ncatt_put(ncout, "TEMP", "summary","The weighted cruise means for the months within each decade is averaged.")
ncatt_put(ncout, "TEMP", "valid_min", -100.)
ncatt_put(ncout, "TEMP", "valid_max", 100.)
ncatt_put(ncout, "TEMP", "coordinates","longitude latitude time")


ncatt_put(ncout, "PSAL", "standard_name","sea_water_salinity")
ncatt_put(ncout, "PSAL", "long_name","sea_water_salinity")
ncatt_put(ncout, "PSAL", "summary","The weighted cruise means for the months within each decade is averaged.")
ncatt_put(ncout, "PSAL", "valid_min", 0.)
ncatt_put(ncout, "PSAL", "valid_max", 100.)
ncatt_put(ncout, "PSAL", "coordinates","longitude latitude time")

ncatt_put(ncout, "FCO2", "standard_name","fugacity_of_CO2")
ncatt_put(ncout, "FCO2", "long_name","fugacity_of_CO2_in_sea_water")
ncatt_put(ncout, "FCO2", "summary","The weighted cruise means for the months within each decade is averaged.")
ncatt_put(ncout, "FCO2", "valid_min", 0.)
ncatt_put(ncout, "FCO2", "valid_max", 10000.)
ncatt_put(ncout, "FCO2", "coordinates","longitude latitude time")


# Add global attributes
ncatt_put(ncout,0, "data_mode", "D")
ncatt_put(ncout,0, "title", "SOCAT v6 Decadal 1x1 degree gridded dataset")
ncatt_put(ncout,0, "summary", "Surface Ocean Carbon Atlas (SOCAT) Gridded (binned) SOCAT observations, with a spatial grid of 
          1x1 degree and yearly in time. The gridded fields are computed from the monthly 1-degree gridded data, 
          which uses only SOCAT datasets with QC flags of A through D and SOCAT data points flagged with WOCE 
          flag values of 2. This decacal data is computed using data from the start to the end of each decade as 
          described in the summary attribute of each variable. The first decade is 1-Jan-1970 through 31-Dec-1979, 
          the last decade is generally a partial decade. No interpolation was performed. Significant biases are present 
          in these gridded results due to the arbitrar and sparse locations of data values in both space and time")
ncatt_put(ncout,0, "area", "Global")
ncatt_put(ncout,0, "geospatial_lat_min", "-89.5")
ncatt_put(ncout,0, "geospatial_lat_max", "89.5")
ncatt_put(ncout,0, "geospatial_lat_units", "degree_north")
ncatt_put(ncout,0, "geospatial_lon_min", "-179.5")
ncatt_put(ncout,0, "geospatial_lon_max", "179.5")
ncatt_put(ncout,0, "geospatial_lon_units", "degree_east")
ncatt_put(ncout,0, "geospatial_vertical_min", "5.0")
ncatt_put(ncout,0, "geospatial_vertical_max", "5.0")
ncatt_put(ncout,0, "geospatial_vertical_positive", "down")
ncatt_put(ncout,0, "geospatial_vertical_units", "meter")
ncatt_put(ncout,0, "time_coverage_start", "1970-01-01")
ncatt_put(ncout,0, "time_coverage_end", "2018-01-07")
ncatt_put(ncout,0, "cdm_data_type", "trajectory")
ncatt_put(ncout,0, "keywords_vocabulary","SeaDataNet Parameter Discovery Vocabulary")
ncatt_put(ncout,0, "references", "Bakker et al. (2016). A multi-decade record of high-quality fCO2 data in version 3 of the 
          Surface Ocean CO2 Atlas (SOCAT), Earth Syst. Sci. Data, 8, 383-413, https://doi.org/10.5194/essd-8-383-2016. 
          Bakker et al. (2018): Surface Ocean CO2 Atlas (SOCAT) V6. PANGAEA, https://doi.org/10.1594/PANGAEA.890974")
ncatt_put(ncout,0, "contributor_name", "Bakker, Dorothee C E; Lauvset, Siv K; Wanninkhof, Rik; Castaño-Primo, Rocío; 
          Currie, Kim I; Jones, Stephen D; Landa, Camilla S; Metzl, Nicolas; Nakaoka, Shin-Ichiro; Nojiri, Yukihiro; 
          Nonaka, Isao; O'Brien, Kevin M; Olsen, Are; Pfeil, Benjamin; Pierrot, Denis; Schuster, Ute; Smith, Karl; 
          Sullivan, Kevin; Sutton, Adrienne; Tilbrook, Bronte; Alin, Simone; Becker, Meike; Benoit-Cattin, Alice; 
          Bott, Randy; Bozec, Yann; Bozzano, Roberto; Burger, Eugene; Burgers, Tonya; Cai, Wei-Jun; Chen, Liqi; 
          Chierici, Melissa; Corredor, Jorge; Cosca, Catherine E; Cross, Jessica; Dandonneau, Yves; De Carlo, Eric Heinen; 
          Dietrich, Colin; Else, Brent; Emerson, Steven R; Farias, Laura; Fransson, Agneta; Garreaud, René D; Gkritzalis, Thanos; 
          Glockzin, Michael; González-Dávila, Melchor; Gregor, Luke; Hartman, Sue E; Hermes, Rudolf; Hoppema, Mario; Howden, Stephan; 
          Hunt, Christopher W; Hydes, David; Ibánhez, J Severino P; Kitidis, Vassilis; Körtzinger, Arne; Kozyr, Alexander;
          Kuwata, Akira; Lampitt, Richard Stephen; Lefèvre, Nathalie; Lo Monaco, Claire; Maenner, Stacy M; Manke, Ansley; 
          Manzello, Derek P; McGillis, Wade; Mickett, John; Monteiro, Pedro M S; Morell, Julio; Morrison, Ru; Mucci, Alfonso; 
          Munro, David R; Musielewicz, Sylvia; Negri, Ruben M; Newberger, Timothy; Newton, Jan; Noakes, Scott; O'Brien, Chris; 
          Ólafsdóttir, Sólveig Rósa; Ólafsson, Jón; Ono, Tsuneo; Osborne, John; Ouyang, Zhangxian; Padín, Xose Antonio; Papakyriakou, Tim N; 
          Plüddemann, Albert J; Rehder, Gregor; Sabine, Christopher L; Sakurai, Keizo; Salisbury, Joe; Santana-Casiano, Juana Magdalena; 
          Schlitzer, Reiner; Schneider, Bernd; Send, Uwe; Skjelvan, Ingunn; Steinhoff, Tobias; Sulpis, Olivier; Sutherland, Stewart C; 
          Sweeney, Colm; Tadokoro, Kazuaki; Takahashi, Taro; Telszewski, Maciej; Thomas, Helmuth; Tomlinson, Michael; Trull, Tom W; 
          Valdimarsson, Hedinn; van Heuven, Steven; Vandemark, Doug; Wada, Chisato; Wallace, Douglas WR; Watson, Andrew J; Weller, Robert A; Xu, Suqing")
ncatt_put(ncout,0, "data_assembly_center", "University of Bergen, Bjerknes Climate Data Centre")
ncatt_put(ncout,0, "netcdf_version", "4")
ncatt_put(ncout,0, "update_interval", "P1Y")
ncatt_put(ncout,0, "format_version", "1.2")
ncatt_put(ncout,0, "data_type", "OceanSITES trajectory data")
ncatt_put(ncout,0, "platform_code","NA")
ncatt_put(ncout,0, "site_code","NA")
ncatt_put(ncout,0, "date_update","2019-00-00T00:00:00Z")



nc_close(ncout)

#---------
# 
# # Global attributes
# ncatt_put(ncout,0, "site_code", "")
# ncatt_put(ncout,0, "platform_code", "")
# ncatt_put(ncout,0, "naming_authority", )
# ncatt_put(ncout,0, "id", )
# ncatt_put(ncout,0, "wmo_platform_code", )
# ncatt_put(ncout,0, "source", )
# ncatt_put(ncout,0, "network", )
# ncatt_put(ncout,0, "keywords", )
# ncatt_put(ncout,0, "comment", )
# ncatt_put(ncout,0, "time_coverage_duration", ) #Use ISO 8601 P1Y P3M P1Y3M4D
# ncatt_put(ncout,0, "data_type", )
# ncatt_put(ncout,0, "format_version", )
# ncatt_put(ncout,0, "Conventions", )
# ncatt_put(ncout,0, "publisher_name", )
# ncatt_put(ncout,0, "publisher_email", )
# ncatt_put(ncout,0, "publisher_url", )
# ncatt_put(ncout,0, "license", )
# ncatt_put(ncout,0, "citation", )
# ncatt_put(ncout,0, "acknowledgement", )
# ncatt_put(ncout,0, "date_created", )
# ncatt_put(ncout,0, "date_modified", )
# ncatt_put(ncout,0, "history", )
# ncatt_put(ncout,0, "processing_level", )
# ncatt_put(ncout,0, "QC_indicator", )
# ncatt_put(ncout,0, "contributor_role", )
# ncatt_put(ncout,0, "contributor_email", )
# 
# 
# # Variable properties
# ncatt_put(ncout, "TEMP", "long_name", )
# ncatt_put(ncout, "TEMP", "QC_indicator", )
# ncatt_put(ncout, "TEMP", "processing_level", )
# ncatt_put(ncout, "TEMP", "valid_min", )
# ncatt_put(ncout, "TEMP", "valid_max", )
# ncatt_put(ncout, "TEMP", "comment", )
# ncatt_put(ncout, "TEMP", "ancillary_variables", )
# ncatt_put(ncout, "TEMP", "uncertainty", )
# ncatt_put(ncout, "TEMP", "accuracy", )
# ncatt_put(ncout, "TEMP", "precision", )
# ncatt_put(ncout, "TEMP", "resolution", )
# ncatt_put(ncout, "TEMP", "reference_scale", )
# 
# ncatt_put(ncout, "TEMP_QC", "long_name", )
# ncatt_put(ncout, "TEMP_QC", "flag_values", )

