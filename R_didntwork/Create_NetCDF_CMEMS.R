library(ncdf4)

ncpath <- "/Users/rpr061/Dropbox/BCDC_projects/CMEMS_INSTAC/REP_Products/2018/SOCATv6/"
ncname <- "SOCATv6_synthesis_test"
ncfname <- paste(ncpath, ncname, ".nc", sep="")

# Read SOCAT synthesis file
#skiplines <- as.integer(system("grep -n -m1 'Expocode\tversion\tSOCAT_DOI' SOCATv6.tsv | cut -f 1 -d :", intern=TRUE))
# R and dashes don't mix well. Create a SOCAT file without the header lines.
socatsynth <- read_delim("SOCATv6_woheader.tsv","\t",escape_double = FALSE,
  col_types = cols(day = col_character(), mm = col_character(), yr = col_character(), fCO2rec_flag = col_integer()),
  trim_ws = TRUE)

# Change time unit to days since 1950-01-01
socatsynth$OSdate <- as.numeric(difftime(ISOdatetime(socatsynth$yr,socatsynth$mon, socatsynth$day,
                                          socatsynth$hh,socatsynth$mm,socatsynth$ss,tz="UTC"), 
                              ISOdatetime(1950,01,01,00,00,00,tz="UTC")))
# Define dimensions
timedim <- ncdim_def(name="TIME",units="days since 1950-01-01T00:00:00Z", vals=socatsynth$OSdate, unlim=TRUE)
depthdim <- ncdim_def(name="DEPTH", units="meters", vals=c(5)) # for SOCAT; more for GLODAP; different for mapped product
latdim <- ncdim_def(name="LATITUDE", units="degrees_north", vals=socatsynth$`latitude [dec.deg.N]`, unlim=TRUE)
londim <- ncdim_def(name="LONGITUDE", units="degrees_east", vals=socatsynth$`longitude [dec.deg.E]`, unlim=TRUE)

# Define variables
tempvar <- ncvar_def(name="TEMP", units="degree_Celsius", dim=timedim)
salvar <- ncvar_def(name="PSAL", units="PSU", dim=timedim)
pco2var <- ncvar_def(name="PCO2", units="microatmospheres", dim=timedim)

# Create NetCDF files
vars <- list(tempvar, salvar, pco2var)
ncout <- nc_create(ncfname, vars, force_v4=TRUE)

# Put variables values
ncvar_put(ncout, tempvar, socatsynth$`SST [deg.C]`)
ncvar_put(ncout, salvar, socatsynth$sal)
ncvar_put(ncout, pco2var, socatsynth$`fCO2rec [uatm]`)

# Add attributes to dimensions
ncatt_put(ncout, "TIME", "axis", "T")
ncatt_put(ncout, "TIME", "long_name", "time of measurement")
ncatt_put(ncout, "TIME", "valid_min", 0.0)

ncatt_put(ncout, "LATITUDE", "axis", "Y")
ncatt_put(ncout, "LATITUDE", "long_name", "latitude of measurement")
ncatt_put(ncout, "LATITUDE", "valid_min", -90.0)
ncatt_put(ncout, "LATITUDE", "valid_max", 90.0)

# Add attributes to variables
ncatt_put(ncout, "TEMP", "standard_name","sea_water_temperature")
ncatt_put(ncout, "TEMP", "_FillValue", 99999.)
# Add global attributes


# Global attributes
site_code=""
platform_code=""
data_mode="D"
title="SOCATv6 synthesis product"
summary
naming_authority
id
wmo_platform_code
source
network
keywords_vocabulary="SeaDataNet Parameter Discovery Vocabulary"
keywords
comment
----
area. Global?
geospatial_lat_min
geospatial_lat_max
geospatial_lat_units="degree_north"
geospatial_lon_min=
geospatial_lon_max
geospatial_lon_units="degree_east"
geospatial_vertical_min=5.0
geospatial_vertical_max=5.0 #SOCAT
geospatial_vertical_positive="down"
geospatial_vertical_units="meter"
time_coverage_start
time_coverage_end
time_coverage_duration #Use ISO 8601 P1Y P3M P1Y3M4D
cdm_data_type="profile" #?? for GLODAP
data_type
format_version
Conventions
netcdf_version
publisher_name
publisher_email
publisher_url
references
data_assembly_center=
update_interval="P1Y"
license
citation
acknowledgement
date_created
date_modified
history
processing_level
QC_indicator="excellent"
contributor_name="Dorothee C. Bakker; Siv K. Lauvset;..."
contributor_role
contributor_email
  


# Coordinate variables



# Data variables
<PARAM>:
standard_name
units
_FillValue
coordinates
long_name
QC_indicator
processing_level
valid_min
valid_max
comment
ancillary_variables
uncertainty
accuracy
precision
resolution
reference_scale

<PARAM>_QC:long_name ???
flag_values


---
CAPH air_pressure
DEPTH
DOX2 / DOXY
PCO2
PRES
PSAL
TEMP
