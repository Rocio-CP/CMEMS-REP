#!/bin/bash

shopt -s nullglob

# Folder structure
for dataset in */
do

# Skip DNT folders
if [[ ${dataset} == *"DNT"*  ]]; then
continue
fi

# Header of index file
echo "# Title : Carbon in-situ observations catalog" > ${dataset}index_file.txt
echo "# Description : catalog of available in-situ observations per platform." >> ${dataset}index_file.txt
echo "# Project : Copernicus" >> ${dataset}index_file.txt
echo "# Format version : 1.0" >> ${dataset}index_file.txt
datetimenow_index=$(date +%Y%m%d%H%M%S)
echo "# Date of update : ${datetimenow_index}" >> ${dataset}index_file.txt
echo "#" >> ${dataset}index_file.txt
echo "catalog_id,file_name,geospatial_lat_min,geospatial_lat_max,geospatial_lon_min,geospatial_lon_max,time_coverage_start,time_coverage_end,provider,date_update,data_mode,parameters" >> ${dataset}index_file.txt

# List netcdf files. 
for f in ${dataset}*/*.nc
do
catid=COP-GLOBAL-01
latmin=$(ncks -M $f | grep -m 1 "lat_min" | cut -d\" -f2)
lonmin=$(ncks -M $f | grep -m 1 "lon_min" | cut -d\" -f2)
latmax=$(ncks -M $f | grep -m 1 "lat_max" | cut -d\" -f2)
lonmax=$(ncks -M $f | grep -m 1 "lon_max" | cut -d\" -f2)
timestart=$(ncks -M $f | grep -m 1 "time_coverage_start" | cut -d\" -f2)
timeend=$(ncks -M $f | grep -m 1 "time_coverage_end" | cut -d\" -f2)
updated=$(ncks -M $f | grep -m 1 "date_update" | cut -d\" -f2)
if [[ $dataset == *"OBSERVATION"* ]]; then
variables=$(ncks -M $f | grep "_DM" | cut -d_  -f1 | cut -dr -f2 | cut -d\  -f2 | tr '\n' ' ')
elif [[ $dataset == *"GRIDDED"* ]]; then
variables=$(ncks --trd -m $f | grep -E ': type' | grep -v 'lat' | grep -v 'lon' | cut -f 1 -d ' ' | sed 's/://' | tr '\n' ' ')
fi

echo "${catid},ftp://my.cmems-du.eu/Core/INSITU_GLO_CARBON_REP_OBSERVATIONS_013_050/${f},${latmin},${latmax},${lonmin},${lonmax},${timestart},${timeend},University of Bergen Geophysical Institute,${updated},D,${variables}" >> ${dataset}index_file.txt
done

done


