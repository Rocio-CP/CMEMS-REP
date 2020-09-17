#!/bin/bash

shopt -s nullglob

# Folder structure
for dataset in */
do

# Skip DNT folders
if [[ ${dataset} == *"DNT"*  ]]; then
continue
fi
if [[ ${dataset} == *"GRIDDED"*  ]]; then
continue
fi

# Header of index file
echo "# Title : In Situ  platforms catalog" > ${dataset}index_platform.txt
echo "# Description : catalog of available In Situ platforms" >> ${dataset}index_platform.txt
echo "# Project : Copernicus" >> ${dataset}index_platform.txt
echo "# Format version : 1.0" >> ${dataset}index_platform.txt
datetimenow_index=$(date +%Y%m%d%H%M%S)
echo "# Date of update : ${datetimenow_index}" >> ${dataset}index_platform.txt
echo "#" >> ${dataset}index_platform.txt
echo "platform_code,creation_date,update_date,wmo_platform_code,data_source,institution,institution_edmo_code,parameter,last_latitude_observation,last_longitude_observation,last_date_observation" >> ${dataset}index_platform.txt

# List netcdf files. 
for f in ${dataset}*/*.nc
do
platformcode=$(ncks -M $f | grep -m 1 "platform_code" | cut -d\" -f2)
creationdate=$(ncks -M $f | grep -m 1 "time_coverage_start" | cut -d\" -f2)
updated=$(ncks -M $f | grep -m 1 "date_update" | cut -d\" -f2)
wmo=$(ncks -M $f | grep -m 1 "wmo_platform_code" | cut -d\" -f2)
provider=$(ncks -M $f | grep -m 1 "institution " | cut -d\" -f2 | cut -d\" -f1)
provideredmo=$(ncks -M $f | grep -m 1 "institution_edmo_code" | cut -d\" -f2 | cut -d\" -f1)
variables=$(ncks -M $f | grep "_QC" | cut -d_  -f1 | cut -dr -f2 | cut -d\  -f6 | tr '\n' ' ')
last_lat=$(ncks -M $f | grep -m 1 "last_latitude_observation" | cut -d\  -f7)
last_lon=$(ncks -M $f | grep -m 1 "last_longitude_observation" | cut -d\  -f7)
last_date=$(ncks -M $f | grep -m 1 "last_date_observation" | cut -d\  -f7)

echo "${platformcode},${creationdate},${updated},${wmo},${f},${provider},${provideredmo},${variables},${last_lat},${last_lon},${last_date}" >> ${dataset}index_platform.txt
done



done



