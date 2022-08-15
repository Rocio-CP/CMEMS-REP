#!/bin/bash

shopt -s nullglob

# Folder structure
cd /Users/rocio/Documents/templocal/CARBON_REP_202212/INSITU_GLO_BGC_CARBON_DISCRETE_MY_013_050/
for dataset in */
do

# Skip DNT folders
if [[ ${dataset} == *"DNT"*  ]]; then
continue
fi
if [[ ${dataset} == *"gridded"*  ]]; then
continue
fi

# If platform index file exists, remove
if [ -f "${dataset}index_platform.txt" ]; then
    rm "${dataset}index_platform.txt"
fi


# Header of index file
echo "# Title : in-situ platforms catalog" > ${dataset}index_platform.txt
echo "# Description : catalog of available in-situ platforms" >> ${dataset}index_platform.txt
echo "# Project : Copernicus Marine In Situ TAC" >> ${dataset}index_platform.txt
echo "# Format version : 1.4" >> ${dataset}index_platform.txt
datetimenow_index=$(date +%Y%m%d%H%M%S)
echo "# Date of update : ${datetimenow_index}" >> ${dataset}index_platform.txt
echo "platform_code,date_creation,date_update,wmo_platform_code,data_source,institution,institution_edmo_code,parameters,last_latitude_observation,last_longitude_observation,last_date_observation" >> ${dataset}index_platform.txt

# List netcdf files. 
for f in ${dataset}*/*.nc
do
platformcode=$(ncks -M $f | grep -m 1 "platform_code" | cut -d\" -f2)
creationdate=$(ncks -M $f | grep -m 1 "time_coverage_start" | cut -d\" -f2)
updated=$(ncks -M $f | grep -m 1 "date_update" | cut -d\" -f2)
wmo=$(ncks -M $f | grep -m 1 "wmo_platform_code" | cut -d\" -f2)
instit=$(ncks -M $f | grep -w "institution" | cut -d\" -f2 | tr ',' '-')
provideredmo=$(ncks -M $f | grep -m 1 "institution_edmo_code" | cut -d\" -f2 | cut -d\" -f1)
variables=$(ncks -M $f | grep "_QC" | cut -d_  -f1 | awk -F\  '{print $NF}' | grep -vE 'TIME|POSITION|DEPH' | tr '\n' ' ')
last_lat=$(ncks -M $f | grep -m 1 "last_latitude_observation" | cut -d\  -f7)
last_lon=$(ncks -M $f | grep -m 1 "last_longitude_observation" | cut -d\  -f7)
last_date=$(ncks -M $f | grep -m 1 "last_date_observation" | cut -d\" -f2)

echo "${platformcode},${creationdate},${updated},${wmo},${f},${instit},${provideredmo},${variables},${last_lat},${last_lon},${last_date}" >> ${dataset}index_platform.txt
done



done



