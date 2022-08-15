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

# If index file exists, remove
if [ -f "${dataset}index_file.txt" ]; then
    rm "${dataset}index_file.txt"
fi

# Header of index file
echo "# Title : in-situ files catalog" > ${dataset}index_file.txt
echo "# Description : catalog of available in-situ files" >> ${dataset}index_file.txt
echo "# Project : Copernicus Marine In Situ TAC" >> ${dataset}index_file.txt
echo "# Format version : 1.4" >> ${dataset}index_file.txt
datetimenow_index=$(date +%Y-%m-%dT%H:%M:%SZ)
echo "# Date of update : ${datetimenow_index}" >> ${dataset}index_file.txt
echo "product_id,file_name,geospatial_lat_min,geospatial_lat_max,geospatial_lon_min,geospatial_lon_max,time_coverage_start,time_coverage_end,institution,date_update,data_mode,parameters" >> ${dataset}index_file.txt

# List netcdf files. 
for f in ${dataset}*/*.nc
do
#catid=COP-GL-01
latmin=$(ncks -M $f | grep -m 1 "lat_min" | grep -oE "[-]?[0-9]+\.?[0-9]+")
lonmin=$(ncks -M $f | grep -m 1 "lon_min" | grep -oE "[-]?[0-9]+\.?[0-9]+")
latmax=$(ncks -M $f | grep -m 1 "lat_max" | grep -oE "[-]?[0-9]+\.?[0-9]+")
lonmax=$(ncks -M $f | grep -m 1 "lon_max" | grep -oE "[-]?[0-9]+\.?[0-9]+")
timestart=$(ncks -M $f | grep -m 1 "time_coverage_start" | cut -d\" -f2)
timeend=$(ncks -M $f | grep -m 1 "time_coverage_end" | cut -d\" -f2)
updated=$(ncks -M $f | grep -m 1 "date_update" | cut -d\" -f2)
instit=$(ncks -M $f | grep -w "institution" | cut -d\" -f2 | tr ',' '-')
if [[ $dataset == *"obs"* ]]; then
variables=$(ncks -M $f | grep "_QC" | cut -d_  -f1 | awk -F\  '{print $NF}' | grep -vE 'TIME|POSITION|DEPH' | tr '\n' ' ')
#variables=$(ncks -M $f | grep "_QC" | cut -d_  -f1 | cut -dr -f2 | cut -d\  -f2 | tr '\n' ' ')
elif [[ $dataset == *"gridded"* ]]; then
variables=$(ncks --trd -m $f | grep -E ': type' | grep -v 'lat' | grep -v 'lon' | cut -f 1 -d ' ' | sed 's/://' | tr '\n' ' ')
fi

echo "COP-GL-01,ftp://my.cmems-du.eu/Core/INSITU_GLO_BGC_CARBON_DISCRETE_MY_013_050/${f},${latmin},${latmax},${lonmin},${lonmax},${timestart},${timeend},${instit},${updated},D,${variables}" >> ${dataset}index_file.txt
done

done


