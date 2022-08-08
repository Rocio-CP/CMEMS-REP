#!/bin/bash

# skip empty folders
shopt -s nullglob

# Loop counter. Used to add seconds to the DNT timestamp if two are created at the same time 
indcount=0

# Folder structure
cd /Users/rocio/Documents/templocal/CARBON_REP_202212/INSITU_GLO_BGC_CARBON_DISCRETE_MY_013_050/
for dataset in */
do

# Skip DNT folders
if [[ ${dataset} == *"DNT"*  ]]; then
continue
fi

indcount=$(( $indcount + 1 ))

# DNT file header
datetimenow=$(date +%Y%m%dT%H%M%SZ)
datetimestart=$(date -v -2H +%Y%m%dT%H%M%SZ)
datetimestop=$(date -v -10M +%Y%m%dT%H%M%SZ)

# Add seconds to the DNT file name. Rewrite with a while loop
if [ -f ./DNT/INSITU_GLO_BGC_CARBON_DISCRETE_MY_013_050_P${datetimenow}.xml ]; then
datetimenow=$(date -v +${indcount}S +%Y%m%dT%H%M%SZ)
fi

echo "<?xml version=\"1.0\" ?>" > ./DNT/INSITU_GLO_BGC_CARBON_DISCRETE_MY_013_050_P${datetimenow}.xml
echo "<delivery PushingEntity=\"CopernicusMarine-InSitu-Global\" date=\"${datetimenow}\" product=\"INSITU_GLO_BGC_CARBON_DISCRETE_MY_013_050\">" >> ./DNT/INSITU_GLO_BGC_CARBON_DISCRETE_MY_013_050_P${datetimenow}.xml
echo $'\t'"<dataset DatasetName=\"${dataset%%/*}_201904\">" >> ./DNT/INSITU_GLO_BGC_CARBON_DISCRETE_MY_013_050_P${datetimenow}.xml

for f in ${dataset}*/*.nc
do checksum=$(md5 $f | cut -d\  -f4)
echo $'\t'$'\t'"<file Checksum=\"${checksum}\" FileName=\"$(echo $f | cut -d/ -f2-3)\" FinalStatus=\"Delivered\" StartUploadTime=\"${datetimestart}\" StopUploadTime=\"${datetimestop}\"/>" >> ./DNT/INSITU_GLO_BGC_CARBON_DISCRETE_MY_013_050_P${datetimenow}.xml
done

checksumindex=$(md5 ${dataset}/index_file.txt | cut -d\  -f4)
echo $'\t'$'\t'"<file Checksum=\"${checksumindex}\" FileName=\"index_file.txt\" FinalStatus=\"Delivered\" StartUploadTime=\"${datetimestart}\" StopUploadTime=\"${datetimestop}\"/>" >> ./DNT/INSITU_GLO_BGC_CARBON_DISCRETE_MY_013_050_P${datetimenow}.xml

echo $'\t'"</dataset>" >> ./DNT/INSITU_GLO_BGC_CARBON_DISCRETE_MY_013_050_P${datetimenow}.xml
echo "</delivery>" >> ./DNT/INSITU_GLO_BGC_CARBON_DISCRETE_MY_013_050_P${datetimenow}.xml
done

