import netCDF4
import re
import datetime
import os
import shutil
import json


def gridded_INSTAC(listgriddedfiles, output_files_dir):
    for old_nc_filename in listgriddedfiles:
        # Copy the netcdf file, rename and move to destination
        nc_filename = os.path.join(output_files_dir, old_nc_filename.split('/')[-1][0:-3] + "_CMS_INSTAC.nc")
        print(old_nc_filename)
        print(nc_filename)
        shutil.copy(old_nc_filename, nc_filename)

        ncfile = netCDF4.Dataset(nc_filename, mode='r+')

        # Add global attributes
        global_attributes = gridded_global_attributes_dictionary(nc_filename)
        for key, value in global_attributes.items():
            ncfile.setncattr(key, value)
        # Remove unnecesary global attributes (?)

        # Make CF-compliant
        if 'socat'.casefold() in nc_filename:
            # Change salinity units
            if 'qrtrdeg' in nc_filename:
                prefix = 'coast_';suffix = ''
            elif 'monthly' in nc_filename:
                prefix = '';suffix = ''
            elif 'yearly' in nc_filename:
                prefix = '';suffix = '_year'
            elif 'decadal' in nc_filename:
                prefix = '';suffix = '_decade'

            ncfile[prefix + "salinity_ave_weighted" + suffix].setncattr("units", 1)
            ncfile[prefix + "salinity_ave_unwtd" + suffix].setncattr("units", 1)
            ncfile[prefix + "salinity_max_unwtd" + suffix].setncattr("units", 1)
            ncfile[prefix + "salinity_min_unwtd" + suffix].setncattr("units", 1)

        elif 'glodap'.casefold() in nc_filename:
            # Change dimensions variables attributes
            ncfile['lon'].setncattr('standard_name', 'longitude')
            ncfile['lon'].setncattr('units', 'degree_east')
            ncfile['lat'].setncattr('standard_name', 'latitude')
            ncfile['lat'].setncattr('units', 'degree_north')
            ncfile['Depth'].setncattr('standard_name', 'depth')

            # Change units
            for v in ncfile.variables.keys():
                varatts = ncfile[v].ncattrs()
                if 'units' in varatts:
                    if ncfile[v].units == 'micro-mol kg-1':
                        ncfile[v].setncattr('units', 'umol kg-1')
                    elif ncfile[v].units == 'degrees celcius':
                        ncfile[v].setncattr('units', 'degree_C')

        ncfile.close()


def gridded_global_attributes_dictionary(nc_filename):
    with open("CMEMS_INSTAC_metadata.json", "r") as template:
        all_attributes = json.load(template)
    template.close()

    ncfile = netCDF4.Dataset(nc_filename, mode='r')

    all_attributes['globalatt'][0]['id'] = nc_filename.split('/')[-1][0:-3]  # nc_filename comes with path
    all_attributes['globalatt'][0]["data_type"] = "Copernicus Marine gridded data"
    all_attributes['globalatt'][0]["format_version"] = "1.0"
    all_attributes['globalatt'][0]["cdm_data_type"] = "grid"
    all_attributes['globalatt'][0]["netcdf_version"] = "netCDF-4"
    all_attributes['globalatt'][0]["author"] = ""
    all_attributes['globalatt'][0]["naming_authority"] = "Copernicus Marine"
    all_attributes['globalatt'][0]["data_assembly_center"] = "BERGEN"

    # Geospatial and time limits
    all_attributes['globalatt'][0]['geospatial_lat_min'] = -90
    all_attributes['globalatt'][0]['geospatial_lat_max'] = 90
    all_attributes['globalatt'][0]['geospatial_lon_min'] = -180
    all_attributes['globalatt'][0]['geospatial_lon_max'] = 180
    all_attributes['globalatt'][0]['latitude_resolution'] = 1.0
    all_attributes['globalatt'][0]['longitude_resolution'] = 1.0
    if 'coast' in nc_filename:
        all_attributes['globalatt'][0]['latitude_resolution'] = 0.25
        all_attributes['globalatt'][0]['longitude_resolution'] = 0.25

    # History attributes
    all_attributes['globalatt'][0]['date_update'] = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%SZ")

    # Data product-specific attributes
    if 'socat'.casefold() in nc_filename.casefold():

        all_attributes['globalatt'][0]['geospatial_vertical_min'] = 5
        all_attributes['globalatt'][0]['geospatial_vertical_max'] = 5
        all_attributes['globalatt'][0]['time_coverage_start'] = '1970-01-01T00:00:00Z'
        all_attributes['globalatt'][0]['time_coverage_end'] = '2021-02-25T07:57:00Z'  # CHECK!!

        all_attributes['globalatt'][0]["title"] = all_attributes['globalatt'][0]["title"] + " gridded - SOCATv2022"
        all_attributes['globalatt'][0]["references"] = all_attributes['globalatt'][0][
                                                           "references"] + " https://socat.info"
        all_attributes['globalatt'][0]["citation"] = all_attributes['globalatt'][0][
                                                         "citation"] + "  SOCAT is described in Bakker et al. (2016) https://doi.org/10.5194/essd-8-383-2016; traceable citations are essential for justifying and sustaining the effort."
        all_attributes['globalatt'][0]["doi"] = "https://doi.org/10.25921/1h9f-nb73"
        all_attributes['globalatt'][0][
            "summary"] = 'Surface Ocean Carbon Atlas (SOCAT) Gridded (binned) SOCAT observations, with a spatial grid of 1x1 degree and yearly in time. The gridded fields are computed from the monthly 1-degree gridded data, which uses only SOCAT datasets with QC flags of A through D and SOCAT data points flagged with WOCE flag values of 2. This yearly data is computed using data from the start to the end of each year as described in the summary attribute of each variable.'

        # Get date when the grid file was created by SOCAT from the global attributes
        socathistory = ncfile.getncattr('history')
        #match = re.search(r'\d{2}-[a-zA-Z]{3}-\d{2}', socathistory)
        #datetime_obj = datetime.datetime.strptime(match.group(), '%d-%b-%y').date()
        datetime_obj = datetime.datetime.strptime(socathistory.split(' ')[-1], '%d-%b-%y').date()
        datetime_str = datetime_obj.strftime("%Y-%m-%dT%H:%M:%SZ")
        all_attributes['globalatt'][0]['history'] = datetime_str + ": Created by " + socathistory + "\n" \
                                                    + datetime.datetime.now().strftime(
            "%Y-%m-%dT%H:%M:%SZ") + ": Created CMS INSTAC file format"

        socattitle = ncfile.getncattr('title')
        all_attributes['globalatt'][0]["SOCAT_Title"] = socattitle
        all_attributes['globalatt'][0][
            "caution"] = 'NO INTERPOLATION WAS PERFORMED. SIGNIFICANT BIASES ARE PRESENT IN THESE GRIDDED RESULTS DUE TO THE ARBITRARY AND SPARSE LOCATIONS OF DATA VALUES IN BOTH SPACE AND TIME.'

    elif 'glodap'.casefold() in nc_filename.casefold():

        all_attributes['globalatt'][0]['geospatial_vertical_min'] = 0
        all_attributes['globalatt'][0]['geospatial_vertical_max'] = 5500
        all_attributes['globalatt'][0]['time_coverage_start'] = '1972-01-01T00:00:00Z'
        all_attributes['globalatt'][0]['time_coverage_end'] = '2013-12-31T23:59:59Z'

        all_attributes['globalatt'][0]["title"] = all_attributes['globalatt'][0]["title"] + " gridded - GLODAPv2"
        all_attributes['globalatt'][0]["references"] = all_attributes['globalatt'][0][
                                                           "references"] + " https://glodap.info"
        all_attributes['globalatt'][0]["citation"] = all_attributes['globalatt'][0][
                                                         "citation"] + " The GLODAPv2 gridded climatologies are described in Lauvset et al. (2016) doi:10.5194/essd-8-325-2016; traceable citations are essential for justifying and sustaining the effort."
        all_attributes['globalatt'][0]["doi"] = "10.5194/essd-8-325-2016"
        all_attributes['globalatt'][0][
            "summary"] = '1 X 1 global mapped field from the GLODAPv2 data product. Mapping is performed using the DIVA software (Troupin et al., 2012). Error fields are calculated using the clever poor mans error calculation method in DIVA (Beckers et al., 2014). The error fields represent the mapping error only, and does not include measurement or calculation uncertainties in the input data.'
        all_attributes['globalatt'][0]['history'] = "2016-05-12T18:41:01: Created by Siv K. Lauvset \n" \
                                                    + "2019-05-24T11:23:50Z: Created CMEMS INSTAC file \n" \
                                                    + datetime.datetime.now().strftime(
            "%Y-%m-%dT%H:%M:%SZ") + ": Updated to CMS INSTAC file format"

    ncfile.close()

    return all_attributes['globalatt'][0]
