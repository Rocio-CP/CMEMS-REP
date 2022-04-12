# SOCAT import. From local .tsv file or ERDDAP

import os
import pandas as pd
import time
from erddapy import ERDDAP


headerlines = 0
"https://www.ncei.noaa.gov/data/oceans/ncei/ocads/data/0235360/"

##### -----------------------------
# Don't touch from here
file = source+'.tsv'
SOCATcolheadertextshort = 'Expocode\tversion\tSource_DOI\tQC_Flag'
SOCATmetacolheadertextshort = 'Expocode\tversion\tDataset'
# Dictionaries
flagaccuracy = {"A": 2.0, "B": 2.0, "C": 5.0, "D": 5.0, "E": 10.0}

# Read from local files (synthesis + flagE)
if (datafrom == 'local'):
    SOCATmetacolheadertextshort = 'Expocode\tversion\tDataset'
    SOCATcolheadertextshort = 'Expocode\tversion\tSource_DOI\tQC_Flag'

    SOCAT_files = os.listdir(os.path.join(input_dir, 'SOCAT'))

    for file in SOCAT_files:
        filepath = os.path.join(input_dir, 'SOCAT', file)

        # Read metadata header for Cruise Flags and find the number of
        # headerlines before the data
        separator = '\t'
        line = ''
        headerlines = -1
        metaline = ''
        metaheaderlines = -1
        f = open(filepath)

        while SOCATmetacolheadertextshort not in metaline:
            metaline = f.readline()
            metaheaderlines = metaheaderlines + 1  # Where metadata lines start
            # Find SOCAT collection DOI, while reading the file. The DOI is
            # wrong AND have to remove the \n!!
            #if ('DOI of the entire SOCAT collection:' in metaline):
            #  socatdoi = metaline.rsplit(' ', 1)[1]

        endmetaheaderlines = metaheaderlines
        while '\t' in metaline:
            metaline = f.readline()
            endmetaheaderlines = endmetaheaderlines + 1  # End of metadata lines
        # Create metadata dataframe
        metainfoAD = pd.read_csv(filepath, sep = '\t',
          skiprows = metaheaderlines,
          nrows = endmetaheaderlines - metaheaderlines - 1)
        # Find where data columns start
        headerlines = headerlines + endmetaheaderlines
        while SOCATcolheadertextshort not in line:
            line = f.readline()
            headerlines = headerlines + 1
        f.close()

        # Read SOCAT data in dataframe
        print(file + " file has " + str(headerlines) + " header lines")
        start_time = time.time()  # Time the script
        ddtype = {0: str, 2: str}  # add type str to columns 0 and 2
        # Read the SOCAT file into a pandas dataframe
        tempdf1 = pd.read_csv(filepath, sep = separator, skiprows = headerlines,
            dtype = ddtype)
        print("--- %s seconds ---" % (time.time() - start_time))
        print(file + " data frame has " + str(len(tempdf1)) + " lines")

        # Give all the same variable names
        tempdf1.rename(
            columns = {'Expocode': vardict['id'],
              'Source_DOI': vardict['doi'],
              'latitude [dec.deg.N]': vardict['lat'],
              'longitude [dec.deg.E]': vardict['lon'],
              'sample_depth [m]': vardict['dep'],
              'SST [deg.C]': vardict['temp'],
              'sal': vardict['sal'],
              'fCO2rec [uatm]': vardict['fco2w'],
              'fCO2rec_flag': vardict['fco2wf']},
            inplace=True)

        # Create date python object
        tempdtframe = pd.DataFrame(
            {'year': tempdf1['yr'], 'month': tempdf1['mon'],
            'day': tempdf1['day'], 'hour': tempdf1['hh'],
            'minute': tempdf1['mm'], 'seconds': tempdf1['ss']})
        tempdf1['DATEVECTOR1'] = pd.to_datetime(tempdtframe, utc = True)

        # Transform longitude to +-180
        tempdf1.loc[tempdf1[vardict['lon']] > 180, vardict['lon']] = tempdf1[
          vardict['lon']] - 360

        # Subset the dataset HERE (cruise flag assignment takes A LONG TIME).
        if subset :
            tempdf1 = tempdf1[
              (tempdf1['DATEVECTOR1'] >= pd.to_datetime(mindate)) &
              (tempdf1['DATEVECTOR1'] <= pd.to_datetime(maxdate)) &
              (tempdf1[vardict['lat']] >= minlat) &
              (tempdf1[vardict['lat']] <= maxlat) &
              (tempdf1[vardict['lon']] >= minlon) &
              (tempdf1[vardict['lon']] <= maxlon)].copy()
            tempdf1.reset_index(drop = True, inplace = True)
            print(len(tempdf1))

        # Cruise flags / accuracies
        if ('FlagE' in file):  # Pointless to loop when all has the same flag
            tempdf1['Cruise_flag'] = 'E'

        elif ("FlagE" not in file):
            # Read cruise flag from metadata lines in SOCAT synthesis A-D
            tempdf1['Cruise_flag'] = 'X'
            #tempdf[vardict['fco2wac']] = 0.0

            counter = 0
            #allexpoflags = tempdf[['Expocode', 'Cruise_flag']]
            #print(allexpoflags.shape)

            # Assign Cruise flags A-D
            start_time = time.time()
            for expocode in metainfoAD['Expocode']:
                cruiseflag = metainfoAD['QC Flag'][metainfoAD[
                  'Expocode'] == expocode]

                # this line is very slow
                tempdf1['Cruise_flag'].values[
                tempdf1[vardict['id']] == expocode] = cruiseflag

                counter = counter + 1
                # print(counter)
            print("--- %s seconds for cruise flag assignment ---"
             % (time.time() - start_time))

        # Merge synthesis and FlagE dataframes in one
        if 'tempdf' not in globals():
            tempdf = tempdf1
        else:
            tempdf = tempdf.append(tempdf1)

# Create UNIXDATE and ISO DATEVECTOR
tempdf[vardict['unixd']] = tempdf['DATEVECTOR1'].astype('int64') // 10 ** 9
tempdf[vardict['datevec']] = tempdf[
  'DATEVECTOR1'].dt.strftime('%Y-%m-%dT%H:%M:%SZ')

# Assign accuracies following cruise flags
tempdf[vardict['fco2wac']] = 0.0
for key in flagaccuracy:
    tempdf[vardict['fco2wac']].values[
      tempdf['Cruise_flag'] == key] = flagaccuracy[key]
# Flag fco2 as measured
tempdf[vardict['fco2wc']] = 0

# Estimate alkalinity from salinity, and then, estimate ph and dic

# Assign SOCAT DOI if Source DOI is missing
tempdf.loc[tempdf[vardict['doi']].isna(), vardict['doi']] = socatdoi

# Add source (SOCAT, GLODAP, ARGO, etc...)
tempdf['SOURCE'] = source

# Rename and reset indices
printdf = tempdf
printdf.reset_index(drop=True, inplace=True)

print('SOCAT frame size is ')
print(printdf.shape)