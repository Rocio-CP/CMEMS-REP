import os
import read_SOCAT
import read_GLODAP

# Stablish where the files are / should be
input_files_dir = '/Users/rocio/Documents/'
files_path_remote = ('https://www.ncei.noaa.gov'
                     '/data/oceans/ncei/ocads/data/0237935/')
product_dir = '/Users/rocio/Documents/INSITU_GLO_BGC_CARBON_DISCRETE_MY_013_050/'
# Create and store output directory.
if not os.path.isdir(product_dir):
    os.makedirs(os.path.join(product_dir, 'DNT/'))
    os.makedirs(os.path.join(product_dir, 'cmems_obs-ins_glo_bgc-car_my_glodap-gridded_irr/'))
    os.makedirs(os.path.join(product_dir, 'cmems_obs-ins_glo_bgc-car_my_glodap-obs_irr/BO'))
    os.makedirs(os.path.join(product_dir, 'cmems_obs-ins_glo_bgc-car_my_socat-gridded_irr/'))
    os.makedirs(os.path.join(product_dir, 'cmems_obs-ins_glo_bgc-car_my_socat-obs_irr/CO/'))
    os.makedirs(os.path.join(product_dir, 'cmems_obs-ins_glo_bgc-car_my_socat-obs_irr/DB/'))
    os.makedirs(os.path.join(product_dir, 'cmems_obs-ins_glo_bgc-car_my_socat-obs_irr/GL/'))
    os.makedirs(os.path.join(product_dir, 'cmems_obs-ins_glo_bgc-car_my_socat-obs_irr/MO/'))
    os.makedirs(os.path.join(product_dir, 'cmems_obs-ins_glo_bgc-car_my_socat-obs_irr/SD/'))


# Create SOCAT files
SOCAT_files = ['SOCATv2022.tsv','SOCATv2022_FlagE.tsv']
SOCAT_info_file='SOCATv2022CMEMS.tsv'
output_files_dir = os.path.join(product_dir, 'cmems_obs-ins_glo_bgc-car_my_socat-obs_irr/')

read_SOCAT.read_SOCAT_obs(input_files_dir, output_files_dir, SOCAT_files, SOCAT_info_file)

# Create GLODAP files
GLODAP_files = ['GLODAPv2.2021_Merged_Master_File.csv', 'EXPOCODES.txt', 'Dataset_DOIs.txt']
GLODAP_info_file = 'GLODAPv22021CMEMS.tsv'
output_files_dir = os.path.join(product_dir, 'cmems_obs-ins_glo_bgc-car_my_glodap-obs_irr/')

read_GLODAP.read_GLODAP_obs(input_files_dir, output_files_dir, GLODAP_files, GLODAP_info_file)

# Run tests, check if all ok
# Is this way bad? Probably. Will fix it? Eventually

# Remove preexisting reports
reports=[f for f in os.listdir(input_files_dir) if f.__contains__('check')]
for r in reports:
    if os.path.isfile(r):
        os.remove(r)

# File checkers
format_checker_path='/Users/rocio/Dropbox/Projects/CMEMS_INSTAC/current_FormatChecker'
content_checker_path='/Users/rocio/Documents/GitHub/CMEMS-REP/file-content-checker-main/'

for source in ['socat','glodap']:
    # Format checker
    os.system("cd "+ format_checker_path + "; "
              "for f in " + product_dir + "cmems_obs-ins_glo_bgc-car_my_" + source + "-obs_irr/*/*.nc; "
              "do ./control.csh $f"
              " >> " + input_files_dir +"format_check_" + source + "; done")

    # Cleanup the file: it's a bunch of individual xml printed together.
    f = open(input_files_dir +"format_check_" + source, "r")
    lines = f.readlines()
    f.close()

    f = open(input_files_dir +"format_check_" + source, "w")
    f.write('<?xml version="1.0"?>\n')
    f.write('<coriolis_function_report>\n')
    for line in lines:
        if 'xml version' not in line and 'coriolis_function' not in line:
            f.write(line)
    f.write('</coriolis_function_report>')
    f.close()

    # Content checker
    os.system("python3 " + content_checker_path + "Copernicus_InSituTAC_content_checker.py " +
              content_checker_path + "Copernicus_InSituTAC_content_checker.json " +
              product_dir + "/cmems_obs-ins_glo_bgc-car_my_" + source + "-obs_irr/*/*.nc"
              " > " + input_files_dir + "content_check_" + source)

    # CF checker (from IOOS)
    # Can't use CF checker (compliance-checker), because cf-units will not import

reports=[f for f in os.listdir(input_files_dir) if f.__contains__('check')]

# CHECK if it's status ok or something else (file_compliant in the format checker)
# try the xml route
from xml.dom.minidom import parse
for report_file in reports:
    parser = parse(input_files_dir + report_file)
    if 'format' in report_file:
        nc_file_tag=parser.getElementsByTagName('netcdf_file')
        compliant_tag = parser.getElementsByTagName('file_compliant')
    elif 'content' in report_file:
        nc_file_tag=parser.getElementsByTagName('file_name')
        compliant_tag = parser.getElementsByTagName('status')

    fail_count=0
    for n in range(0,compliant_tag.length):
        if compliant_tag[n].firstChild.data in ['no', 'ko']:
            print('FILE NOT INSTAC COMPLIANT '+ nc_file_tag[n].firstChild.data)
            fail_count=fail_count+1

    if not fail_count:
        print(report_file + ' ALL OK')

    # with open(input_files_dir + report_file, 'r') as file:
        # content = file.read()
        # if 'format' in report_file:
        #     # check if string present in a file
        #     if "<file_compliant>yes</file_compliant>" in content:
        #         print('Content checker OK in ' + report_file.split('/')[-1])
        #     else:
        #         print('Content checker NOT OK. Check the report ' + report_file)
        # elif 'content' in report_file:
        # # check if string present in a file
        #     if "<status>ok</status>" in content:
        #         print('Content checker OK in '+ report_file.split('/')[-1])
        #     else:
        #         print('Content checker NOT OK. Check the report '+ report_file)

#import sys
#sys.path.insert(0, '/Users/rocio/Documents/GitHub/CMEMS-REP/compliance-checker-master/compliance_checker')
#from compliance_checker.runner import ComplianceChecker

#
#
# os.system("cd "+ format_checker_path + "; "
#           "for f in /Users/rocio/Documents/INSITU_GLO_BGC_CARBON_DISCRETE_MY_013_050/cmems_obs-ins_glo_bgc-car_my_glodap-obs_irr/*/*.nc; do ./control.csh $f"
#           " >> " + input_files_dir +"formatcheckoutGLODAP"+
#           "; done")
#
# # Content checker
# os.system("python3 /Users/rocio/Documents/GitHub/CMEMS-REP/file-content-checker-main/Copernicus_InSituTAC_content_checker.py "
#           "/Users/rocio/Documents/GitHub/CMEMS-REP/file-content-checker-main/Copernicus_InSituTAC_content_checker.json "
#           "/Users/rocio/Documents/INSITU_GLO_BGC_CARBON_DISCRETE_MY_013_050/cmems_obs-ins_glo_bgc-car_my_socat-obs_irr/**/*.nc"
#           " > " + input_files_dir +"socat_report.csv")
#
# os.system("python3 /Users/rocio/Documents/GitHub/CMEMS-REP/file-content-checker-main/Copernicus_InSituTAC_content_checker.py "
#           "/Users/rocio/Documents/GitHub/CMEMS-REP/file-content-checker-main/Copernicus_InSituTAC_content_checker.json "
#           "/Users/rocio/Documents/INSITU_GLO_BGC_CARBON_DISCRETE_MY_013_050/cmems_obs-ins_glo_bgc-car_my_glodap-obs_irr/**/*.nc"
#           " > " + input_files_dir +"glodap_report.csv")
#
# # CHECK if it's status ok or something else (file_compliant in the format checker)
# for report_file in reports:
#     with open(report_file, 'r') as file:
#         content = file.read()
#         if 'format' in report_file:
#             # check if string present in a file
#             if "<file_compliant>no</file_compliant>" in content:
#                 print('Content checker OK in ' + report_file.split('/')[-1])
#             else:
#                 print('Content checker NOT OK. Check the report ' + report_file)
#         elif 'content' in report_file:
#         # check if string present in a file
#             if "<status>ok</status>" in content:
#                 print('Content checker OK in '+ report_file.split('/')[-1])
#             else:
#                 print('Content checker NOT OK. Check the report '+ report_file)
#
# # Create index and platform files
#
# # Create DNT notes
#
# Push??