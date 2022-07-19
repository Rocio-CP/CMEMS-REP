import os
# Stablish where the files are / should be
input_files_dir = '/Users/rocio/Documents/'
files_path_remote = ('https://www.ncei.noaa.gov'
                     '/data/oceans/ncei/ocads/data/0237935/')
product_dir='/Users/rocio/Documents/INSITU_GLO_BGC_CARBON_DISCRETE_MY_013_050/'
# Create and store output directory.
if not os.path.isdir(product_dir):
    os.makedirs(os.path.join(product_dir,'DNT/'))
    os.makedirs(os.path.join(product_dir,'cmems_obs-ins_glo_bgc-car_my_glodap-gridded_irr/'))
    os.makedirs(os.path.join(product_dir,'cmems_obs-ins_glo_bgc-car_my_glodap-obs_irr/VESSEL'))
    os.makedirs(os.path.join(product_dir,'cmems_obs-ins_glo_bgc-car_my_glodap-obs_irr/ETC'))
    os.makedirs(os.path.join(product_dir,'cmems_obs-ins_glo_bgc-car_my_socat-gridded_irr/'))
    os.makedirs(os.path.join(product_dir,'cmems_obs-ins_glo_bgc-car_my_socat-obs_irr/VESSEL/'))
    os.makedirs(os.path.join(product_dir,'cmems_obs-ins_glo_bgc-car_my_socat-obs_irr/MOORING/'))
    os.makedirs(os.path.join(product_dir,'cmems_obs-ins_glo_bgc-car_my_socat-obs_irr/DRIFTER/'))
    os.makedirs(os.path.join(product_dir,'cmems_obs-ins_glo_bgc-car_my_socat-obs_irr/AUV/'))

output_files_dir = os.path.join(product_dir,'cmems_obs-ins_glo_bgc-car_my_glodap-obs_irr/')
# Create SOCAT
# Create GLODAP
# Create index and platform files
# Create DNT notes
# Push??