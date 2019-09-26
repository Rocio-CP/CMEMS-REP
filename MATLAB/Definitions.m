clc; clear all; close all;
workrootdir='/Users/rpr061/Dropbox/BCDC_Projects/CMEMS_INSTAC/REP_Products/'

inputproducts={'SOCATv6','GLODAPv2'};
inputtypes={'observations','gridded'};

inputfiles


%% SOCAT
if inputproduct=='SOCATv6'
    qcmanual='SOCAT Quality Control Cookbook v3: https://www.socat.info/wp-content/uploads/2017/04/2015_SOCAT_QC_Cookbook_v3.pdf';
    globcomment= ['Further details can be found in Bakker et al. (2016). A multi-decade record of high-quality fCO2 data in version 3 of the ',...
        'Surface Ocean CO2 Atlas (SOCAT), Earth Syst. Sci. Data, 8, 383-413, https://doi.org/10.5194/essd-8-383-2016. ',...
        'and Bakker et al. (2018): Surface Ocean CO2 Atlas (SOCAT) V6. PANGAEA, https://doi.org/10.1594/PANGAEA.890974'];
    refwebsite='https://www.socat.info';
        
    if inputtype=='observations'
        inputdir=[workrootdir, 'SOCATv6/workspace/'];
        workdir=inputdir;
        fileload=[inputdir,'SOCATv6.mat'];
        
    % Which variables are included in the product
    invars={''};
    %vars={'bath','temp','sal','oxy','nitra','nitri','phos','sili','ph','ph25','dic','alk','doc','don','tdn','chla'};
    osvars={'BATH','TEMP','PSAL','DOX2','NTAW','NTIW','PHOW','SLCW','PHPH','PH25','TICW','ALKW','CORG','NODW','NT1D','CPHL'};
    osunits={'meters','degrees_C','0.001','micromol_per_kg',...
        'micromol_per_kg','micromol_per_kg','micromol_per_kg',...
        'micromol_per_kg','','',...
        'micromol_per_kg','micromol_per_kg','micromol_per_kg','micromol_per_kg','micromol_per_kg','microgram_per_l'};
    oslongname={'Bottom depth','Sea temperature','Practical salinity','Dissolved Oxygen',...
        'Nitrate (NO3-N)','Nitrite (NO2-N)','Phosphate (PO-P)',...
        'Silicate (SiO4-Si)','Ph','Ph at 25 degrees and 0 dbar',...
        'Dissolved inorganic carbon','Total alkalinity','Dissolved organic carbon','Dissolved organic nitrogen',...
        'Total dissolved nitrogen','Chlorophyll-a'};
    osstandardname={'bottom_depth','sea_water_temperature','sea_water_practical_salinity','moles_of_oxygen_per_unit_mass_in_sea_water',...
        'moles_of_nitrate_per_unit_mass_in_sea_water','moles_of_nitrite_per_unit_mass_in_sea_water','moles_of_phosphate_per_unit_mass_in_sea_water',...
        'moles_of_silicate_per_unit_mass_in_sea_water','sea_water_ph_reported_on_total_scale','sea_water_ph_reported_on_total_scale_at_25_degrees_and_0_dbar',...
        'moles_of_dissolved_inorganic_carbon_per_unit_mass_in_sea_water','sea_water_alkalinity_per_unit_mass',...
        'moles_of_dissolved_organic_carbon_per_unit_mass_in_sea_water','moles_of_dissolved_organic_nitrogen_per_unit_mass_in_sea_water',...
        'moles_of_dissolved_total_nitrogen_per_unit_mass_in_sea_water','mass_concentration_of_chlorophyll_in_sea_water'};
    
        
    elseif inputtype=='grid'
        inputdir=[workrootdir, 'SOCATv6/original/'];
        workdir=[workrootdir, 'SOCATv6/workspace/'];
        
    end
    
    
    %% GLODAP
elseif inputproduct == 'GLODAPv2'
    refwebsite='https://www.glodap.info';
    
    
    if inputtype =='observations'
        qcmanual='Olsen et. al (2016), ESSD';
        globcomment={['Further details can be found in Olsen et al. (2016).',...
            'The Global Ocean Data Analysis Project version 2 (GLODAPv2) – an internally consistent data product for the world ocean, ',...
            'Earth Syst. Sci. Data, 8, 297-323, https://doi.org/10.5194/essd-8-297-2016']};
        
        
        % Which variables are included in the product
        invars={'bottomdepth','temperature','salinity','oxygen','nitrate','nitrite','phosphate','silicate','phtsinsitutp','phts25p0','tco2','talk','doc1','don','tdn','chla'};
        %vars={'bath','temp','sal','oxy','nitra','nitri','phos','sili','ph','ph25','dic','alk','doc','don','tdn','chla'};
        osvars={'BATH','TEMP','PSAL','DOX2','NTAW','NTIW','PHOW','SLCW','PHPH','PH25','TICW','ALKW','CORG','NODW','NT1D','CPHL'};
        osunits={'meters','degrees_C','0.001','micromol_per_kg',...
            'micromol_per_kg','micromol_per_kg','micromol_per_kg',...
            'micromol_per_kg','','',...
            'micromol_per_kg','micromol_per_kg','micromol_per_kg','micromol_per_kg','micromol_per_kg','microgram_per_l'};
        oslongname={'Bottom depth','Sea temperature','Practical salinity','Dissolved Oxygen',...
            'Nitrate (NO3-N)','Nitrite (NO2-N)','Phosphate (PO-P)',...
            'Silicate (SiO4-Si)','Ph','Ph at 25 degrees and 0 dbar',...
            'Dissolved inorganic carbon','Total alkalinity','Dissolved organic carbon','Dissolved organic nitrogen',...
            'Total dissolved nitrogen','Chlorophyll-a'};
        osstandardname={'bottom_depth','sea_water_temperature','sea_water_practical_salinity','moles_of_oxygen_per_unit_mass_in_sea_water',...
            'moles_of_nitrate_per_unit_mass_in_sea_water','moles_of_nitrite_per_unit_mass_in_sea_water','moles_of_phosphate_per_unit_mass_in_sea_water',...
            'moles_of_silicate_per_unit_mass_in_sea_water','sea_water_ph_reported_on_total_scale','sea_water_ph_reported_on_total_scale_at_25_degrees_and_0_dbar',...
            'moles_of_dissolved_inorganic_carbon_per_unit_mass_in_sea_water','sea_water_alkalinity_per_unit_mass',...
            'moles_of_dissolved_organic_carbon_per_unit_mass_in_sea_water','moles_of_dissolved_organic_nitrogen_per_unit_mass_in_sea_water',...
            'moles_of_dissolved_total_nitrogen_per_unit_mass_in_sea_water','mass_concentration_of_chlorophyll_in_sea_water'};
        
        
    elseif inputtype == 'gridded'
        qcmanual='Lauvset et. al (2016), ESSD';
        globcomment=['Further details can be found in Lauvset et al. (2016).',...
            'A new global interior ocean mapped climatology: the 1°x1° GLODAP version 2, ',...
            'Earth Syst. Sci. Data, 8, 325–340, 2016, doi:10.5194/essd-8-325-2016'];
    end
end



%OStypes={'profile','time-series','trajectory', 'profile-trajectory'};
%sources={'research vessel','vessel of opportunity','moored surface buoy', 'drifting buoy'};
%sourceplatformcategorycodes={'31','32','41'};
%cdmtypes={'point','profile', 'section', 'station', 'station_profile', ...
%    'trajectory', 'grid', 'radial', 'swath', 'image'};
%filetitles={};

%platformcode=;
%platformname=;
%sitecode=;

