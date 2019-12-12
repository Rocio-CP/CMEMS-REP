%the horrific QC setting

% vars={'bath','temp','sal','oxy','nitra',...
%     'nitri','phos','sili','ph','ph25','dic',...
%     'alk','doc','don','tdn','chla'};


netcdf.putAtt(ncid, bathvarqc, 'valid_min', int8(0));
netcdf.putAtt(ncid, bathvarqc, 'valid_max', int8(9));
netcdf.putAtt(ncid, bathvarqc, 'flag_values',  [int8(0)  int8(1)  int8(2)  int8(3)  int8(4)  int8(5)  int8(6)  int8(7)  int8(8)  int8(9)]);

netcdf.putAtt(ncid, tempvarqc, 'valid_min', int8(0));
netcdf.putAtt(ncid, tempvarqc, 'valid_max', int8(9));
netcdf.putAtt(ncid, tempvarqc, 'flag_values',  [int8(0)  int8(1)  int8(2)  int8(3)  int8(4)  int8(5)  int8(6)  int8(7)  int8(8)  int8(9)]);

netcdf.putAtt(ncid, salvarqc, 'valid_min', int8(0));
netcdf.putAtt(ncid, salvarqc, 'valid_max', int8(9));
netcdf.putAtt(ncid, salvarqc, 'flag_values',  [int8(0)  int8(1)  int8(2)  int8(3)  int8(4)  int8(5)  int8(6)  int8(7)  int8(8)  int8(9)]);

netcdf.putAtt(ncid, oxyvarqc, 'valid_min', int8(0));
netcdf.putAtt(ncid, oxyvarqc, 'valid_max', int8(9));
netcdf.putAtt(ncid, oxyvarqc, 'flag_values',  [int8(0)  int8(1)  int8(2)  int8(3)  int8(4)  int8(5)  int8(6)  int8(7)  int8(8)  int8(9)]);

netcdf.putAtt(ncid, nitravarqc, 'valid_min', int8(0));
netcdf.putAtt(ncid, nitravarqc, 'valid_max', int8(9));
netcdf.putAtt(ncid, nitravarqc, 'flag_values',  [int8(0)  int8(1)  int8(2)  int8(3)  int8(4)  int8(5)  int8(6)  int8(7)  int8(8)  int8(9)]);

netcdf.putAtt(ncid, nitrivarqc, 'valid_min', int8(0));
netcdf.putAtt(ncid, nitrivarqc, 'valid_max', int8(9));
netcdf.putAtt(ncid, nitrivarqc, 'flag_values',  [int8(0)  int8(1)  int8(2)  int8(3)  int8(4)  int8(5)  int8(6)  int8(7)  int8(8)  int8(9)]);

netcdf.putAtt(ncid, phosvarqc, 'valid_min', int8(0));
netcdf.putAtt(ncid, phosvarqc, 'valid_max', int8(9));
netcdf.putAtt(ncid, phosvarqc, 'flag_values',  [int8(0)  int8(1)  int8(2)  int8(3)  int8(4)  int8(5)  int8(6)  int8(7)  int8(8)  int8(9)]);

netcdf.putAtt(ncid, silivarqc, 'valid_min', int8(0));
netcdf.putAtt(ncid, silivarqc, 'valid_max', int8(9));
netcdf.putAtt(ncid, silivarqc, 'flag_values',  [int8(0)  int8(1)  int8(2)  int8(3)  int8(4)  int8(5)  int8(6)  int8(7)  int8(8)  int8(9)]);

netcdf.putAtt(ncid, phvarqc, 'valid_min', int8(0));
netcdf.putAtt(ncid, phvarqc, 'valid_max', int8(9));
netcdf.putAtt(ncid, phvarqc, 'flag_values',  [int8(0)  int8(1)  int8(2)  int8(3)  int8(4)  int8(5)  int8(6)  int8(7)  int8(8)  int8(9)]);

netcdf.putAtt(ncid, ph25varqc, 'valid_min', int8(0));
netcdf.putAtt(ncid, ph25varqc, 'valid_max', int8(9));
netcdf.putAtt(ncid, ph25varqc, 'flag_values',  [int8(0)  int8(1)  int8(2)  int8(3)  int8(4)  int8(5)  int8(6)  int8(7)  int8(8)  int8(9)]);

netcdf.putAtt(ncid, dicvarqc, 'valid_min', int8(0));
netcdf.putAtt(ncid, dicvarqc, 'valid_max', int8(9));
netcdf.putAtt(ncid, dicvarqc, 'flag_values',  [int8(0)  int8(1)  int8(2)  int8(3)  int8(4)  int8(5)  int8(6)  int8(7)  int8(8)  int8(9)]);

netcdf.putAtt(ncid, alkvarqc, 'valid_min', int8(0));
netcdf.putAtt(ncid, alkvarqc, 'valid_max', int8(9));
netcdf.putAtt(ncid, alkvarqc, 'flag_values',  [int8(0)  int8(1)  int8(2)  int8(3)  int8(4)  int8(5)  int8(6)  int8(7)  int8(8)  int8(9)]);

netcdf.putAtt(ncid, docvarqc, 'valid_min', int8(0));
netcdf.putAtt(ncid, docvarqc, 'valid_max', int8(9));
netcdf.putAtt(ncid, docvarqc, 'flag_values',  [int8(0)  int8(1)  int8(2)  int8(3)  int8(4)  int8(5)  int8(6)  int8(7)  int8(8)  int8(9)]);

netcdf.putAtt(ncid, donvarqc, 'valid_min', int8(0));
netcdf.putAtt(ncid, donvarqc, 'valid_max', int8(9));
netcdf.putAtt(ncid, donvarqc, 'flag_values',  [int8(0)  int8(1)  int8(2)  int8(3)  int8(4)  int8(5)  int8(6)  int8(7)  int8(8)  int8(9)]);

netcdf.putAtt(ncid, tdnvarqc, 'valid_min', int8(0));
netcdf.putAtt(ncid, tdnvarqc, 'valid_max', int8(9));
netcdf.putAtt(ncid, tdnvarqc, 'flag_values',  [int8(0)  int8(1)  int8(2)  int8(3)  int8(4)  int8(5)  int8(6)  int8(7)  int8(8)  int8(9)]);

netcdf.putAtt(ncid, chlavarqc, 'valid_min', int8(0));
netcdf.putAtt(ncid, chlavarqc, 'valid_max', int8(9));
netcdf.putAtt(ncid, chlavarqc, 'flag_values',  [int8(0)  int8(1)  int8(2)  int8(3)  int8(4)  int8(5)  int8(6)  int8(7)  int8(8)  int8(9)]);
