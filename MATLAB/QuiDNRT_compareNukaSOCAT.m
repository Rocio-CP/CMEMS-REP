clc
clear all
close all
workrootdir=...
    '/Users/rpr061/Documents/localtestarea/CARBON-REP-112020/';
 load([workrootdir, 'SOCATv2020_synthAE.mat']);
 
 % 2019 nuka data
 moreworkrootdir='/Users/rpr061/Documents/localtestarea/CARBON-REP-112020/latest/'
 file=dir([moreworkrootdir,'*OXYH2_20190*']);
 cd(moreworkrootdir)
 %
for a = 1:length(file)
time{a,1}=ncread(file(a).name, 'TIME');
fco2{a,1}=transpose(ncread(file(a).name, 'FCO2'));
fco2qc{a,1}=transpose(ncread(file(a).name, 'FCO2_QC'));
lat{a,1}=ncread(file(a).name, 'LATITUDE');
lon{a,1}=ncread(file(a).name, 'LONGITUDE');
end
fco2qc=cell2mat(fco2qc);

 
 %% SOCAT pre-2019 data
SOCAT.Expocode=cellstr(SOCAT.Expocode);
 for a=1:length(SOCAT.Expocode);
 platcodes(a,1)=strcmp(SOCAT.Expocode{a}(1:4),'26NA'); end
 %%
 %platcodes=platcodes';
 nukafco2=SOCAT.fCO2recuatm(platcodes & SOCAT.yr<2019);
 fakeyear=2019*ones(size(SOCAT.yr));

 %%
% time=cell2mat(time);
%time=time+712224;
% fco2=cell2mat(fco2);
figure
 hold on
 plot(nukatime, nukafco2, '.g');
 plot(time(fco2qc==2), fco2(fco2qc==2),'.r');
  plot(time(fco2qc==4), fco2(fco2qc==4),'.b');

xlabel('2019')
ylabel('fugacity of CO_2')
legend('SOCAT data','QCflag=2', 'QCflag=4')
  dateaxis('x')
datetick('x', 'mmm')
set(gca, 'ylim', [-10 700],...
    'xlim',[datenum('2019-04-01')-5, datenum('2019-06-05')+5]);
%%

 alltime=datenum(double([SOCAT.yr,...
     SOCAT.mon,...
     SOCAT.day1,...
     SOCAT.hh,...
     SOCAT.mm,...
     SOCAT.ss]));
 