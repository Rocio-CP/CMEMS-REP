clc
clear all
close all

%load('/Users/rpr061/Dropbox/BCDC_projects/CMEMS_INSTAC/Releases/03_2020/CARBON-REP-032020/original_files/GLODAPv2.2019_Merged_Master_File.mat')
load('/Users/rpr061/Downloads/GLODAPv2.2020_Merged_Master_File.mat');
%%
% For GLODAPv2.2019
v2plot={'temperature','salinity','oxygen','nitrate',...
    'nitrite','phosphate','silicate','phtsinsitutp','phts25p0','tco2',...
    'talk','doc','don','tdn','chla'};
osvars={'TEMP','PSAL','DOX2','NTAW',...
    'NTIW','PHOW','SLCW','PHPH','PH25','TICW',...
    'ALKW','CORG','NODW','NT1D','CPHL'};
% INSTAC long name
oslongname={'Sea temperature','Practical salinity',...
    'Dissolved oxygen','Nitrate (NO3-N)',...
    'Nitrite (NO2-N)','Phosphate (PO-P)','Silicate (SiO4-Si)',...
    'pH','Ph at 25 degrees and 0 dbar','Dissolved inorganic carbon',...
    'Total alkalinity','Dissolved organic carbon','Dissolved organic nitrogen',...
    'Total dissolved nitrogen','Chlorophyll-a'};

G2temperaturef=ones(size(G2temperature))*2;

% All sampling stations.
G2cruisestat=[G2cruise,G2station];
[C,IA,IC]=unique(G2cruisestat,'row');
stG2lat=G2latitude(IA);
stG2lon=G2longitude(IA);

stG2botdep=G2bottomdepth(IA); % 1 bottom depth per universal statin
h=histogram(stG2botdep, [0:100:7000,10000],...
        'FaceColor','none','Orientation','horizontal')%,'visible','off');
    bottomvals=(h.BinCounts);
    bottombins=(h.BinEdges);
    % "Cumulative" potential sampling depth,
    for h2=1:length(bottomvals); hyps(h2)=sum(bottomvals(h2:end));end

 %   for h2=1:length(h.BinCounts); hyps(h2)=sum(h.BinCounts(h2:end));end
    close

G2time=datenum([G2year,G2month,G2day,G2hour,G2minute,zeros(size(G2minute))]);

for vv=[8,9,12,13,14]%1:length(v2plot)
    variable=eval(['G2',v2plot{vv}]);
    v2filter=eval(['G2',v2plot{vv},'f==2']);
    
%     % Sampling map
     figure
 set(gcf,'Units','centimeters','Position',[0,0,20,20])

     subplot(2,2, [1,2])
    worldmap('world')
    hold on
    geoshow(stG2lat,stG2lon,...
        'displaytype','point','marker','o','markeredgecolor','k','markersize',5);
    geoshow(G2latitude(v2filter),G2longitude(v2filter),...
        'displaytype','point','marker','*',...
        'markeredgecolor','r','markersize',1);
    geoshow('landareas.shp')
    title ([oslongname{vv},' (',osvars{vv},')']);
    set(gca,'fontsize',24,'FontName','Times New Roman')
    tightmap
    
%     saveas(gcf,['/Users/rpr061/Desktop/',osvars{vv},'_map.png'])
    
    % Depth histogram    
    vardepths=G2depth(v2filter);
    [Z,Y,W]=histcounts(vardepths,[0:100:7000,10000]);
    depstat=[IC(v2filter),W]; % [stationindex, histogramindex]
    [V,U,T]=unique(depstat,'rows');
    
subplot(2,2,3)
    pos = get(gca, 'Position');

        pos(1) = 0.15;
    pos(3) = 0.27;
    pos(2) = 0.12;
    pos(4) = 0.43;
    
%     pos(1) = 0.15;
%     pos(3) = 0.3;
%     pos(2) = 0.12;
%     pos(4) = 0.43;
    set(gca, 'Position', pos)
%    set(gcf,'Units','centimeters','Position',[0,0,10,15])
    histogram('Bincounts',hyps,'Binedges',bottombins,...
        'FaceColor','none','Orientation','horizontal')
    hold on
    histogram(vardepths(U),[0:100:7000,10000],...
        'Orientation','horizontal',...
        'facealpha',0.5,'facecolor','blue');
        axis tight
    set(gca,'ydir','reverse','ylim',[0,6500],...
        'fontsize',20,'FontName','Times New Roman')
    ylabel('Depth (m)')
    xlabel('# Samples');
%    title ([oslongname{vv},' (',osvars{vv},')']);
%    saveas(gcf,['/Users/rpr061/Desktop/',osvars{vv},'_depth.png'])
    %clear hyps
    
%     % Time distribution
%     figure
%     set(gcf,'Units','centimeters','Position',[0,0,15,10]);
subplot(2,2,4)
    pos = get(gca, 'Position');

    pos(1) = 0.57;
    pos(3) = 0.37;
    pos(2) = 0.09;
    pos(4) = 0.43;
    
%     pos(1) = 0.55;
%     pos(3) = 0.4;
%     pos(2) = 0.09;
%     pos(4) = 0.43;
    set(gca, 'Position', pos)
set(gca,'fontsize',20,'FontName','Times New Roman');
    histogram(G2year(v2filter), [1970:1:2020], 'FaceColor',[0.2 0.2 0.5]);
%    title ([oslongname{vv},' (',osvars{vv},')']);
    set(gca,'fontsize',20,'FontName','Times New Roman',...
        'xlim',[1969,2021])
    ylabel('# Samples');
    xlabel('Year')
   axis tight
    
    
    saveas(gcf,['/Users/rpr061/Desktop/',osvars{vv},'.png']);
%     
    %close all
end

% %%
% plot(IC(10000:10500),G2depth(10000:10500), '*')
% hold on; plot(IC(10000:10500),G2bottomdepth(10000:10500),'o');
% %%
% depths=[0:100:1000];
% 
% for ii=1:max(IC);
% for iii=2:length(depths)
% if isempty(G2depth(G2depth(IC==ii)<depths(iii) & G2depth(IC==ii)>=depths(iii-1)));
%     nosurf(ii,iii)=1; end
% end
% end
%     