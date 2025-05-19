% simple_maps.m
%
% Starter maps for seeing the HOBITSS pressure data distribution
%

clear; close all

addpath('m_map')

% set up projection
figure(54); clf; hold on
lat1=-40; lat2=-37.5;
lon1=176; lon2=180;
m_proj('mercator','longitudes',[lon1 lon2],'latitudes',[lat1 lat2]);

% bathymetry
X=linspace(174.5,180,1320);
Y=linspace(-37.5,-42,1080);
Z=imread('../gshhg-bin-2.3.6/NZ_bathymetry.tiff');
Z(Z>=0)=NaN;

% coastline
m_gshhs_i('patch',[.7 .7 .7]);

% contour shelf and slope
m_contour(X,Y,Z,[-150 -150],'k','linewidth',1) % shelf at 150 m?
m_text(177.05,-39.85,'150 m','fontsize',12,'rotation',55)
% I think the slope base is ~3000 m in the north and ~2500 m in the south
m_contour(X,Y,Z,[-2750 -2750],'k','linewidth',1) % 2750 m does a decent job visually
m_text(178.4,-39.9,'2750 m','fontsize',12,'rotation',80)

% format axes
colormap(gca,cmocean('-deep'))
cb1 = colorbar(gca,'eastoutside');
ylabel(cb1,'Depth (m)')
set(cb1,'fontsize',14)

m_grid('xlabeldir','end','fontsize',14);

%--- instruments
% GPS stations
gps=load('../processed_data/GPS.mat');
cond_gps= gps.stalon>=lon1 & gps.stalon<=lon2 & gps.stalat>=lat1 & gps.stalat<=lat2;
h_gps=m_scatter(gps.stalon(cond_gps),gps.stalat(cond_gps),100,'ys','filled','markeredgecolor','k','linewidth',1);
h_gpst=m_text(gps.stalon(cond_gps)-0.15,gps.stalat(cond_gps)+0.1,gps.staname(cond_gps),'fontsize',12);

% APG stations
yr_str={'I','III','IV','VI','VIII'};
for i=1:length(yr_str)
    load(['../processed_data/HOBITSS_' yr_str{i} '.mat'])

    h1=m_scatter(stalon,stalat,200,'c^','filled','markeredgecolor','k','linewidth',1);
    h1t=m_text(stalon+0.05,stalat+0.1,staname,'fontsize',12);

    if i>3
        cork=load('../processed_data/CORK.mat');
        hck=m_scatter(cork.stalon,cork.stalat,200,'ro','filled','markeredgecolor','k','linewidth',1);
        hckt=m_text(cork.stalon+0.05,cork.stalat+0.1,cork.staname,'fontsize',12);

        hl=legend([h1,hck,h_gps],'APG','CORK','GPS','location','northeast');
    elseif i>2
        cork=load('../processed_data/CORK.mat');
        hck=m_scatter(cork.stalon(2:3),cork.stalat(2:3),200,'ro','filled','markeredgecolor','k','linewidth',1);
        hckt=m_text(cork.stalon(2:3)+0.05,cork.stalat(2:3)+0.1,cork.staname(2:3),'fontsize',12);
        hl=legend([h1,hck,h_gps],'APG','CORK','GPS','location','northeast');
    else
        hl=legend([h1,h_gps],'APG','GPS','location','northeast');
    end
    set(hl,'fontsize',14)

    % save
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 11 8.5];
    print(['../figures/exploratory/maps/HOBBITS_' yr_str{i}],'-dtiff','-r300')
    print(['../figures/exploratory/maps/HOBBITS_' yr_str{i}],'-depsc','-vector')

    % remove APG markers
    delete(h1)
    delete(h1t)
    delete(hl)
    if i>2
        delete(hck)
        delete(hckt)
    end
end

%% HOBITSS V includes stations off Hawke's Bay

% set up projection
figure(55); clf; hold on
lat1=-40.5; lat2=-38;
lon1=176; lon2=180;
m_proj('mercator','longitudes',[lon1 lon2],'latitudes',[lat1 lat2]);

% bathymetry
[~,hc]=m_contourf(X,Y,Z,-5000:100:0,'linecolor','none');

% coastline
m_gshhs_i('patch',[.7 .7 .7]);

% contour shelf and slope
m_contour(X,Y,Z,[-150 -150],'k','linewidth',1) % shelf at 150 m?
m_text(176.9,-40.25,'150 m','fontsize',12,'rotation',60)
% I think the slope base is ~3000 m in the north and ~2500 m in the south
m_contour(X,Y,Z,[-2750 -2750],'k','linewidth',1) % 2750 m does a decent job visually
m_text(178,-40.45,'2750 m','fontsize',12,'rotation',45)

% format axes
colormap(gca,cmocean('-deep'))
cb1 = colorbar(gca,'eastoutside');
ylabel(cb1,'Depth (m)')
set(cb1,'fontsize',14)

m_grid('xlabeldir','end','fontsize',14);

%--- instruments
% GPS stations
gps=load('../processed_data/GPS.mat');
cond_gps= gps.stalon>=lon1 & gps.stalon<=lon2 & gps.stalat>=lat1 & gps.stalat<=lat2;
h_gps=m_scatter(gps.stalon(cond_gps),gps.stalat(cond_gps),100,'ys','filled','markeredgecolor','k','linewidth',1);
h_gpst=m_text(gps.stalon(cond_gps)-0.15,gps.stalat(cond_gps)+0.1,gps.staname(cond_gps),'fontsize',12);

% CORK line
cork=load('../processed_data/CORK.mat');
hck=m_scatter(cork.stalon,cork.stalat,200,'ro','filled','markeredgecolor','k','linewidth',1);
hckt=m_text(cork.stalon+0.05,cork.stalat+0.1,cork.staname,'fontsize',12);

% APG stations
load('../processed_data/HOBITSS_V.mat')

h2=m_scatter(stalon,stalat,200,'c^','filled','markeredgecolor','k','linewidth',1);
h2t=m_text(stalon+0.05,stalat+0.1,staname,'fontsize',12);

hl=legend([h2,hck,h_gps],'APG','CORK','GPS','location','northeast');
set(hl,'fontsize',14)

% save
fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../figures/exploratory/maps/HOBBITS_V','-dtiff','-r300')

%% HOBITSS VII is much further south

% set up projection
figure(56); clf; hold on
lat1=-41; lat2=-38.5;
lon1=176; lon2=180;
m_proj('mercator','longitudes',[lon1 lon2],'latitudes',[lat1 lat2]);

% bathymetry
[~,hc]=m_contourf(X,Y,Z,-5000:100:0,'linecolor','none');

% coastline
m_gshhs_i('patch',[.7 .7 .7]);

% contour shelf and slope
m_contour(X,Y,Z,[-150 -150],'k','linewidth',1) % shelf at 150 m?
m_text(176.5,-40.8,'150 m','fontsize',12,'rotation',45)
% I think the slope base is ~3000 m in the north and ~2500 m in the south
m_contour(X,Y,Z,[-2750 -2750],'k','linewidth',1) % 2750 m does a decent job visually
m_text(177.8,-40.95,'2750 m','fontsize',12,'rotation',45)

% format axes
colormap(gca,cmocean('-deep'))
cb1 = colorbar(gca,'eastoutside');
ylabel(cb1,'Depth (m)')
set(cb1,'fontsize',14)

m_grid('xlabeldir','end','fontsize',14);

%--- instruments
% CORK line
cork=load('../processed_data/CORK.mat');
hck=m_scatter(cork.stalon,cork.stalat,200,'ro','filled','markeredgecolor','k','linewidth',1);
hckt=m_text(cork.stalon+0.05,cork.stalat+0.1,cork.staname,'fontsize',12);

% GPS stations
gps=load('../processed_data/GPS.mat');
cond_gps= gps.stalon>=lon1 & gps.stalon<=lon2 & gps.stalat>=lat1 & gps.stalat<=lat2;
h_gps=m_scatter(gps.stalon(cond_gps),gps.stalat(cond_gps),100,'ys','filled','markeredgecolor','k','linewidth',1);
h_gpst=m_text(gps.stalon(cond_gps)-0.15,gps.stalat(cond_gps)+0.1,gps.staname(cond_gps),'fontsize',12);

% APG stations
load('../processed_data/HOBITSS_VII.mat')

h3=m_scatter(stalon,stalat,200,'c^','filled','markeredgecolor','k','linewidth',1);
h3t=m_text(stalon+0.05,stalat+0.1,staname,'fontsize',12);

hl=legend([h3,hck,h_gps],'APG','CORK','GPS','location','southeast');
set(hl,'fontsize',14)

% save
fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../figures/exploratory/maps/HOBBITS_VII','-dtiff','-r300')