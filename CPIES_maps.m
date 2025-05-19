% CPIES_maps.m
%
% Make some maps and additional plots to assess which CPIES expected to
% yield the most benefit for A-0-A analysis
%

clear; close all

%% station geometry (from cruise reports)

clat=[-38.9074161;-38.7860271;-39.0689236;-39.0261897;-38.9316017;...
    -39.0312884;-38.8299051;-38.649488;-38.9671271;-39.128068;...
    -39.0270269;-38.7432545];
clon=[178.587845;178.670576;178.561354;178.509784;178.681372;178.650462;...
    178.872485;178.896993;178.880384;178.685767;178.840153;178.884469];
cdepth=[784;1051;1396;1444;1456;1585;2114;2257;2266;2324;2339;2436];
cname={'URI22-CF';'URI22-CD';'URI22-CK';'URI22-CA';'URI22-CG';'URI22-CJ';...
    'URI22-CE';'URI22-CB';'URI22-CH';'URI22-CL';'URI22-CI';'URI22-CC'};

%% network map

% addpath('m_map')
% 
% % setup base map
% figure(54); clf; hold on
% lat1=-39.5; lat2=-38.5;
% lon1=177.5; lon2=179.5;
% m_proj('mercator','longitudes',[lon1 lon2],'latitudes',[lat1 lat2]);
% 
% X=linspace(174.5,180,1320);
% Y=linspace(-37.5,-42,1080);
% Z=imread('../../hikurangi/gshhg-bin-2.3.6/NZ_bathymetry.tiff');
% Z(Z>=0)=NaN;
% [~,hc]=m_contourf(X,Y,Z,-5000:100:0,'linecolor','none');
% 
% m_gshhs_i('patch',[.7 .7 .7]);
% 
% m_contour(X,Y,Z,[-150 -150],'k','linewidth',1) % shelf at 150 m?
% m_text(177.95,-39.45,'150 m','fontsize',12,'rotation',35)
% m_contour(X,Y,Z,[-2750 -2750],'k','linewidth',1) % 2750 m does a decent job visually
% m_text(178.65,-39.45,'2750 m','fontsize',12,'rotation',90)
% 
% colormap(gca,cmocean('-deep'))
% cb1 = colorbar(gca,'eastoutside');
% ylabel(cb1,'Depth (m)')
% set(cb1,'fontsize',14)
% m_grid('xlabeldir','end','fontsize',14);
% 
% % station markers
% m_plot(clon,clat,'sk','linewidth',1,'markerfacecolor','r','markersize',10)
% m_text(clon-0.2,clat-0.02,cname,'fontsize',12)
% 
% fh=gcf;
% fh.PaperUnits='inches';
% fh.PaperPosition=[0 0 11 8.5];
% print('../figures/exploratory/CPIES/CPIES_map','-dpng','-r300')

%% pseudo across slope profile

% figure(55); clf; hold on
% plot(1:12,-cdepth,'ko-','markersize',10,'linewidth',1)
% set(gca,'fontsize',14)
% ylabel('Depth (m)')
% xlabel('~Distance from Shore')
% box on; grid on
% 
% fh=gcf;
% fh.PaperUnits='inches';
% fh.PaperPosition=[0 0 11 8.5];
% print('../figures/exploratory/CPIES/CPIES_profile','-dpng','-r300')

%% ECCO2 correlation map(s)

addpath('m_map')

% setup base map
figure(5); clf; hold on
lat1=-40.5; lat2=-37.5;
lon1=177.25; lon2=179.75;
m_proj('mercator','longitudes',[lon1 lon2],'latitudes',[lat1 lat2]);

E2=load('../ocean_data/ECCO2/ecco2_extraction_full.mat');
[E2.x,E2.y]=meshgrid(E2.lon,E2.lat);

std_mat=NaN(10,10,100);
n=0;
for i=1:length(E2.lon)
    for j=1:length(E2.lat)
        n=n+1;
        p1=squeeze(E2.bpa(i,j,:));
        reflon=E2.lon(i);
        reflat=E2.lat(j);
        if sum(p1)==0
            std_mat(:,:,n)=NaN;
            continue
        end
        for k=1:length(E2.lon)
            for l=1:length(E2.lat)
                p2=squeeze(E2.bpa(k,l,:));
                if sum(p2)==0
                    p2=NaN(size(p2));
                end
                std_mat(k,l,n)=std(p1-p2);
            end
        end
        % contour on a map
        [~,h1]=m_contourf(E2.x,E2.y,squeeze(std_mat(:,:,n))');
        h2=m_gshhs_i('patch',[.7 .7 .7]);
        h3=m_plot(reflon,reflat,'pk','markersize',10','markerfacecolor','y');
        h4=m_plot(clon,clat,'sk','linewidth',1,'markerfacecolor','r');
        m_grid('xlabeldir','end','fontsize',14);
        colorbar; clim([0 4])

        fh=gcf;
        fh.PaperUnits='inches';
        fh.PaperPosition=[0 0 8.5 11];
        print(['../figures/exploratory/CPIES/ECCO/map' num2str(n)],'-dpng','-r100')

        delete(h1); delete(h2); delete(h3); delete(h4);
    end
end

%% GLORYS correlation map(s)
% need to calculate pressure from model over broader region if I want to do
% this comparison (only have POBS-proximal pressures)

% addpath('m_map')
% 
% % setup base map
% figure(6); clf; hold on
% lat1=-40.5; lat2=-37.5;
% lon1=177.25; lon2=179.75;
% m_proj('mercator','longitudes',[lon1 lon2],'latitudes',[lat1 lat2]);
% 
% load('../ocean_data/GLORYS/glorys_pressure.mat');
% G=p_glo;
% [G.x,G.y]=meshgrid(G.lon,G.lat);
% 
% std_mat=NaN(10,10,100);
% n=0;
% for i=1:length(G.lon)
%     for j=1:length(G.lat)
%         n=n+1;
%         p1=squeeze(G.bpa(i,j,:));
%         reflon=G.lon(i);
%         reflat=G.lat(j);
%         if sum(p1)==0
%             std_mat(:,:,n)=NaN;
%             continue
%         end
%         for k=1:length(G.lon)
%             for l=1:length(G.lat)
%                 p2=squeeze(G.bpa(k,l,:));
%                 if sum(p2)==0
%                     p2=NaN(size(p2));
%                 end
%                 std_mat(k,l,n)=std(p1-p2);
%             end
%         end
%         % contour on a map
%         [~,h1]=m_contourf(G.x,G.y,squeeze(std_mat(:,:,n))');
%         h2=m_gshhs_i('patch',[.7 .7 .7]);
%         h3=m_plot(reflon,reflat,'pk','markersize',10','markerfacecolor','y');
%         h4=m_plot(clon,clat,'sk','linewidth',1,'markerfacecolor','r');
%         m_grid('xlabeldir','end','fontsize',14);
%         colorbar; clim([0 4])
% 
%         fh=gcf;
%         fh.PaperUnits='inches';
%         fh.PaperPosition=[0 0 8.5 11];
%         print(['../figures/exploratory/CPIES/ECCO/map' num2str(n)],'-dpng','-r100')
% 
%         delete(h1); delete(h2); delete(h3); delete(h4);
%     end
% end