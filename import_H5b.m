% import_H5b.m
%
% Read HOBITSS V (2018-2019) data from Katie's public files
% (https://doi.org/10.5281/zenodo.5834879)and reformat for consistency 
% with other years
% Note that 'import_H5.m' instead imports the data directly from txt/csv
%

clear; close all

%% HOBITSS V (2018-2019)

topdir='/Volumes/Gorgoroth/apg_database/original/2018-2019_HOBITSS-V/';
wrtdir='/Volumes/Gorgoroth/apg_database/processed/2018-2019_HOBITSS-Vb/';
figdir='../figures/exploratory/HOBITSS_V/hobitss5b_stack';
svdir='../processed_data/HOBITSS_Vb';

load([topdir 'APG_hikurangi_2019.mat'])

% station info
staname={'GNS18-1','GNS18-3','GNS18-7','LBPR18-4','LBPR18-5','KU18-1','KU18-2',...
    'KU18-3','KU18-4','POBS18-1','POBS18-2','POBS18-3','U1518','U1519'};
stalat=[-39.876,-39.647,-40.283,-39.767,-39.693,-38.906,-38.876,...
    -38.893,-38.794,-39.852,-39.892,-39.919,-38.859,-38.727];
stalon=[178.666,178.149,177.523,178.485,178.692,178.982,178.845,...
    178.755,178.671,177.918,178.088,178.352,178.896,178.615];
stadepth=[3309,779,1895,2207,3311,3483,1912,...
    1357,1083,1225,1455,1930,2631,1000];

altname={seafloor_pressure.site};

tf={};
pf={};
for i=1:length(staname)
    ii=find(strcmp(staname(i),altname));
    tf{i}=datenum(seafloor_pressure(ii).dates);
    pf{i}=seafloor_pressure(ii).pressure;
end

% get dimensions consistent
tf{10}=tf{10}';
tf{11}=tf{11}';
tf{12}=tf{12}';
pf{10}=pf{10}';
pf{11}=pf{11}';
pf{12}=pf{12}';

%-----PLOTTING-----%
% sort by depth
[depth,id]=sort(stadepth,'descend');
tf=tf(id);
pf=pf(id);
sta=staname(id);

figure(2); clf; hold on
for j=1:length(pf)
    t_plot = tf{j};
    p_plot = detrend(pf{j});
    plot(t_plot,p_plot+(j-1)*20,'linewidth',1)
    text(t_plot(end)+10,p_plot(end-100)+(j-1)*20,{sta{j};[num2str(depth(j)) ' m']},'fontsize',14)
end
xlim([datenum(2018,12,01) datenum(2019,12,01)])
datetick('x','keeplimits')
set(gca,'fontsize',14)
ylabel('\DeltaP (hPa)')
title('HOBITSS 2018-2019')
box on; grid on

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print(figdir,'-dpng','-r300')

% save data
save(svdir,'tf','pf','staname','stalat','stalon','stadepth')

% write to text files
for i=1:length(tf)
    % decimated, detided
    writematrix([tf{i}',pf{i}'],[wrtdir staname{i} '_1hr_detided'])
end