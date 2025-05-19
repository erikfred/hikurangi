% import_GPS.m
%

clear; close all

%% GPS time series

% station info
staname={'ANAU','CNST','MAKO','PARI','MAHI','CKID','KAHU','PAWA','AKTO',...
    'CAST','OROA','PORA','BIRF'};
stalat=[-38.2682,-38.4880,-38.6438,-38.9226,-39.1526,-39.6579,-39.7938,-40.0331,-40.5398,...
    -40.9098,-40.1044,-40.2664,-40.6798];
stalon=[178.2912,178.2111,178.1291,177.8833,177.907,177.0764,176.8639,176.8763,176.4612,...
    176.2016,176.6807,176.6352,176.2461];
stadepth=[0,0,0,0,0,0,0,0,0,...
    0,0,0,0];

% ANAU
fid = '../apg_data/gps/ANAU_e.txt';
ANAU = readtable(fid);
temp = table2array(ANAU(:,1));
temp1 = cell2mat(temp);
temp2 = datenum(temp1(:,1:10),'yyyy-mm-dd') + datenum(temp1(:,12:19),'HH:MM:SS')-datenum(2023,1,1);
ANAU = [temp2,table2array(ANAU(:,2:3))]';

% CNST
fid = '../apg_data/gps/CNST_e.txt';
CNST = readtable(fid);
temp = table2array(CNST(:,1));
temp1 = cell2mat(temp);
temp2 = datenum(temp1(:,1:10),'yyyy-mm-dd') + datenum(temp1(:,12:19),'HH:MM:SS')-datenum(2023,1,1);
CNST = [temp2,table2array(CNST(:,2:3))]';

% MAKO
fid = '../apg_data/gps/MAKO_e.txt';
MAKO = readtable(fid);
temp = table2array(MAKO(:,1));
temp1 = cell2mat(temp);
temp2 = datenum(temp1(:,1:10),'yyyy-mm-dd') + datenum(temp1(:,12:19),'HH:MM:SS')-datenum(2023,1,1);
MAKO = [temp2,table2array(MAKO(:,2:3))]';

% PARI
fid = '../apg_data/gps/PARI_e.txt';
PARI = readtable(fid);
temp = table2array(PARI(:,1));
temp1 = cell2mat(temp);
temp2 = datenum(temp1(:,1:10),'yyyy-mm-dd') + datenum(temp1(:,12:19),'HH:MM:SS')-datenum(2023,1,1);
PARI = [temp2,table2array(PARI(:,2:3))]';

% MAHI
fid = '../apg_data/gps/MAHI_e.txt';
MAHI = readtable(fid);
temp = table2array(MAHI(:,1));
temp1 = cell2mat(temp);
temp2 = datenum(temp1(:,1:10),'yyyy-mm-dd') + datenum(temp1(:,12:19),'HH:MM:SS')-datenum(2023,1,1);
MAHI = [temp2,table2array(MAHI(:,2:3))]';

% CKID
fid = '../apg_data/gps/CKID_e.txt';
CKID = readtable(fid);
temp = table2array(CKID(:,1));
temp1 = cell2mat(temp);
temp2 = datenum(temp1(:,1:10),'yyyy-mm-dd') + datenum(temp1(:,12:19),'HH:MM:SS')-datenum(2023,1,1);
CKID = [temp2,table2array(CKID(:,2:3))]';

% KAHU
fid = '../apg_data/gps/CAST_e.txt';
KAHU = readtable(fid);
temp = table2array(KAHU(:,1));
temp1 = cell2mat(temp);
temp2 = datenum(temp1(:,1:10),'yyyy-mm-dd') + datenum(temp1(:,12:19),'HH:MM:SS')-datenum(2023,1,1);
KAHU = [temp2,table2array(KAHU(:,2:3))]';

% PAWA
fid = '../apg_data/gps/PAWA_e.txt';
PAWA = readtable(fid);
temp = table2array(PAWA(:,1));
temp1 = cell2mat(temp);
temp2 = datenum(temp1(:,1:10),'yyyy-mm-dd') + datenum(temp1(:,12:19),'HH:MM:SS')-datenum(2023,1,1);
PAWA = [temp2,table2array(PAWA(:,2:3))]';

% AKTO
fid = '../apg_data/gps/AKTO_e.txt';
AKTO = readtable(fid);
temp = table2array(AKTO(:,1));
temp1 = cell2mat(temp);
temp2 = datenum(temp1(:,1:10),'yyyy-mm-dd') + datenum(temp1(:,12:19),'HH:MM:SS')-datenum(2023,1,1);
AKTO = [temp2,table2array(AKTO(:,2:3))]';

% CAST
fid = '../apg_data/gps/CAST_e.txt';
CAST = readtable(fid);
temp = table2array(CAST(:,1));
temp1 = cell2mat(temp);
temp2 = datenum(temp1(:,1:10),'yyyy-mm-dd') + datenum(temp1(:,12:19),'HH:MM:SS')-datenum(2023,1,1);
CAST = [temp2,table2array(CAST(:,2:3))]';

% OROA
fid = '../apg_data/gps/CAST_e.txt';
OROA = readtable(fid);
temp = table2array(OROA(:,1));
temp1 = cell2mat(temp);
temp2 = datenum(temp1(:,1:10),'yyyy-mm-dd') + datenum(temp1(:,12:19),'HH:MM:SS')-datenum(2023,1,1);
OROA = [temp2,table2array(OROA(:,2:3))]';

% PORA
fid = '../apg_data/gps/CAST_e.txt';
PORA = readtable(fid);
temp = table2array(PORA(:,1));
temp1 = cell2mat(temp);
temp2 = datenum(temp1(:,1:10),'yyyy-mm-dd') + datenum(temp1(:,12:19),'HH:MM:SS')-datenum(2023,1,1);
PORA = [temp2,table2array(PORA(:,2:3))]';

% BIRF
fid = '../apg_data/gps/CAST_e.txt';
BIRF = readtable(fid);
temp = table2array(BIRF(:,1));
temp1 = cell2mat(temp);
temp2 = datenum(temp1(:,1:10),'yyyy-mm-dd') + datenum(temp1(:,12:19),'HH:MM:SS')-datenum(2023,1,1);
BIRF = [temp2,table2array(BIRF(:,2:3))]';

%-----PLOTTING-----%
%--- make adjustments as needed for targeted figs
% stations north of Hawke's Bay
figure(2); clf; hold on
plot(ANAU(1,ANAU(1,:)>=datenum(2021,11,01) & ANAU(1,:)<=datenum(2022,10,01)),...
    detrend(ANAU(2,ANAU(1,:)>=datenum(2021,11,01) & ANAU(1,:)<=datenum(2022,10,01))),...
    'linewidth',2)
plot(CNST(1,CNST(1,:)>=datenum(2021,11,01) & CNST(1,:)<=datenum(2022,10,01)),...
    detrend(CNST(2,CNST(1,:)>=datenum(2021,11,01) & CNST(1,:)<=datenum(2022,10,01)))-20,...
    'linewidth',2)
plot(MAKO(1,MAKO(1,:)>=datenum(2021,11,01) & MAKO(1,:)<=datenum(2022,10,01)),...
    detrend(MAKO(2,MAKO(1,:)>=datenum(2021,11,01) & MAKO(1,:)<=datenum(2022,10,01)))-40,...
    'linewidth',2)
plot(PARI(1,PARI(1,:)>=datenum(2021,11,01) & PARI(1,:)<=datenum(2022,10,01)),...
    detrend(PARI(2,PARI(1,:)>=datenum(2021,11,01) & PARI(1,:)<=datenum(2022,10,01)))-60,...
    'linewidth',2)
plot(MAHI(1,MAHI(1,:)>=datenum(2021,11,01) & MAHI(1,:)<=datenum(2022,10,01)),...
    detrend(MAHI(2,MAHI(1,:)>=datenum(2021,11,01) & MAHI(1,:)<=datenum(2022,10,01)))-80,...
    'linewidth',2)
plot(CKID(1,CKID(1,:)>=datenum(2021,11,01) & CKID(1,:)<=datenum(2022,10,01)),...
    detrend(CKID(2,CKID(1,:)>=datenum(2021,11,01) & CKID(1,:)<=datenum(2022,10,01)))-100,...
    'linewidth',2)
plot(KAHU(1,KAHU(1,:)>=datenum(2021,11,01) & KAHU(1,:)<=datenum(2022,10,01)),...
    detrend(KAHU(2,KAHU(1,:)>=datenum(2021,11,01) & KAHU(1,:)<=datenum(2022,10,01)))-120,...
    'linewidth',2)
% plot(PAWA(1,PAWA(1,:)>=datenum(2021,11,01) & PAWA(1,:)<=datenum(2022,10,01)),...
%     detrend(PAWA(2,PAWA(1,:)>=datenum(2021,11,01) & PAWA(1,:)<=datenum(2022,10,01)))-140,...
%     'linewidth',2)
% plot(AKTO(1,AKTO(1,:)>=datenum(2021,11,01) & AKTO(1,:)<=datenum(2022,10,01)),...
%     detrend(AKTO(2,AKTO(1,:)>=datenum(2021,11,01) & AKTO(1,:)<=datenum(2022,10,01)))-160,...
%     'linewidth',2)
% plot(CAST(1,CAST(1,:)>=datenum(2021,11,01) & CAST(1,:)<=datenum(2022,10,01)),...
%     detrend(CAST(2,CAST(1,:)>=datenum(2021,11,01) & CAST(1,:)<=datenum(2022,10,01)))-180,...
%     'linewidth',2)
% plot(OROA(1,OROA(1,:)>=datenum(2021,11,01) & OROA(1,:)<=datenum(2022,10,01)),...
%     detrend(OROA(2,OROA(1,:)>=datenum(2021,11,01) & OROA(1,:)<=datenum(2022,10,01)))-200,...
%     'linewidth',2)
% plot(PORA(1,PORA(1,:)>=datenum(2021,11,01) & PORA(1,:)<=datenum(2022,10,01)),...
%     detrend(PORA(2,PORA(1,:)>=datenum(2021,11,01) & PORA(1,:)<=datenum(2022,10,01)))-220,...
%     'linewidth',2)
% plot(BIRF(1,BIRF(1,:)>=datenum(2021,11,01) & BIRF(1,:)<=datenum(2022,10,01)),...
%     detrend(BIRF(2,BIRF(1,:)>=datenum(2021,11,01) & BIRF(1,:)<=datenum(2022,10,01)))-240,...
%     'linewidth',2)
legend('ANAU','CNST','MAKO','PARI','MAHI','CKID','KAHU','location','northwest')
ylim([-140 20])
datetick('x',3)
ylabel('\DeltaP (cm)')
set(gca,'fontsize',14)
box on; grid on

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print('../figures/exploratory/HOBITSS_VIII/differences/coastal_GPS_21-22_N','-dpng','-r300')

% stations south of Hawke's Bay
figure(3); clf; hold on
% plot(ANAU(1,ANAU(1,:)>=datenum(2021,11,01) & ANAU(1,:)<=datenum(2022,10,01)),...
%     detrend(ANAU(2,ANAU(1,:)>=datenum(2021,11,01) & ANAU(1,:)<=datenum(2022,10,01))),...
%     'linewidth',2)
% plot(CNST(1,CNST(1,:)>=datenum(2021,11,01) & CNST(1,:)<=datenum(2022,10,01)),...
%     detrend(CNST(2,CNST(1,:)>=datenum(2021,11,01) & CNST(1,:)<=datenum(2022,10,01)))-20,...
%     'linewidth',2)
% plot(MAKO(1,MAKO(1,:)>=datenum(2021,11,01) & MAKO(1,:)<=datenum(2022,10,01)),...
%     detrend(MAKO(2,MAKO(1,:)>=datenum(2021,11,01) & MAKO(1,:)<=datenum(2022,10,01)))-40,...
%     'linewidth',2)
% plot(PARI(1,PARI(1,:)>=datenum(2021,11,01) & PARI(1,:)<=datenum(2022,10,01)),...
%     detrend(PARI(2,PARI(1,:)>=datenum(2021,11,01) & PARI(1,:)<=datenum(2022,10,01)))-60,...
%     'linewidth',2)
% plot(MAHI(1,MAHI(1,:)>=datenum(2021,11,01) & MAHI(1,:)<=datenum(2022,10,01)),...
%     detrend(MAHI(2,MAHI(1,:)>=datenum(2021,11,01) & MAHI(1,:)<=datenum(2022,10,01)))-80,...
%     'linewidth',2)
% plot(CKID(1,CKID(1,:)>=datenum(2021,11,01) & CKID(1,:)<=datenum(2022,10,01)),...
%     detrend(CKID(2,CKID(1,:)>=datenum(2021,11,01) & CKID(1,:)<=datenum(2022,10,01)))-100,...
%     'linewidth',2)
% plot(KAHU(1,KAHU(1,:)>=datenum(2021,11,01) & KAHU(1,:)<=datenum(2022,10,01)),...
%     detrend(KAHU(2,KAHU(1,:)>=datenum(2021,11,01) & KAHU(1,:)<=datenum(2022,10,01)))-120,...
%     'linewidth',2)
plot(PAWA(1,PAWA(1,:)>=datenum(2021,11,01) & PAWA(1,:)<=datenum(2022,10,01)),...
    detrend(PAWA(2,PAWA(1,:)>=datenum(2021,11,01) & PAWA(1,:)<=datenum(2022,10,01)))-140,...
    'linewidth',2)
plot(AKTO(1,AKTO(1,:)>=datenum(2021,11,01) & AKTO(1,:)<=datenum(2022,10,01)),...
    detrend(AKTO(2,AKTO(1,:)>=datenum(2021,11,01) & AKTO(1,:)<=datenum(2022,10,01)))-160,...
    'linewidth',2)
plot(CAST(1,CAST(1,:)>=datenum(2021,11,01) & CAST(1,:)<=datenum(2022,10,01)),...
    detrend(CAST(2,CAST(1,:)>=datenum(2021,11,01) & CAST(1,:)<=datenum(2022,10,01)))-180,...
    'linewidth',2)
plot(OROA(1,OROA(1,:)>=datenum(2021,11,01) & OROA(1,:)<=datenum(2022,10,01)),...
    detrend(OROA(2,OROA(1,:)>=datenum(2021,11,01) & OROA(1,:)<=datenum(2022,10,01)))-200,...
    'linewidth',2)
plot(PORA(1,PORA(1,:)>=datenum(2021,11,01) & PORA(1,:)<=datenum(2022,10,01)),...
    detrend(PORA(2,PORA(1,:)>=datenum(2021,11,01) & PORA(1,:)<=datenum(2022,10,01)))-220,...
    'linewidth',2)
plot(BIRF(1,BIRF(1,:)>=datenum(2021,11,01) & BIRF(1,:)<=datenum(2022,10,01)),...
    detrend(BIRF(2,BIRF(1,:)>=datenum(2021,11,01) & BIRF(1,:)<=datenum(2022,10,01)))-240,...
    'linewidth',2)
legend('PAWA','AKTO','CAST','OROA','PORA','BIRF','location','northwest')
ylim([-260 -120])
datetick('x',3)
ylabel('\DeltaP (cm)')
set(gca,'fontsize',14)
box on; grid on

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print('../figures/exploratory/HOBITSS_VIII/differences/coastal_GPS_21-22_S','-dpng','-r300')

save('../processed_data/GPS','ANAU','CNST','MAKO','PARI','MAHI','CKID',...
    'PAWA','AKTO','CAST','KAHU','OROA','PORA','BIRF','staname','stalat',...
    'stalon','stadepth')