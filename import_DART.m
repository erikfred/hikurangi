% import_DART.m
%

clear; close all

%% DART time series

topdir='/Volumes/Gorgoroth/apg_database/original/DART/';
wrtdir='/Volumes/Gorgoroth/apg_database/processed/DART/';
figdir='../figures/exploratory/DART/dart_stack';
svdir='../processed_data/DART';

% station info
staname={'NZA','NZB','NZC','NZD','NZE','NZF','NZG','NZH','NZI'};
stalat=[-42.3690,-40.5979,-38.1969,-36.1000,-36.0500,-29.6826,-23.3517,-20.0885,-16.8890];
stalon=[176.9120,179.1005,179.7968,178.6009,177.6986,175.0125,173.4018,171.8630,171.1905];
stadepth=[2600,3200,3500,2500,5800,5100,5700,5500,5200];

%--- NZA
fid = [topdir 'NZA_water_height.csv'];
NZA = readtable(fid);
temp = table2array(NZA(:,7));
temp1 = cell2mat(temp);
temp2 = datenum([temp1(:,1:10),repmat(' ',length(temp1),1),temp1(:,12:19)]);
nza.t = temp2;
nza.p = table2array(NZA(:,8));
disp(median(nza.p))

% empirically-identified offsets (instrument changes?)
nza.i_offset=[69652,133231];

% fit offset(s) with linear model
tinv=nza.t-nza.t(1);
for i=1:length(nza.i_offset)
    io=nza.i_offset(i);
    G=[[tinv(1:io);zeros(size(tinv(io+1:end)))],[ones(size(tinv(1:io)));zeros(size(tinv(io+1:end)))],...
        [zeros(size(tinv(1:io)));tinv(io+1:end)],[zeros(size(tinv(1:io)));ones(size(tinv(io+1:end)))]];
    m=inv(G'*G)*G'*nza.p;
    fit_test=G*m;
    nza.p=nza.p-fit_test;
end

% decimation and tidal filter
[nza.tf,nza.pf,~,~] = downsample_uneven(nza.t,nza.p,1/24);
nza.pf=Z_godin(nza.pf)*100; % [hPa]
nza.tf(isnan(nza.pf))=[];
nza.pf(isnan(nza.pf))=[];

% interpolate onto monotonic time basis
tff = (nza.tf(1)+datenum(0,0,0,0,30,0):1/24:nza.tf(end)-datenum(0,0,0,0,30,0))';
nza.pf = interp1(nza.tf,nza.pf,tff);
nza.tf = tff;

%--- NZB
fid = [topdir 'NZB_water_height.csv'];
NZB = readtable(fid);
temp = table2array(NZB(:,7));
temp1 = cell2mat(temp);
temp2 = datenum([temp1(:,1:10),repmat(' ',length(temp1),1),temp1(:,12:19)]);
nzb.t = temp2;
nzb.p = table2array(NZB(:,8));
disp(median(nzb.p))

% empirically-identified offsets (instrument changes?)
nzb.i_offset=[66566,10748,112154];

% fit offset(s) with linear model
tinv=nzb.t-nzb.t(1);
for i=1:length(nzb.i_offset)
    io=nzb.i_offset(i);
    G=[[tinv(1:io);zeros(size(tinv(io+1:end)))],[ones(size(tinv(1:io)));zeros(size(tinv(io+1:end)))],...
        [zeros(size(tinv(1:io)));tinv(io+1:end)],[zeros(size(tinv(1:io)));ones(size(tinv(io+1:end)))]];
    m=inv(G'*G)*G'*nzb.p;
    fit_test=G*m;
    nzb.p=nzb.p-fit_test;
end

% decimation and tidal filter
[nzb.tf,nzb.pf,~,~] = downsample_uneven(nzb.t,nzb.p,1/24);
nzb.pf=Z_godin(nzb.pf)*100; % [hPa]
nzb.tf(isnan(nzb.pf))=[];
nzb.pf(isnan(nzb.pf))=[];

% interpolate onto monotonic time basis
tff = (nzb.tf(1)+datenum(0,0,0,0,30,0):1/24:nzb.tf(end)-datenum(0,0,0,0,30,0))';
nzb.pf = interp1(nzb.tf,nzb.pf,tff);
nzb.tf = tff;

%--- NZC
fid = [topdir 'NZC_water_height.csv'];
NZC = readtable(fid);
temp = table2array(NZC(:,7));
temp1 = cell2mat(temp);
temp2 = datenum([temp1(:,1:10),repmat(' ',length(temp1),1),temp1(:,12:19)]);
nzc.t = temp2;
nzc.p = table2array(NZC(:,8));
disp(median(nzc.p))

% empirically-identified offsets (instrument changes?)
nzc.i_offset=[69229,132355];

% fit offset(s) with linear model
tinv=nzc.t-nzc.t(1);
for i=1:length(nzc.i_offset)
    io=nzc.i_offset(i);
    G=[[tinv(1:io);zeros(size(tinv(io+1:end)))],[ones(size(tinv(1:io)));zeros(size(tinv(io+1:end)))],...
        [zeros(size(tinv(1:io)));tinv(io+1:end)],[zeros(size(tinv(1:io)));ones(size(tinv(io+1:end)))]];
    m=inv(G'*G)*G'*nzc.p;
    fit_test=G*m;
    nzc.p=nzc.p-fit_test;
end
% fit improves with second iteration
for i=1:length(nzc.i_offset)
    io=nzc.i_offset(i);
    G=[[tinv(1:io);zeros(size(tinv(io+1:end)))],[ones(size(tinv(1:io)));zeros(size(tinv(io+1:end)))],...
        [zeros(size(tinv(1:io)));tinv(io+1:end)],[zeros(size(tinv(1:io)));ones(size(tinv(io+1:end)))]];
    m=inv(G'*G)*G'*nzc.p;
    fit_test=G*m;
    nzc.p=nzc.p-fit_test;
end

% decimation and tidal filter
[nzc.tf,nzc.pf,~,~] = downsample_uneven(nzc.t,nzc.p,1/24);
nzc.pf=Z_godin(nzc.pf)*100; % [hPa]
nzc.tf(isnan(nzc.pf))=[];
nzc.pf(isnan(nzc.pf))=[];

% interpolate onto monotonic time basis
tff = (nzc.tf(1)+datenum(0,0,0,0,30,0):1/24:nzc.tf(end)-datenum(0,0,0,0,30,0))';
nzc.pf = interp1(nzc.tf,nzc.pf,tff);
nzc.tf = tff;

%--- NZD
fid = [topdir 'NZD_water_height.csv'];
NZD = readtable(fid);
temp = table2array(NZD(:,7));
temp1 = cell2mat(temp);
temp2 = datenum([temp1(:,1:10),repmat(' ',length(temp1),1),temp1(:,12:19)]);
nzd.t = temp2;
nzd.p = table2array(NZD(:,8));
disp(median(nzd.p))

% empirically-identified offsets (instrument changes?)
nzd.t(66076:66079)=[];
nzd.p(66076:66079)=[];
nzd.i_offset=66075;

% fit offset(s) with linear model
tinv=nzd.t-nzd.t(1);
for i=1:length(nzd.i_offset)
    io=nzd.i_offset(i);
    G=[[tinv(1:io);zeros(size(tinv(io+1:end)))],[ones(size(tinv(1:io)));zeros(size(tinv(io+1:end)))],...
        [zeros(size(tinv(1:io)));tinv(io+1:end)],[zeros(size(tinv(1:io)));ones(size(tinv(io+1:end)))]];
    m=inv(G'*G)*G'*nzd.p;
    fit_test=G*m;
    nzd.p=nzd.p-fit_test;
end

% decimation and tidal filter
[nzd.tf,nzd.pf,~,~] = downsample_uneven(nzd.t,nzd.p,1/24);
nzd.pf=Z_godin(nzd.pf)*100; % [hPa]
nzd.tf(isnan(nzd.pf))=[];
nzd.pf(isnan(nzd.pf))=[];

% interpolate onto monotonic time basis
tff = (nzd.tf(1)+datenum(0,0,0,0,30,0):1/24:nzd.tf(end)-datenum(0,0,0,0,30,0))';
nzd.pf = interp1(nzd.tf,nzd.pf,tff);
nzd.tf = tff;

%--- NZE
fid = [topdir 'NZE_water_height.csv'];
NZE = readtable(fid);
temp = table2array(NZE(:,7));
temp1 = cell2mat(temp);
temp2 = datenum([temp1(:,1:10),repmat(' ',length(temp1),1),temp1(:,12:19)]);
nze.t = temp2;
nze.p = table2array(NZE(:,8));
disp(median(nze.p))

% empirically-identified offsets (instrument changes?)
nze.t(137011:137014)=[];
nze.p(137011:137014)=[];
nze.i_offset=[67556,137010];

% fit offset(s) with linear model
tinv=nze.t-nze.t(1);
for i=1:length(nze.i_offset)
    io=nze.i_offset(i);
    G=[[tinv(1:io);zeros(size(tinv(io+1:end)))],[ones(size(tinv(1:io)));zeros(size(tinv(io+1:end)))],...
        [zeros(size(tinv(1:io)));tinv(io+1:end)],[zeros(size(tinv(1:io)));ones(size(tinv(io+1:end)))]];
    m=inv(G'*G)*G'*nze.p;
    fit_test=G*m;
    nze.p=nze.p-fit_test;
end

% decimation and tidal filter
[nze.tf,nze.pf,~,~] = downsample_uneven(nze.t,nze.p,1/24);
nze.pf=Z_godin(nze.pf)*100; % [hPa]
nze.tf(isnan(nze.pf))=[];
nze.pf(isnan(nze.pf))=[];

% interpolate onto monotonic time basis
tff = (nze.tf(1)+datenum(0,0,0,0,30,0):1/24:nze.tf(end)-datenum(0,0,0,0,30,0))';
nze.pf = interp1(nze.tf,nze.pf,tff);
nze.tf = tff;

%--- NZF
fid = [topdir 'NZF_water_height.csv'];
NZF = readtable(fid);
temp = table2array(NZF(:,7));
temp1 = cell2mat(temp);
temp2 = datenum([temp1(:,1:10),repmat(' ',length(temp1),1),temp1(:,12:19)]);
nzf.t = temp2;
nzf.p = table2array(NZF(:,8));
disp(median(nzf.p))

% empirically-identified offsets (instrument changes?)
nzf.t(1:11625)=[];
nzf.p(1:11625)=[];
nzf.i_offset=36669;

% fit offset(s) with linear model
tinv=nzf.t-nzf.t(1);
for i=1:length(nzf.i_offset)
    io=nzf.i_offset(i);
    G=[[tinv(1:io);zeros(size(tinv(io+1:end)))],[ones(size(tinv(1:io)));zeros(size(tinv(io+1:end)))],...
        [zeros(size(tinv(1:io)));tinv(io+1:end)],[zeros(size(tinv(1:io)));ones(size(tinv(io+1:end)))]];
    m=inv(G'*G)*G'*nzf.p;
    fit_test=G*m;
    nzf.p=nzf.p-fit_test;
end

% decimation and tidal filter
[nzf.tf,nzf.pf,~,~] = downsample_uneven(nzf.t,nzf.p,1/24);
nzf.pf=Z_godin(nzf.pf)*100; % [hPa]
nzf.tf(isnan(nzf.pf))=[];
nzf.pf(isnan(nzf.pf))=[];

% interpolate onto monotonic time basis
tff = (nzf.tf(1)+datenum(0,0,0,0,30,0):1/24:nzf.tf(end)-datenum(0,0,0,0,30,0))';
nzf.pf = interp1(nzf.tf,nzf.pf,tff);
nzf.tf = tff;

%--- NZG
fid = [topdir 'NZG_water_height.csv'];
NZG = readtable(fid);
temp = table2array(NZG(:,7));
temp1 = cell2mat(temp);
temp2 = datenum([temp1(:,1:10),repmat(' ',length(temp1),1),temp1(:,12:19)]);
nzg.t = temp2;
nzg.p = table2array(NZG(:,8));
disp(median(nzg.p))

% empirically-identified offsets (instrument changes?)
nzg.i_offset=53094;

% fit offset(s) with linear model
tinv=nzg.t-nzg.t(1);
for i=1:length(nzg.i_offset)
    io=nzg.i_offset(i);
    G=[[tinv(1:io);zeros(size(tinv(io+1:end)))],[ones(size(tinv(1:io)));zeros(size(tinv(io+1:end)))],...
        [zeros(size(tinv(1:io)));tinv(io+1:end)],[zeros(size(tinv(1:io)));ones(size(tinv(io+1:end)))]];
    m=inv(G'*G)*G'*nzg.p;
    fit_test=G*m;
    nzg.p=nzg.p-fit_test;
end

% decimation and tidal filter
[nzg.tf,nzg.pf,~,~] = downsample_uneven(nzg.t,nzg.p,1/24);
nzg.pf=Z_godin(nzg.pf)*100; % [hPa]
nzg.tf(isnan(nzg.pf))=[];
nzg.pf(isnan(nzg.pf))=[];

% interpolate onto monotonic time basis
tff = (nzg.tf(1)+datenum(0,0,0,0,30,0):1/24:nzg.tf(end)-datenum(0,0,0,0,30,0))';
nzg.pf = interp1(nzg.tf,nzg.pf,tff);
nzg.tf = tff;

%--- NZH
fid = [topdir 'NZH_water_height.csv'];
NZH = readtable(fid);
temp = table2array(NZH(:,7));
temp1 = cell2mat(temp);
temp2 = datenum([temp1(:,1:10),repmat(' ',length(temp1),1),temp1(:,12:19)]);
nzh.t = temp2;
nzh.p = table2array(NZH(:,8));
disp(median(nzh.p))

% empirically-identified offsets (instrument changes?)
nzh.t(1:9847)=[];
nzh.p(1:9847)=[];
nzh.i_offset=[]; % none!

% fit with linear model (no offsets)
tinv=nzh.t-nzh.t(1);
G=[tinv,ones(size(tinv))];
m=inv(G'*G)*G'*nzh.p;
fit_test=G*m;
nzh.p=nzh.p-fit_test;

% decimation and tidal filter
[nzh.tf,nzh.pf,~,~] = downsample_uneven(nzh.t,nzh.p,1/24);
nzh.pf=Z_godin(nzh.pf)*100; % [hPa]
nzh.tf(isnan(nzh.pf))=[];
nzh.pf(isnan(nzh.pf))=[];

% interpolate onto monotonic time basis
tff = (nzh.tf(1)+datenum(0,0,0,0,30,0):1/24:nzh.tf(end)-datenum(0,0,0,0,30,0))';
nzh.pf = interp1(nzh.tf,nzh.pf,tff);
nzh.tf = tff;

%--- NZI
fid = [topdir 'NZI_water_height.csv'];
NZI = readtable(fid);
temp = table2array(NZI(:,7));
temp1 = cell2mat(temp);
temp2 = datenum([temp1(:,1:10),repmat(' ',length(temp1),1),temp1(:,12:19)]);
nzi.t = temp2;
nzi.p = table2array(NZI(:,8));
disp(median(nzi.p))

% empirically-identified offsets (instrument changes?)
nzi.t(1:23145)=[];
nzi.p(1:23145)=[];
nzi.i_offset=[]; % none!

% fit with linear model (no offsets)
tinv=nzi.t-nzi.t(1);
G=[tinv,ones(size(tinv))];
m=inv(G'*G)*G'*nzi.p;
fit_test=G*m;
nzi.p=nzi.p-fit_test;

% decimation and tidal filter
[nzi.tf,nzi.pf,~,~] = downsample_uneven(nzi.t,nzi.p,1/24);
nzi.pf=Z_godin(nzi.pf)*100; % [hPa]
nzi.tf(isnan(nzi.pf))=[];
nzi.pf(isnan(nzi.pf))=[];

% interpolate onto monotonic time basis
tff = (nzi.tf(1)+datenum(0,0,0,0,30,0):1/24:nzi.tf(end)-datenum(0,0,0,0,30,0))';
nzi.pf = interp1(nzi.tf,nzi.pf,tff);
nzi.tf = tff;

%-----PLOTTING-----%

figure(2); clf; hold on
plot(nza.tf,nza.pf,'linewidth',1)
plot(nzb.tf,nzb.pf+10,'linewidth',1)
plot(nzc.tf,nzc.pf+20,'linewidth',1)
plot(nzd.tf,nzd.pf+30,'linewidth',1)
plot(nze.tf,nze.pf+40,'linewidth',1)
plot(nzf.tf,nzf.pf+50,'linewidth',1)
plot(nzg.tf,nzg.pf+60,'linewidth',1)
plot(nzh.tf,nzh.pf+70,'linewidth',1)
plot(nzi.tf,nzi.pf+80,'linewidth',1)
legend('NZA','NZB','NZC','NZD','NZE','NZF','NZG','NZH','NZI','location','northwest')
datetick('x',3)
ylabel('\DeltaP (cm)')
set(gca,'fontsize',14)
box on; grid on

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print(figdir,'-dpng','-r300')

% combine and save data
t={nza.t,nzb.t,nzc.t,nzd.t,nze.t,nzf.t,nzg.t,nzh.t,nzi.t};
p={nza.p,nzb.p,nzc.p,nzd.p,nze.p,nzf.p,nzg.p,nzh.p,nzi.p};
tf={nza.tf,nzb.tf,nzc.tf,nzd.tf,nze.tf,nzf.tf,nzg.tf,nzh.tf,nzi.tf};
pf={nza.pf,nzb.pf,nzc.pf,nzd.pf,nze.pf,nzf.pf,nzg.pf,nzh.pf,nzi.pf};
save(svdir,'t','p','tf','pf','staname','stalat','stalon','stadepth')

% write to text files
for i=1:length(t)
    % decimated
    writematrix([t{i},p{i}],[wrtdir staname{i} '_1hr'])
    % decimated, detided
    writematrix([tf{i},pf{i}],[wrtdir staname{i} '_1hr_detided'])
end

%% optionally, some simple comparisons with other datasets

if true
    return
end

% load up abyssal POBS data from 2022-2023
load('../../A-0-A/stitched_data/POBS2.mat')
A2=dataf;
load('../../A-0-A/stitched_data/POBS3.mat')
A3=dataf;
load('../../A-0-A/stitched_data/POBS4.mat')
A4=dataf;

%----- differences with NZB
[t{1},ia{1},ib{1}]=intersect(round(nzb.tf,6),round(A2.tf,6));
p1{1}=detrend(nzb.pf(ia{1}));
p2{1}=detrend(A2.p2f(ib{1}));
[t{2},ia{2},ib{2}]=intersect(round(nzb.tf,6),round(A3.tf,6));
p1{2}=detrend(nzb.pf(ia{2}));
p2{2}=detrend(A3.p2f(ib{2}));
[t{3},ia{3},ib{3}]=intersect(round(nzb.tf,6),round(A4.tf,6));
p1{3}=detrend(nzb.pf(ia{3}));
p2{3}=detrend(A4.p2f(ib{3}));

figure(9); clf; hold on
plot(t{1},p1{1},'b')
plot(t{1},p2{1},'r')
plot(t{1},p1{1}-p2{1},'k','linewidth',1)
text(t{1}(end)+10,p1{1}(end),'POBS2')
plot(t{2},p1{2}+10,'b')
plot(t{2},p2{2}+10,'r')
plot(t{2},p1{2}-p2{2}+10,'k','linewidth',1)
text(t{2}(end)+10,p1{2}(end)+10,'POBS3')
plot(t{3},p1{3}+20,'b')
plot(t{3},p2{3}+20,'r')
plot(t{3},p1{3}-p2{3}+20,'k','linewidth',1)
text(t{3}(end)+10,p1{3}(end)+20,'POBS4')
ylabel('P (hPa)')
legend('DART','POBS','difference','location','northeast')
datetick('x',6)
title('NZB vs. abyssal POBS')
set(gca,'fontsize',14)
box on; grid on

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print('../figures/exploratory/DART/NZB_comps','-dpng','-r300')

%----- differences with NZC
[t{1},ia{1},ib{1}]=intersect(round(nzc.tf,6),round(A2.tf,6));
p1{1}=detrend(nzc.pf(ia{1}));
p2{1}=detrend(A2.p2f(ib{1}));
[t{2},ia{2},ib{2}]=intersect(round(nzc.tf,6),round(A3.tf,6));
p1{2}=detrend(nzc.pf(ia{2}));
p2{2}=detrend(A3.p2f(ib{2}));
[t{3},ia{3},ib{3}]=intersect(round(nzc.tf,6),round(A4.tf,6));
p1{3}=detrend(nzc.pf(ia{3}));
p2{3}=detrend(A4.p2f(ib{3}));

figure(8); clf; hold on
plot(t{1},p1{1},'b')
plot(t{1},p2{1},'r')
plot(t{1},p1{1}-p2{1},'k','linewidth',1)
text(t{1}(end)+10,p1{1}(end),'POBS2')
plot(t{2},p1{2}+10,'b')
plot(t{2},p2{2}+10,'r')
plot(t{2},p1{2}-p2{2}+10,'k','linewidth',1)
text(t{2}(end)+10,p1{2}(end)+10,'POBS3')
plot(t{3},p1{3}+20,'b')
plot(t{3},p2{3}+20,'r')
plot(t{3},p1{3}-p2{3}+20,'k','linewidth',1)
text(t{3}(end)+10,p1{3}(end)+20,'POBS4')
ylabel('P (hPa)')
legend('DART','POBS','difference','location','northeast')
datetick('x',6)
title('NZC vs. abyssal POBS')
set(gca,'fontsize',14)
box on; grid on

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print('../figures/exploratory/DART/NZC_comps','-dpng','-r300')
