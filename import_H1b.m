% import_H1b.m
%
% Read HOBITSS I (2014-2015) data matfiles
%

clear; close all

%% HOBITSS I (2014-2015)

topdir='/Volumes/Gorgoroth/apg_database/original/2014-2015_HOBITSS-I/';
wrtdir='/Volumes/Gorgoroth/apg_database/processed/2014-2015_HOBITSS-Ib/';
figdir='../figures/exploratory/HOBITSS_I/hobitss1b_stack';
svdir='../processed_data/HOBITSS_Ib';

% file info
dirinfo=dir([topdir 'mat']);
dirinfo(1:2)=[]; % these are just the hidden files
staname={dirinfo.name};
for j=1:length(staname)
    staname{j}=staname{j}(9:end-4);
end

% station info (need to get this from cruise logs, etc.)
stalat=[];
stalon=[];
stadepth=[];

t={}; p={}; tf={}; pf={}; T={}; Tf={};
for i=1:length(dirinfo)
    load([dirinfo(i).folder '/' dirinfo(i).name])
    
    t{i}=times;
    p{i}=pressures;
    T{i}=temperatures;

    % tidal filter
    pf{i}=Z_godin(pressures);
    Tf{i}=Z_godin(temperatures);

    % remove NaNs from tidal filter
    t{i}(isnan(pf{i}))=[];
    T{i}(isnan(pf{i}))=[];
    Tf{i}(isnan(pf{i}))=[];
    p{i}(isnan(pf{i}))=[];
    pf{i}(isnan(pf{i}))=[];

    % interpolate onto monotonic time basis
    tf{i} = t{i}(1)+datenum(0,0,0,0,30,0):1/24:t{i}(end)-datenum(0,0,0,0,30,0);
    pf{i} = interp1(t{i},pf{i},tf{i});
    Tf{i} = interp1(t{i},Tf{i},tf{i});
end

save(svdir,'t','p','tf','pf','T','Tf','staname','stalat','stalon','stadepth')

% write to text files
for i=1:length(t)
    % decimated
    writematrix([t{i},p{i},T{i}],[wrtdir staname{i} '_1hr'])
    % decimated, detided
    writematrix([tf{i}',pf{i}',Tf{i}'],[wrtdir staname{i} '_1hr_detided'])
end


%-----PLOTTING-----%
figure(2); clf;
subplot(121); hold on
for k=1:10
    p_plot = detrend(pf{k});
    plot(tf{k},p_plot+20*(k-1),'linewidth',1)
    text(tf{k}(end)+10,p_plot(end-100)+20*(k-1),staname{k},'fontsize',14)
end
xlim([datenum(2014,05,01) datenum(2015,09,01)])
ylim([-20 200])
datetick('x',6,'keeplimits')
set(gca,'fontsize',14)
ylabel('\DeltaP (hPa)')
title('HOBITSS 2014-2015')
box on; grid on

subplot(122); hold on
for k=11:20
    p_plot = detrend(pf{k});
    plot(tf{k},p_plot+20*(k-11),'linewidth',1)
    text(tf{k}(end)+10,p_plot(end-100)+20*(k-11),staname{k},'fontsize',14)
end
xlim([datenum(2014,05,01) datenum(2015,09,01)])
ylim([-20 200])
datetick('x',6,'keeplimits')
set(gca,'fontsize',14)
ylabel('\DeltaP (hPa)')
title('HOBITSS 2014-2015')
box on; grid on

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print(figdir,'-dpng','-r300')