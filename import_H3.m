% import_H3.m
%
% Read HOBITSS III (2016-2017) data from text/csv files
%

clear; close all

%% HOBITSS III (2016-2017)

topdir='/Volumes/Gorgoroth/apg_database/original/2016-2017_HOBITSS-III/';
wrtdir='/Volumes/Gorgoroth/apg_database/processed/2016-2017_HOBITSS-III/';
figdir='../figures/exploratory/HOBITSS_III/hobitss3_stack';
svdir='../processed_data/HOBITSS_III';

% station info
staname={'KU16-2','KU16-3','KU16-4','KU16-5'};
stalat=[-38.873,-38.888,-38.708,-38.722];
stalon=[178.873,178.757,178.666,178.900];
stadepth=[2146.6,1363.0,1051.6,2468.6];

%---KU16-2---%
% read in data
fid = fopen([topdir 'data_out_KU16-2.txt'],'r');
fspec = '%f %f %f'; % [P T t]
sizeA = [3 Inf];
A = fscanf(fid,fspec,sizeA);
fclose(fid);

% time is as yymmddHHMMSS
temp.base = A(3,:);
temp.yr = floor(temp.base/10^10);
temp.base = temp.base-temp.yr*10^10;
temp.mnth = floor(temp.base/10^8);
temp.base = temp.base-temp.mnth*10^8;
temp.dy = floor(temp.base/10^6);
temp.base = temp.base-temp.dy*10^6;
temp.hr = floor(temp.base/10^4);
temp.base = temp.base-temp.hr*10^4;
temp.min = floor(temp.base/10^2);
temp.sec = temp.base-temp.min*10^2;
A(3,:) = datenum(2000+temp.yr,temp.mnth,temp.dy,temp.hr,temp.min,temp.sec);

% empirical trimming
A = A(:,153355:32007424);

% remove duplicate/nonmonotonic timesteps
ilist=find(diff(A(1,:))<=0);
if ~isempty(ilist)
    A(:,ilist+1)=[];
end

% interpolate to 1 Hz
AA=A(1,1):1/86400:A(1,end);
AA(2,:)=interp1(A(1,:),A(2,:),AA); % P
AA(3,:)=interp1(A(1,:),A(3,:),AA(1,:)); % T

% write to text file, with pressure as Pa
writematrix([AA(3,:);AA(2,:)*100;AA(1,:)]',[wrtdir 'KU16-2_1Hz']) % [t P T]

% decimation loop
ta=[];
a=[];
Ta=[];
i1 = 1;
d2 = floor(A(3,1))+1;
while i1<length(A)
    i2 = find(A(3,:)>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segA,~,~] = downsample_uneven(A(3,i1:i2-1),A(1:2,i1:i2-1),1/24);
    if length(segt)>24
        keyboard
    end
    ta=cat(2,ta,segt);
    a=cat(2,a,segA(1,:));
    Ta=cat(2,Ta,segA(2,:));
    i1=i2;
    d2=floor(A(3,i2))+1;
end

% tidal filter
af=Z_godin(a);
Taf=Z_godin(Ta);

% remove NaNs from tidal filter
ta(isnan(af))=[];
Ta(isnan(af))=[];
Taf(isnan(af))=[];
a(isnan(af))=[];
af(isnan(af))=[];

% interpolate onto monotonic time basis on the hour
taf = ta(1)+datenum(0,0,0,0,30,0):1/24:ta(end)-datenum(0,0,0,0,30,0);
af = interp1(ta,af,taf);
Taf = interp1(ta,Taf,taf);

% this yields ~369 days of data, which matches cruise reports
clearvars('temp','A','AA')

%---KU16-3---%
% read in data
fid = fopen([topdir 'data_out_KU16-3.txt'],'r');
fspec = '%f %f %f'; % [P T t]
sizeB = [3 Inf];
B = fscanf(fid,fspec,sizeB);
fclose(fid);

% time is as yymmddHHMMSS
temp.base = B(3,:);
temp.yr = floor(temp.base/10^10);
temp.base = temp.base-temp.yr*10^10;
temp.mnth = floor(temp.base/10^8);
temp.base = temp.base-temp.mnth*10^8;
temp.dy = floor(temp.base/10^6);
temp.base = temp.base-temp.dy*10^6;
temp.hr = floor(temp.base/10^4);
temp.base = temp.base-temp.hr*10^4;
temp.min = floor(temp.base/10^2);
temp.sec = temp.base-temp.min*10^2;
B(3,:) = datenum(2000+temp.yr,temp.mnth,temp.dy,temp.hr,temp.min,temp.sec);

% empirical trimming
B = B(:,134274:15936673);

% remove duplicate/nonmonotonic timesteps
ilist=find(diff(B(1,:))<=0);
if ~isempty(ilist)
    B(:,ilist+1)=[];
end

% interpolate to 1 Hz
BB=B(1,1):1/86400:B(1,end);
BB(2,:)=interp1(B(1,:),B(2,:),BB); % P
BB(3,:)=interp1(B(1,:),B(3,:),BB(1,:)); % T

% write to text file, with pressure as Pa
writematrix([BB(3,:);BB(2,:)*100;BB(1,:)]',[wrtdir 'KU16-3_1Hz']) % [t P T]

% decimation loop
tb=[];
b=[];
Tb=[];
i1 = 1;
d2 = floor(B(3,1))+1;
while i1<length(B)
    i2 = find(B(3,:)>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segB,~,~] = downsample_uneven(B(3,i1:i2-1),B(1:2,i1:i2-1),1/24);
    if length(segt)>24
        keyboard
    end
    tb=cat(2,tb,segt);
    b=cat(2,b,segB(1,:));
    Tb=cat(2,Tb,segB(2,:));
    i1=i2;
    d2=floor(B(3,i2))+1;
end

% tidal filter
bf=Z_godin(b);
Tbf=Z_godin(Tb);

% remove NaNs from tidal filter
tb(isnan(bf))=[];
Tb(isnan(bf))=[];
Tbf(isnan(bf))=[];
b(isnan(bf))=[];
bf(isnan(bf))=[];

% interpolate onto monotonic time basis
tbf = tb(1)+datenum(0,0,0,0,30,0):1/24:tb(end)-datenum(0,0,0,0,30,0);
bf = interp1(tb,bf,tbf);
Tbf = interp1(tb,Tbf,tbf);

% this yields ~365 days of data
clearvars('temp','B','BB')

%---KU16-4---%
% read in data
fid = fopen([topdir 'data_out_KU16-4.txt'],'r');
fspec = '%f %f %f'; % [P T t]
sizeC = [3 Inf];
C = fscanf(fid,fspec,sizeC);
fclose(fid);

% time is as yymmddHHMMSS
temp.base = C(3,:);
temp.yr = floor(temp.base/10^10);
temp.base = temp.base-temp.yr*10^10;
temp.mnth = floor(temp.base/10^8);
temp.base = temp.base-temp.mnth*10^8;
temp.dy = floor(temp.base/10^6);
temp.base = temp.base-temp.dy*10^6;
temp.hr = floor(temp.base/10^4);
temp.base = temp.base-temp.hr*10^4;
temp.min = floor(temp.base/10^2);
temp.sec = temp.base-temp.min*10^2;
C(3,:) = datenum(2000+temp.yr,temp.mnth,temp.dy,temp.hr,temp.min,temp.sec);

% empirical trimming
C = C(:,134013:16021708);

% remove duplicate/nonmonotonic timesteps
ilist=find(diff(C(1,:))<=0);
if ~isempty(ilist)
    C(:,ilist+1)=[];
end

% interpolate to 1 Hz
CC=C(1,1):1/86400:C(1,end);
CC(2,:)=interp1(C(1,:),C(2,:),CC); % P
CC(3,:)=interp1(C(1,:),C(3,:),CC(1,:)); % T

% write to text file, with pressure as Pa
writematrix([CC(3,:);CC(2,:)*100;CC(1,:)]',[wrtdir 'KU16-4_1Hz']) % [t P T]

% decimation loop
tc=[];
c=[];
Tc=[];
i1 = 1;
d2 = floor(C(3,1))+1;
while i1<length(C)
    i2 = find(C(3,:)>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segC,~,~] = downsample_uneven(C(3,i1:i2-1),C(1:2,i1:i2-1),1/24);
    if length(segt)>24
        keyboard
    end
    tc=cat(2,tc,segt);
    c=cat(2,c,segC(1,:));
    Tc=cat(2,Tc,segC(2,:));
    i1=i2;
    d2=floor(C(3,i2))+1;
end

% tidal filter
cf=Z_godin(c);
Tcf=Z_godin(Tc);

% remove NaNs from tidal filter
tc(isnan(cf))=[];
Tc(isnan(cf))=[];
Tcf(isnan(cf))=[];
c(isnan(cf))=[];
cf(isnan(cf))=[];

% interpolate onto monotonic time basis
tcf = tc(1)+datenum(0,0,0,0,30,0):1/24:tc(end)-datenum(0,0,0,0,30,0);
cf = interp1(tc,cf,tcf);
Tcf = interp1(tc,Tcf,tcf);

% this yields ~368 days of data
clearvars('temp','C','CC')

%---KU16-5---%
% read in data
fid = fopen([topdir 'data_out_KU16-5.txt'],'r');
fspec = '%f %f %f'; % [P T t]
sizeD = [3 Inf];
D = fscanf(fid,fspec,sizeD);
D = D(:,136:end); % data file skips ahead in time
fclose(fid);

% time is as yymmddHHMMSS
temp.base = D(3,:);
temp.yr = floor(temp.base/10^10);
temp.base = temp.base-temp.yr*10^10;
temp.mnth = floor(temp.base/10^8);
temp.base = temp.base-temp.mnth*10^8;
temp.dy = floor(temp.base/10^6);
temp.base = temp.base-temp.dy*10^6;
temp.hr = floor(temp.base/10^4);
temp.base = temp.base-temp.hr*10^4;
temp.min = floor(temp.base/10^2);
temp.sec = temp.base-temp.min*10^2;
D(3,:) = datenum(2000+temp.yr,temp.mnth,temp.dy,temp.hr,temp.min,temp.sec);

% empirical trimming
D = D(:,115884:16020088);

% remove duplicate/nonmonotonic timesteps
ilist=find(diff(D(1,:))<=0);
if ~isempty(ilist)
    D(:,ilist+1)=[];
end

% interpolate to 1 Hz
DD=D(1,1):1/86400:D(1,end);
DD(2,:)=interp1(D(1,:),D(2,:),DD); % P
DD(3,:)=interp1(D(1,:),D(3,:),DD(1,:)); % T

% write to text file, with pressure as Pa
writematrix([DD(3,:);DD(2,:)*100;DD(1,:)]',[wrtdir 'KU16-5_1Hz']) % [t P T]

% decimation loop
td=[];
d=[];
Td=[];
i1 = 1;
d2 = floor(D(3,1))+1;
while i1<length(D)
    i2 = find(D(3,:)>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segD,~,~] = downsample_uneven(D(3,i1:i2-1),D(1:2,i1:i2-1),1/24);
    if length(segt)>24
        keyboard
    end
    if length(segt)~=length(segD)
        keyboard
    end
    td=cat(2,td,segt);
    d=cat(2,d,segD(1,:));
    Td=cat(2,Td,segD(2,:));
    i1=i2;
    d2=floor(D(3,i2))+1;
end

% tidal filter
df=Z_godin(d);
Tdf=Z_godin(Td);

% remove NaNs from tidal filter
td(isnan(df))=[];
Td(isnan(df))=[];
Tdf(isnan(df))=[];
d(isnan(df))=[];
df(isnan(df))=[];

% interpolate onto monotonic time basis
tdf = td(1)+datenum(0,0,0,0,30,0):1/24:td(end)-datenum(0,0,0,0,30,0);
df = interp1(td,df,tdf);
Tdf = interp1(td,Tdf,tdf);

% this yields ~368 days of data
clearvars('temp','D','DD')

%-----PLOTTING-----%
figure(1); clf; hold on
af_plot = detrend(af);
plot(taf,af_plot,'linewidth',1)
text(taf(end)+10,af_plot(end-100),{staname{1};[num2str(stadepth(1)) ' m']},'fontsize',14)
df_plot = detrend(df);
plot(tdf,df_plot-20,'linewidth',1)
text(tdf(end)+10,df_plot(end-100)-20,{staname{4};[num2str(stadepth(4)) ' m']},'fontsize',14)
bf_plot = detrend(bf);
plot(tbf,bf_plot+20,'linewidth',1)
text(tbf(end)+10,bf_plot(end-100)+20,{staname{2};[num2str(stadepth(2)) ' m']},'fontsize',14)
cf_plot = detrend(cf);
plot(tcf,cf_plot+40,'linewidth',1)
text(tcf(end)+10,cf_plot(end-100)+40,{staname{3};[num2str(stadepth(3)) ' m']},'fontsize',14)
xlim([datenum(2016,06,01) datenum(2017,09,01)])
datetick('x','keeplimits')
set(gca,'fontsize',14)
ylabel('\DeltaP (hPa)')
title('HOBITSS 2016-2017')
box on; grid on

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print(pltdir,'-dpng','-r300')

% combine and save data
t={ta,tb,tc,td};
p={a,b,c,d};
T={Ta,Tb,Tc,Td};
tf={taf,tbf,tcf,tdf};
pf={af,bf,cf,df};
Tf={Taf,Tbf,Tcf,Tdf};
save(svdir,'t','p','T','tf','pf','Tf','staname','stalat','stalon','stadepth')

% write to text files
for i=1:length(t)
    % decimated
    writematrix([t{i}',p{i}',T{i}'],[wrtdir staname{i} '_1hr'])
    % decimated, detided
    writematrix([tf{i}',pf{i}',Tf{i}'],[wrtdir staname{i} '_1hr_detided'])
end