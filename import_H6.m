% import_H6.m
%
% Read HOBITSS VI (2019-2020) data from text/csv files
%

clear; close all

%% HOBITSS VI (2019-2020)

topdir='/Volumes/Gorgoroth/apg_database/original/2019-2020_HOBITSS-VI/';
wrtdir='/Volumes/Gorgoroth/apg_database/processed/2019-2020_HOBITSS-VI/';
figdir='../figures/exploratory/HOBITSS_VI/hobitss6_stack';
svdir='../processed_data/HOBITSS_VI';

% station info
staname={'KU19-1','KU19-2','KU19-3','KU19-4','KU19-5'};
stalat=[-38.75366,-38.87643,-38.89205,-38.79347,-38.72192];
stalon=[179.01582,178.84522,178.75519,178.67079,178.89374];
stadepth=[3537.3,1913.6,1358,1077.8,2453.6];

%---KU19-1---%
% read in data
fid = fopen([topdir 'KU19-1_126646_RAW.dat'],'r');
fspec = '%s %s %f %f'; % [date time P T]
A = textscan(fid,fspec);
fclose(fid);
t_str=cat(2,cat(1,A{1}{1:end-1}),repmat(' ',length(A{1})-1,1),cat(1,A{2}{1:end-1}));
tA = datenum(t_str,'yyyy/mm/dd HH:MM:SS');
A = [tA';A{3}';A{4}']; % [t P T]

% empirical trimming
A = A(:,381798:32864821);

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
writematrix([AA(1,:);AA(2,:)*100;AA(3,:)]',[wrtdir 'KU19-1_1Hz']) % [t P T]

% decimation loop
ta=[];
a=[];
Ta=[];
i1 = 1;
while i1<length(A)
    [~,i2] = min(abs(A(1,:)-(floor(A(1,i1))+1)));
    [segt,segA,~,~] = downsample_uneven(A(1,i1:i2-1),A(2:3,i1:i2-1),1/24);
    ta=cat(2,ta,segt);
    a=cat(2,a,segA(1,:));
    Ta=cat(2,Ta,segA(2,:));
    i1=i2;
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

% interpolate onto monotonic time basis
taf = ta(1)+datenum(0,0,0,0,30,0):1/24:ta(end)-datenum(0,0,0,0,30,0);
af = interp1(ta,af,taf);
Taf = interp1(ta,Taf,taf);

clearvars('t_str','tA','A','AA')

%---KU19-2---%
% read in data
fid = fopen([topdir 'KU19-2_131477_RAW.dat'],'r');
fspec = '%f %f %f'; % [P T t]
B = textscan(fid,fspec);
fclose(fid);
B = [B{3},B{1},B{2}]'; % [t P T]

% time is as yymmddHHMMSS
temp.base = B(1,:);
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
B(1,:) = datenum(2000+temp.yr,temp.mnth,temp.dy,temp.hr,temp.min,temp.sec);

% empirical trimming
B = B(:,17035:16273967);

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
writematrix([BB(1,:);BB(2,:)*100;BB(3,:)]',[wrtdir 'KU19-2_1Hz']) % [t P T]

% decimation loop
tb=[];
b=[];
Tb=[];
i1 = 1;
d2 = floor(B(1,1))+1;
while i1<length(B)
    i2 = find(B(1,:)>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segB,~,~] = downsample_uneven(B(1,i1:i2-1),B(2:3,i1:i2-1),1/24);
    if length(segt)>24
        keyboard
    end
    tb=cat(2,tb,segt);
    b=cat(2,b,segB(1,:));
    Tb=cat(2,Tb,segB(2,:));
    i1=i2;
    d2=floor(B(1,i2))+1;
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

clearvars('temp','B','BB')

%---KU19-3---%
% read in data
fid = fopen([topdir 'KU19-3_131479_RAW.dat'],'r');
fspec = '%f %f %f'; % [P T t]
C = textscan(fid,fspec);
fclose(fid);
C = [C{3},C{1},C{2}]'; % [t P T]

% time is as yymmddHHMMSS
temp.base = C(1,:);
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
C(1,:) = datenum(2000+temp.yr,temp.mnth,temp.dy,temp.hr,temp.min,temp.sec);

% empirical trimming
C = C(:,21653:16408572);

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
writematrix([CC(1,:);CC(2,:)*100;CC(3,:)]',[wrtdir 'KU19-3_1Hz']) % [t P T]

% decimation loop
tc=[];
c=[];
Tc=[];
i1 = 1;
d2 = floor(C(1,1))+1;
while i1<length(C)
    i2 = find(C(1,:)>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segC,~,~] = downsample_uneven(C(1,i1:i2-1),C(2:3,i1:i2-1),1/24);
    if length(segt)>24
        keyboard
    end
    tc=cat(2,tc,segt);
    c=cat(2,c,segC(1,:));
    Tc=cat(2,Tc,segC(2,:));
    i1=i2;
    d2=floor(C(1,i2))+1;
end

% tidal filter
cf=Z_godin(c);
Tcf=Z_godin(Tc);

% remove NaNs from tidal filter
tc(isnan(cf))=[];
c(isnan(cf))=[];
Tc(isnan(cf))=[];
Tcf(isnan(cf))=[];
cf(isnan(cf))=[];

% interpolate onto monotonic time basis
tcf = tc(1)+datenum(0,0,0,0,30,0):1/24:tc(end)-datenum(0,0,0,0,30,0);
cf = interp1(tc,cf,tcf);
Tcf = interp1(tc,Tcf,tcf);

clearvars('temp','C','CC')

%---KU19-4---%
% read in data
fid = fopen([topdir 'KU19-4_125142_RAW.dat'],'r');
fspec = '%s %s %f %f'; % [date time P T]
D = textscan(fid,fspec);
fclose(fid);
t_str=cat(2,cat(1,D{1}{1:end-1}),repmat(' ',length(D{1})-1,1),cat(1,D{2}{1:end-1}));
tD = datenum(t_str,'yyyy/mm/dd HH:MM:SS');
D = [tD';D{3}';D{4}'];

% empirical trimming
D = D(:,53440:32883871);

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
writematrix([DD(1,:);DD(2,:)*100;DD(3,:)]',[wrtdir 'KU19-4_1Hz']) % [t P T]

% decimation loop
td=[];
d=[];
Td=[];
i1 = 1;
while i1<length(D)
    [~,i2] = min(abs(D(1,:)-(floor(D(1,i1))+1)));
    [segt,segD,~,~] = downsample_uneven(D(1,i1:i2-1),D(2:3,i1:i2-1),1/24);
    td=cat(2,td,segt);
    d=cat(2,d,segD(1,:));
    Td=cat(2,Td,segD(2,:));
    i1=i2;
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

clearvars('t_str','tD','D','DD')

%---KU19-5---%
% read in data
fid = fopen([topdir 'KU19-5_126652_RAW.dat'],'r');
fspec = '%f %f %f'; % [P T t]
E = textscan(fid,fspec);
fclose(fid);
E = [E{3},E{1},E{2}]'; % [t P T]

% time is as yymmddHHMMSS
temp.base = E(1,:);
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
E(1,:) = datenum(2000+temp.yr,temp.mnth,temp.dy,temp.hr,temp.min,temp.sec);

% empirical trimming
E = E(:,62962:16403366);

% remove duplicate/nonmonotonic timesteps
ilist=find(diff(E(1,:))<=0);
if ~isempty(ilist)
    E(:,ilist+1)=[];
end

% interpolate to 1 Hz
EE=E(1,1):1/86400:E(1,end);
EE(2,:)=interp1(E(1,:),E(2,:),EE); % P
EE(3,:)=interp1(E(1,:),E(3,:),EE(1,:)); % T

% write to text file, with pressure as Pa
writematrix([EE(1,:);EE(2,:)*100;EE(3,:)]',[wrtdir 'KU19-5_1Hz']) % [t P T]

% decimation loop
te=[];
e=[];
Te=[];
i1 = 1;
d2 = floor(E(1,1))+1;
while i1<length(E)
    i2 = find(E(1,:)>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segE,~,~] = downsample_uneven(E(1,i1:i2-1),E(2:3,i1:i2-1),1/24);
    if length(segt)>24
        keyboard
    end
    te=cat(2,te,segt);
    e=cat(2,e,segE(1,:));
    Te=cat(2,Te,segE(2,:));
    i1=i2;
    d2=floor(E(1,i2))+1;
end

% tidal filter
ef=Z_godin(e);
Tef=Z_godin(Te);

% remove NaNs from tidal filter
te(isnan(ef))=[];
Te(isnan(ef))=[];
Tef(isnan(ef))=[];
e(isnan(ef))=[];
ef(isnan(ef))=[];

% interpolate onto monotonic time basis
tef = te(1)+datenum(0,0,0,0,30,0):1/24:te(end)-datenum(0,0,0,0,30,0);
ef = interp1(te,ef,tef);
Tef = interp1(te,Tef,tef);

clearvars('temp','E','EE')

%-----PLOTTING-----%
figure(2); clf; hold on
af_plot = detrend(af);
plot(taf,af_plot,'linewidth',1)
text(taf(end)+10,af_plot(end-50),{staname{1};[num2str(stadepth(1)) ' m']},'fontsize',14)
ef_plot = detrend(ef);
plot(tef,ef_plot+20,'linewidth',1)
text(tef(end)+10,ef_plot(end-50)+20,{staname{5};[num2str(stadepth(5)) ' m']},'fontsize',14)
bf_plot = detrend(bf);
plot(tbf,bf_plot+40,'linewidth',1)
text(tbf(end)+10,bf_plot(end-50)+40,{staname{2};[num2str(stadepth(2)) ' m']},'fontsize',14)
cf_plot = detrend(cf);
plot(tcf,cf_plot+60,'linewidth',1)
text(tcf(end)+10,cf_plot(end-50)+60,{staname{3};[num2str(stadepth(3)) ' m']},'fontsize',14)
df_plot = detrend(df);
plot(tdf,df_plot+80,'linewidth',1)
text(tdf(end)+10,df_plot(end-50)+80,{staname{4};[num2str(stadepth(4)) ' m']},'fontsize',14)
xlim([datenum(2019,11,01) datenum(2021,02,01)])
datetick('x','keeplimits')
set(gca,'fontsize',14)
ylabel('\DeltaP (hPa)')
title('HOBITSS 2019-2020')
box on; grid on

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print(figdir,'-dpng','-r300')

% combine and save data
t={ta,tb,tc,td,te};
p={a,b,c,d,e};
T={Ta,Tb,Tc,Td,Te};
tf={taf,tbf,tcf,tdf,tef};
pf={af,bf,cf,df,ef};
Tf={Taf,Tbf,Tcf,Tdf,Tef};
save(svdir,'t','p','T','tf','pf','Tf','staname','stalat','stalon','stadepth')

% write to text files
for i=1:length(t)
    % decimated
    writematrix([t{i}',p{i}',T{i}'],[wrtdir staname{i} '_1hr'])
    % decimated, detided
    writematrix([tf{i}',pf{i}',Tf{i}'],[wrtdir staname{i} '_1hr_detided'])
end