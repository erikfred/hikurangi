% import_H5.m
%
% Read HOBITSS V (2018-2019) data directly from text/csv files
% Note that 'import_H5b.m' simply uses Katie Woods' pre-processed data
%

clear; close all

%% HOBITSS V (2018-2019)

topdir='/Volumes/Gorgoroth/apg_database/original/2018-2019_HOBITSS-V/';
wrtdir='/Volumes/Gorgoroth/apg_database/processed/2018-2019_HOBITSS-V/';
figdir='../figures/exploratory/HOBITSS_V/hobitss5_stack';
svdir='../processed_data/HOBITSS_V';

% station info
staname={'KU18-1','KU18-2','KU18-3','KU18-4','GNS18-01','GNS18-03','GNS18-05','LBPR18-05','GNS18-07','LBPR18-04'};
stalat=[-38.906,-38.876,-38.893,-38.794,-39.876,-39.647,-39.507,-39.693,-40.283,-39.767];
stalon=[178.982,178.845,178.755,178.671,178.666,178.149,178.530,178.692,177.523,178.485];
stadepth=[3483,1912,1357,1083,3309,779,2259,3311,1895,2207];

%---KU18-1---%
% read in data
fid = fopen([topdir 'apg_KU/KU18-1_181/181_126899_RAW.dat'],'r');
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
A(:,16928767:end)=[];
A(:,1:108350)=[];

% remove duplicate/nonmonotonic timesteps
ilist=find(diff(A(3,:))<=0);
if ~isempty(ilist)
    A(:,ilist+1)=[];
end

% interpolate to 1 Hz
AA=A(3,1):1/86400:A(3,end);
AA(2,:)=interp1(A(3,:),A(2,:),AA); % P
AA(3,:)=interp1(A(3,:),A(1,:),AA(1,:)); % T

% write to text file, with pressure as Pa
writematrix([AA(1,:);AA(2,:)*100;AA(3,:)]',[wrtdir 'KU18-1_1Hz']) % [t P T]

% decimation loop
ta=[];
Ta=[];
a=[];
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
Taf= interp1(ta,Taf,taf);

clearvars('temp','A','AA')

%---KU18-2---%
% read in data
fid = fopen([topdir 'apg_KU/KU18-2_182/182_126896_RAW.dat'],'r');
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
B(:,16834042:end)=[];
B(:,1:113411)=[];

% remove duplicate/nonmonotonic timesteps
ilist=find(diff(B(3,:))<=0);
if ~isempty(ilist)
    B(:,ilist+1)=[];
end

% interpolate to 1 Hz
BB=B(3,1):1/86400:B(3,end);
BB(2,:)=interp1(B(3,:),B(2,:),BB); % P
BB(3,:)=interp1(B(3,:),B(1,:),BB(1,:)); % T

% write to text file, with pressure as Pa
writematrix([BB(1,:);BB(2,:)*100;BB(3,:)]',[wrtdir 'KU18-2_1Hz']) % [t P T]

% decimation loop
tb=[];
Tb=[];
b=[];
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

clearvars('temp','B','BB')

%---KU18-3---%
% read in data
fid = fopen([topdir 'apg_KU/KU18-3_183/183_126897_RAW.dat'],'r');
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
C(:,16836948:end)=[];
C(:,1:105779)=[];

% remove duplicate/nonmonotonic timesteps
ilist=find(diff(C(3,:))<=0);
if ~isempty(ilist)
    C(:,ilist+1)=[];
end

% interpolate to 1 Hz
CC=C(3,1):1/86400:C(3,end);
CC(2,:)=interp1(C(3,:),C(2,:),CC); % P
CC(3,:)=interp1(C(3,:),C(1,:),CC(1,:)); % T

% write to text file, with pressure as Pa
writematrix([CC(1,:);CC(2,:)*100;CC(3,:)]',[wrtdir 'KU18-3_1Hz']) % [t P T]

% decimation loop
tc=[];
Tc=[];
c=[];
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

clearvars('temp','C','CC')

%---KU18-4---%
% read in data
fid = fopen([topdir 'apg_KU/KU18-4_184/184_125143_RAW.dat'],'r');
fspec = '%s %s %f %f'; % [date time P T]
D = textscan(fid,fspec);
fclose(fid);
t_str=cat(2,cat(1,D{1}{:}),repmat(' ',length(D{1}),1),cat(1,D{2}{:}));
tD = datenum(t_str,'yyyy/mm/dd HH:MM:SS');

% empirical trimming
tD=tD(127724:33622203);
D{3}=D{3}(127724:33622203);
D{4}=D{4}(127724:33622203);

% remove duplicate/nonmonotonic timesteps
ilist=find(diff(tD)<=0);
if ~isempty(ilist)
    tD(ilist+1)=[];
end

% interpolate to 1 Hz
DD=tD(1):1/86400:tD(end);
DD(2,:)=interp1(tD,D{3},DD); % P
DD(3,:)=interp1(tD,D{4},DD(1,:)); % T

% write to text file, with pressure as Pa
writematrix([DD(1,:);DD(2,:)*100;DD(3,:)]',[wrtdir 'KU18-4_1Hz']) % [t P T]

% decimation loop
td=[];
Td=[];
d=[];
i1 = 1;
while i1<length(D{3})
    [~,i2] = min(abs(tD-(floor(tD(i1))+1)));
    [segt,segD,~,~] = downsample_uneven(tD(i1:i2-1)',[D{3}(i1:i2-1)';D{4}(i1:i2-1)'],1/24);
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

clearvars('tD','D','DD')

%---GNS18-01---%
% read in data
fid = [topdir 'apg_GNS/out/GNS18-01.csv'];
E = readtable(fid); % [t P T]
E = [datenum(table2array(E(:,1))),table2array(E(:,2:3))]';

% empirical trimming
E(:,end-9150:end) = [];
E(:,1:12000) = [];

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
writematrix(EE',[wrtdir 'GNS18-01_1Hz']) % [t P T]

% decimation loop
te=[];
Te=[];
e=[];
i1 = 1;
d2 = floor(E(1,1)) + 1;
while i1<length(E)
    i2 = find(E(1,:)>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segE,~,~] = downsample_uneven(E(1,i1:i2-1),E(2:3,i1:i2-1),1/24);
    te=cat(2,te,segt);
    e=cat(2,e,segE(1,:));
    Te=cat(2,Te,segE(2,:));
    i1=i2;
    d2=floor(E(1,i2))+1;
end

% scale to hPa, filter tides
e=e/100;
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
Tef = interp1(te,Tef,tef);
ef = interp1(te,ef,tef);

clearvars('E','EE')

%---GNS18-03---%
% read in data
fid = [topdir 'apg_GNS/out/GNS18-03.csv'];
F = readtable(fid); % [t P T]
F = [datenum(table2array(F(:,1))),table2array(F(:,2:3))]';

% cheat that helps me not have to adjust all the indices in the trimming segment below
F2=F;
F2(:,1:9340) = [];

% remove duplicate/nonmonotonic timesteps
ilist=find(diff(F2(1,:))<=0);
if ~isempty(ilist)
    F2(:,ilist+1)=[];
end

% interpolate to 1 Hz
FF=F2(1,1):1/86400:F2(1,end);
FF(2,:)=interp1(F2(1,:),F2(2,:),FF); % P
FF(3,:)=interp1(F2(1,:),F2(3,:),FF(1,:)); % T

% write to text file, with pressure as Pa
writematrix(FF',[wrtdir 'GNS18-03_1Hz']) % [t P T]

% empirical trimming (annoying, but pays off)
F(2,2.3694e6:2.3771e6) = linspace(F(2,2.3694e6),F(2,2.3771e6),2.3771e6-2.3694e6+1);
F(2,2.4965e6:2.5045e6) = linspace(F(2,2.4965e6),F(2,2.5045e6),2.5045e6-2.4965e6+1);
F(2,2.5062e6:2.5112e6) = linspace(F(2,2.5062e6),F(2,2.5112e6),2.5112e6-2.5062e6+1);
F(2,2.517078e6:2.517148e6) = linspace(F(2,2.517078e6),F(2,2.517148e6),2.517148e6-2.517078e6+1);
F(2,2.519828e6:2.519873e6) = linspace(F(2,2.519828e6),F(2,2.519873e6),2.519873e6-2.519828e6+1);
F(2,2.5301e6:2.5306e6) = linspace(F(2,2.5301e6),F(2,2.5306e6),2.5306e6-2.5301e6+1);
F(2,2.5315e6:2.5365e6) = linspace(F(2,2.5315e6),F(2,2.5365e6),2.5365e6-2.5315e6+1);
F(2,2.78e6:2.7915e6) = linspace(F(2,2.78e6),F(2,2.7915e6),2.7915e6-2.78e6+1);
F(2,2.7995e6:2.8045e6) = linspace(F(2,2.7995e6),F(2,2.8045e6),2.8045e6-2.7995e6+1);
F(2,3.5859e6:3.58691e6) = linspace(F(2,3.5859e6),F(2,3.58691e6),3.58691e6-3.5859e6+1);
F(2,3.5894e6:3.5912e6) = linspace(F(2,3.5894e6),F(2,3.5912e6),3.5912e6-3.5894e6+1);
F(2,3.5917e6:3.5925e6) = linspace(F(2,3.5917e6),F(2,3.5925e6),3.5925e6-3.5917e6+1);
F(2,3.835e6:3.8412e6) = linspace(F(2,3.835e6),F(2,3.8412e6),3.8412e6-3.835e6+1);
F(2,3.8632e6:3.8673e6) = linspace(F(2,3.8632e6),F(2,3.8673e6),3.8673e6-3.8632e6+1);
F(2,3.9268e6:3.9305e6) = linspace(F(2,3.9268e6),F(2,3.9305e6),3.9305e6-3.9268e6+1);
F(2,4.1103e6:4.1135e6) = linspace(F(2,4.1103e6),F(2,4.1135e6),4.1135e6-4.1103e6+1);
F(2,2.65615e7:2.65631e7) = linspace(F(2,2.65615e7),F(2,2.65631e7),2.65631e7-2.65615e7+1);
F(3,2.65615e7:2.65631e7) = linspace(F(3,2.65615e7),F(3,2.65631e7),2.65631e7-2.65615e7+1); % T
F(2,2.7458e7:2.7474e7) = linspace(F(2,2.7458e7),F(2,2.7474e7),2.7474e7-2.7458e7+1);
F(3,2.7458e7:2.7474e7) = linspace(F(3,2.7458e7),F(3,2.7474e7),2.7474e7-2.7458e7+1); % T
F(2,3.1851e7:3.18525e7) = linspace(F(2,3.1851e7),F(2,3.18525e7),3.18525e7-3.1851e7+1);
F(3,3.1851e7:3.18525e7) = linspace(F(3,3.1851e7),F(3,3.18525e7),3.18525e7-3.1851e7+1); % T
F(2,2.632e7:2.641e7) = median(F(2,2.58e7:2.69e7)); % gap too large to linearly interpolate
F(3,2.632e7:2.641e7) = median(F(3,2.58e7:2.69e7)); % T
F(:,1:9340) = [];

% decimation loop
tf=[];
Tf=[];
f=[];
i1 = 1;
d2 = floor(F(1,1)) + 1;
while i1<length(F)
    i2 = find(F(1,:)>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segF,~,~] = downsample_uneven(F(1,i1:i2-1),F(2:3,i1:i2-1),1/24);
    tf=cat(2,tf,segt);
    f=cat(2,f,segF(1,:));
    Tf=cat(2,Tf,segF(2,:));
    i1=i2;
    d2=floor(F(1,i2))+1;
end

% scale to hPa, filter tides
f=f/100;
ff=Z_godin(f);
Tff=Z_godin(Tf);

% remove NaNs from tidal filter
tf(isnan(ff))=[];
Tf(isnan(ff))=[];
Tff(isnan(ff))=[];
f(isnan(ff))=[];
ff(isnan(ff))=[];

% interpolate onto monotonic time basis
tff = tf(1)+datenum(0,0,0,0,30,0):1/24:tf(end)-datenum(0,0,0,0,30,0);
Tff = interp1(tf,Tff,tff);
ff = interp1(tf,ff,tff);

clearvars('F','F2','FF')

%---GNS18-05---%
% read in data
%%%%%-----This data file will need to be updated once clock drift issues sorted out-----%%%%%
fid = [topdir 'apg_GNS/out/GNS18-05_nosync.csv'];
%%%%%-----This data file will need to be updated once clock drift issues sorted out-----%%%%%
G = readtable(fid); % [t P T]
G = [datenum(2018,10,12,01,55,00)+table2array(G(:,1))/60/60/24,table2array(G(:,2:3))]';

% empirical trimming (instrument flatlines midway through deployment)
G=G(:,9820:14092570);

% remove duplicate/nonmonotonic timesteps
ilist=find(diff(G(1,:))<=0);
if ~isempty(ilist)
    G(:,ilist+1)=[];
end

% interpolate to 1 Hz
GG=G(1,1):1/86400:G(1,end);
GG(2,:)=interp1(G(1,:),G(2,:),GG); % P
GG(3,:)=interp1(G(1,:),G(3,:),GG(1,:)); % T

% write to text file, with pressure as Pa
writematrix(GG',[wrtdir 'GNS18-05_1Hz']) % [t P T]

% decimation loop
tg=[];
Tg=[];
g=[];
i1 = 1;
d2 = floor(G(1,1)) + 1;
while i1<length(G)
    i2 = find(G(1,:)>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segG,~,~] = downsample_uneven(G(1,i1:i2-1),G(2:3,i1:i2-1),1/24);
    tg=cat(2,tg,segt);
    g=cat(2,g,segG(1,:));
    Tg=cat(2,Tg,segG(2,:));
    i1=i2;
    d2=floor(G(1,i2))+1;
end

% scale to hPa, filter tides
g=g/100;
gf=Z_godin(g);
Tgf=Z_godin(Tg);

% remove NaNs from tidal filter
tg(isnan(gf))=[];
Tg(isnan(gf))=[];
Tgf(isnan(gf))=[];
g(isnan(gf))=[];
gf(isnan(gf))=[];

% interpolate onto monotonic time basis
tgf = tg(1)+datenum(0,0,0,0,30,0):1/24:tg(end)-datenum(0,0,0,0,30,0);
Tgf = interp1(tg,Tgf,tgf);
gf = interp1(tg,gf,tgf);

clearvars('G','GG')

%---LBPR18-05---%
% read in data
fid = [topdir 'apg_LDEO/out/LBPR-05.csv'];
H = readtable(fid); % [t P T]
H = [datenum(table2array(H(:,1))),table2array(H(:,2:3))]';

% empirical trimming
H=H(:,84200:31797773);

% remove duplicate/nonmonotonic timesteps
ilist=find(diff(H(1,:))<=0);
if ~isempty(ilist)
    H(:,ilist+1)=[];
end

% interpolate to 1 Hz
HH=H(1,1):1/86400:H(1,end);
HH(2,:)=interp1(H(1,:),H(2,:),HH); % P
HH(3,:)=interp1(H(1,:),H(3,:),HH(1,:)); % T

% write to text file, with pressure as Pa
writematrix(HH',[wrtdir 'LBPR18-05_1Hz']) % [t P T]

% decimation loop
th=[];
Th=[];
h=[];
i1 = 1;
d2 = floor(H(1,1)) + 1;
while i1<length(H)
    i2 = find(H(1,:)>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segH,~,~] = downsample_uneven(H(1,i1:i2-1),H(2:3,i1:i2-1),1/24);
    th=cat(2,th,segt);
    h=cat(2,h,segH(1,:));
    Th=cat(2,Th,segH(2,:));
    i1=i2;
    d2=floor(H(1,i2))+1;
end

% scale to hPa, filter tides
h=h/100;
hf=Z_godin(h);
Thf=Z_godin(Th);

% remove NaNs from tidal filter
th(isnan(hf))=[];
Th(isnan(hf))=[];
Thf(isnan(hf))=[];
h(isnan(hf))=[];
hf(isnan(hf))=[];

% interpolate onto monotonic time basis
thf = th(1)+datenum(0,0,0,0,30,0):1/24:th(end)-datenum(0,0,0,0,30,0);
Thf = interp1(th,Thf,thf);
hf = interp1(th,hf,thf);

clearvars('H','HH')

%---GNS18-07---%
% read in data
fid = [topdir 'apg_GNS/out/GNS18-07.csv'];
J = readtable(fid); % [t P T]
J = [datenum(table2array(J(:,1))),table2array(J(:,2:3))]';

% empirical trimming
J=J(:,496907:14092570);

% remove duplicate/nonmonotonic timesteps
ilist=find(diff(J(1,:))<=0);
if ~isempty(ilist)
    J(:,ilist+1)=[];
end

% interpolate to 1 Hz
JJ=J(1,1):1/86400:J(1,end);
JJ(2,:)=interp1(J(1,:),J(2,:),JJ); % P
JJ(3,:)=interp1(J(1,:),J(3,:),JJ(1,:)); % T

% write to text file, with pressure as Pa
writematrix(JJ',[wrtdir 'GNS18-07_1Hz']) % [t P T]

% decimation loop
tj=[];
Tj=[];
j=[];
i1 = 1;
d2 = floor(J(1,1)) + 1;
while i1<length(J)
    i2 = find(J(1,:)>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segJ,~,~] = downsample_uneven(J(1,i1:i2-1),J(2:3,i1:i2-1),1/24);
    tj=cat(2,tj,segt);
    j=cat(2,j,segJ(1,:));
    Tj=cat(2,Tj,segJ(2,:));
    i1=i2;
    d2=floor(J(1,i2))+1;
end

% scale to hPa, filter tides
j=j/100;
jf=Z_godin(j);
Tjf=Z_godin(Tj);

% remove NaNs from tidal filter
tj(isnan(jf))=[];
Tj(isnan(jf))=[];
Tjf(isnan(jf))=[];
j(isnan(jf))=[];
jf(isnan(jf))=[];

% interpolate onto monotonic time basis
tjf = tj(1)+datenum(0,0,0,0,30,0):1/24:tj(end)-datenum(0,0,0,0,30,0);
Tjf = interp1(tj,Tjf,tjf);
jf = interp1(tj,jf,tjf);

clearvars('J','JJ')

%---LBPR18-04---%
% read in data
fid = [topdir 'apg_LDEO/out/LBPR-04.csv'];
K = readtable(fid); % [t P T]
K = [datenum(table2array(K(:,1))),table2array(K(:,2:3))]';

% remove duplicate/nonmonotonic timesteps
ilist=find(diff(K(1,:))<=0);
if ~isempty(ilist)
    K(:,ilist+1)=[];
end

% interpolate to 1 Hz
KK=K(1,1):1/86400:K(1,end);
KK(2,:)=interp1(K(1,:),K(2,:),KK); % P
KK(3,:)=interp1(K(1,:),K(3,:),KK(1,:)); % T

% write to text file, with pressure as Pa
writematrix(KK',[wrtdir 'LBPR18-04_1Hz']) % [t P T]

% decimation loop
tk=[];
Tk=[];
k=[];
i1 = 1;
d2 = floor(K(1,1)) + 1;
while i1<length(K)
    i2 = find(K(1,:)>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segK,~,~] = downsample_uneven(K(1,i1:i2-1),K(2:3,i1:i2-1),1/24);
    tk=cat(2,tk,segt);
    k=cat(2,k,segK(1,:));
    Tk=cat(2,Tk,segK(2,:));
    i1=i2;
    d2=floor(K(1,i2))+1;
end

% scale to hPa, filter tides
k=k/100;
kf=Z_godin(k);
Tkf=Z_godin(Tk);

% remove NaNs from tidal filter
tk(isnan(kf))=[];
Tk(isnan(kf))=[];
Tkf(isnan(kf))=[];
k(isnan(kf))=[];
kf(isnan(kf))=[];

% interpolate onto monotonic time basis
tkf = tk(1)+datenum(0,0,0,0,30,0):1/24:tk(end)-datenum(0,0,0,0,30,0);
Tkf = interp1(tk,Tkf,tkf);
kf = interp1(tk,kf,tkf);

clearvars('K','KK')

%-----PLOTTING-----%
figure(1); clf; hold on
af_plot = detrend(af);
plot(taf,af_plot,'linewidth',1)
text(taf(end)+10,af_plot(end-100),{staname{1};[num2str(stadepth(1)) ' m']},'fontsize',14)
hf_plot = detrend(hf);
plot(thf,hf_plot+20,'linewidth',1)
text(thf(end)+10,hf_plot(end-100)+20,{staname{8};[num2str(stadepth(8)) ' m']},'fontsize',14)
ef_plot = detrend(ef);
plot(tef,ef_plot+40,'linewidth',1)
text(tef(end)+10,ef_plot(end-100)+40,{staname{5};[num2str(stadepth(5)) ' m']},'fontsize',14)
gf_plot = detrend(gf);
plot(tgf,gf_plot+60,'linewidth',1)
text(tgf(end)+10,gf_plot(end-100)+60,{staname{7};[num2str(stadepth(7)) ' m']},'fontsize',14)
kf_plot = detrend(kf);
plot(tkf,kf_plot+80,'linewidth',1)
text(tkf(end)+10,kf_plot(end-100)+80,{staname{10};[num2str(stadepth(10)) ' m']},'fontsize',14)
bf_plot = detrend(bf);
plot(tbf,bf_plot+100,'linewidth',1)
text(tbf(end)+10,bf_plot(end-100)+100,{staname{2};[num2str(stadepth(2)) ' m']},'fontsize',14)
jf_plot = detrend(jf);
plot(tjf,jf_plot+120,'linewidth',1)
text(tjf(end)+10,jf_plot(end-100)+120,{staname{9};[num2str(stadepth(9)) ' m']},'fontsize',14)
cf_plot = detrend(cf);
plot(tcf,cf_plot+140,'linewidth',1)
text(tcf(end)+10,cf_plot(end-100)+140,{staname{3};[num2str(stadepth(3)) ' m']},'fontsize',14)
df_plot = detrend(df);
plot(tdf,df_plot+160,'linewidth',1)
text(tdf(end)+10,df_plot(end-100)+160,{staname{4};[num2str(stadepth(4)) ' m']},'fontsize',14)
ff_plot = detrend(ff);
plot(tff,ff_plot+180,'linewidth',1)
text(tff(end)+10,ff_plot(end-100)+180,{staname{6};[num2str(stadepth(6)) ' m']},'fontsize',14)
xlim([datenum(2018,10,01) datenum(2020,01,01)])
datetick('x','keeplimits')
set(gca,'fontsize',14)
ylabel('\DeltaP (hPa)')
title('HOBITSS 2018-2019')
box on; grid on

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print(figdir,'-dpng','-r300')

% combine and save data
t={ta,tb,tc,td,te,tf,tg,th,tj,tk};
p={a,b,c,d,e,f,g,h,j,k};
T={Ta,Tb,Tc,Td,Te,Tf,Tg,Th,Tj,Tk};
tf={taf,tbf,tcf,tdf,tef,tff,tgf,thf,tjf,tkf};
pf={af,bf,cf,df,ef,ff,gf,hf,jf,kf};
Tf={Taf,Tbf,Tcf,Tdf,Tef,Tff,Tgf,Thf,Tjf,Tkf};
save(svdir,'t','p','T','tf','pf','Tf','staname','stalat','stalon','stadepth')

% write to text files
for i=1:length(t)
    % decimated
    writematrix([t{i}',p{i}',T{i}'],[wrtdir staname{i} '_1hr'])
    % decimated, detided
    writematrix([tf{i}',pf{i}',Tf{i}'],[wrtdir staname{i} '_1hr_detided'])
end