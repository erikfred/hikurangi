% import_H7.m
%
% Read HOBITSS VII (2020-2021) data from text/csv files
%

clear; close all

%% HOBITSS VII (2020-2021)

topdir='/Volumes/Gorgoroth/apg_database/original/2020-2021_HOBITSS-VII/';
wrtdir='/Volumes/Gorgoroth/apg_database/processed/2020-2021_HOBITSS-VII/';
figdir='../figures/exploratory/HOBITSS_VII/hobitss7_stack';
svdir='../processed_data/HOBITSS_VII';

% station info
staname={'TU20-PN','TU20-PO','TU20-PP','GNS20-PF','GNS20-PH','GNS20-PC','GNS20-PD',...
    'GNS20-PE','GNS20-PG','GNS20-PI','GNS20-PK','GNS20-PL'};
stalat=[-40.83932,-40.40646,-40.24851,-40.71846,-40.61130,-40.94346,-40.89617,...
    -40.69812,-40.71058,-40.58172,-40.44469,-40.27102];
stalon=[177.64973,177.56482,177.52576,177.07040,177.48167,177.11340,177.29208,...
    177.33006,177.60205,177.16853,177.76283,177.69919];
stadepth=[2417,1877,1925,1475,1877,2020,2022,...
    1991,2198,1642,2037,1969];

%---TU20-PM---%
% Skip bc it's only a few days and data is sparse even within those days

%---TU20-PN---%
% read in data
fid = fopen([topdir 'Tohoku_Kyoto_OBP_data/TU20_PN_126897_ASCII.dat'],'r');
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
B = B(:,235512:13734302);

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
writematrix([BB(1,:);BB(2,:)*100;BB(3,:)]',[wrtdir 'TU20-PN_1Hz']) % [t P T]

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

%---TU20-PO---%
% read in data
fid = fopen([topdir 'Tohoku_Kyoto_OBP_data/TU20_PO_126896_ASCII.dat'],'r');
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
C = C(:,201175:13775717);

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
writematrix([CC(1,:);CC(2,:)*100;CC(3,:)]',[wrtdir 'TU20-PO_1Hz']) % [t P T]

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
Tc(isnan(cf))=[];
Tcf(isnan(cf))=[];
c(isnan(cf))=[];
cf(isnan(cf))=[];

% interpolate onto monotonic time basis
tcf = tc(1)+datenum(0,0,0,0,30,0):1/24:tc(end)-datenum(0,0,0,0,30,0);
cf = interp1(tc,cf,tcf);
Tcf = interp1(tc,Tcf,tcf);

clearvars('temp','C','CC')

%---TU20-PP---%
% read in data
fid = fopen([topdir 'Tohoku_Kyoto_OBP_data/TU20_PP_125143_ASCII.dat'],'r');
fspec = '%s %s %f %f'; % [date time P T]
D = textscan(fid,fspec);
fclose(fid);
t_str=cat(2,cat(1,D{1}{1:end}),repmat(' ',length(D{1}),1),cat(1,D{2}{1:end}));
tD = datenum(t_str,'yyyy/mm/dd HH:MM:SS');
D = [tD,D{3},D{4}]'; % [t P T]

% empirical trimming
D = D(:,380984:27583911);

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
writematrix([DD(1,:);DD(2,:)*100;DD(3,:)]',[wrtdir 'TU20-PP_1Hz']) % [t P T]

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

%---GNS20-PF---%
% read in data
fid = [topdir 'TN415_2023/out/GNS20-PF.csv'];
F = readtable(fid); % [t P T]
F = [datenum(table2array(F(:,1))),table2array(F(:,2:3))]';

% empirical trimming
F=F(:,46195:22049359);

% remove duplicate/nonmonotonic timesteps
ilist=find(diff(F(1,:))<=0);
if ~isempty(ilist)
    F(:,ilist+1)=[];
end

% interpolate to 1 Hz
FF=F(1,1):1/86400:F(1,end);
FF(2,:)=interp1(F(1,:),F(2,:),FF); % P
FF(3,:)=interp1(F(1,:),F(3,:),FF(1,:)); % T

% write to text file, with pressure as Pa
writematrix(FF',[wrtdir 'GNS20-PF_1Hz']) % [t P T]

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

% scale to hPa, remove anomaly at end, filter tides
tf=tf(1:5797);
f=f(1:5797)/100;
Tf=Tf(1:5797);
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

clearvars('F','FF')

%---GNS20-PH---%
% read in data
fid = [topdir 'TN415_2023/out/GNS20-PH.csv'];
H = readtable(fid); % [t P T]
H = [datenum(table2array(H(:,1))),table2array(H(:,2:3))]';

% empirical trimming
H=H(:,8089:end);

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
writematrix(HH',[wrtdir 'GNS20-PH_1Hz']) % [t P T]

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

%---GNS20-PC---%
% read in data
fid = [topdir 'TN415_2023/out/GNS20-PC.csv'];
M = readtable(fid); % [t P T]
M = [datenum(table2array(M(:,1))),table2array(M(:,2:3))]';

% empirical trimming
M=M(:,17385:26856856);

% remove duplicate/nonmonotonic timesteps
ilist=find(diff(M(1,:))<=0);
if ~isempty(ilist)
    M(:,ilist+1)=[];
end

% interpolate to 1 Hz
MM=M(1,1):1/86400:M(1,end);
MM(2,:)=interp1(M(1,:),M(2,:),MM); % P
MM(3,:)=interp1(M(1,:),M(3,:),MM(1,:)); % T

% write to text file, with pressure as Pa
writematrix(MM',[wrtdir 'GNS20-PC_1Hz']) % [t P T]

% decimation loop
tm=[];
Tm=[];
m=[];
i1 = 1;
d2 = floor(M(1,1)) + 1;
while i1<length(M)
    i2 = find(M(1,:)>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segM,~,~] = downsample_uneven(M(1,i1:i2-1),M(2:3,i1:i2-1),1/24);
    tm=cat(2,tm,segt);
    m=cat(2,m,segM(1,:));
    Tm=cat(2,Tm,segM(2,:));
    i1=i2;
    d2=floor(M(1,i2))+1;
end

% scale to hPa, filter tides
m=m/100;
mf=Z_godin(m);
Tmf=Z_godin(Tm);

% remove NaNs from tidal filter
tm(isnan(mf))=[];
Tm(isnan(mf))=[];
Tmf(isnan(mf))=[];
m(isnan(mf))=[];
mf(isnan(mf))=[];

% interpolate onto monotonic time basis
tmf = tm(1)+datenum(0,0,0,0,30,0):1/24:tm(end)-datenum(0,0,0,0,30,0);
Tmf = interp1(tm,Tmf,tmf);
mf = interp1(tm,mf,tmf);

clearvars('M','MM')

%---GNS20-PD---%
% read in data
fid = [topdir 'TN415_2023/out/GNS20-PD.csv'];
N = readtable(fid); % [t P T]
N = [datenum(table2array(N(:,1))),table2array(N(:,2:3))]';

% empirical trimming
N=N(:,30000:26820350);

% remove duplicate/nonmonotonic timesteps
ilist=find(diff(N(1,:))<=0);
if ~isempty(ilist)
    N(:,ilist+1)=[];
end

% interpolate to 1 Hz
NN=N(1,1):1/86400:N(1,end);
NN(2,:)=interp1(N(1,:),N(2,:),NN); % P
NN(3,:)=interp1(N(1,:),N(3,:),NN(1,:)); % T

% write to text file, with pressure as Pa
writematrix(NN',[wrtdir 'GNS20-PD_1Hz']) % [t P T]

% decimation loop
tn=[];
Tn=[];
n=[];
i1 = 1;
d2 = floor(N(1,1)) + 1;
while i1<length(N)
    i2 = find(N(1,:)>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segN,~,~] = downsample_uneven(N(1,i1:i2-1),N(2:3,i1:i2-1),1/24);
    tn=cat(2,tn,segt);
    n=cat(2,n,segN(1,:));
    Tn=cat(2,Tn,segN(2,:));
    i1=i2;
    d2=floor(N(1,i2))+1;
end

% scale to hPa, filter tides
n=n/100;
nf=Z_godin(n);
Tnf=Z_godin(Tn);

% remove NaNs from tidal filter
tn(isnan(nf))=[];
Tn(isnan(nf))=[];
Tnf(isnan(nf))=[];
n(isnan(nf))=[];
nf(isnan(nf))=[];

% interpolate onto monotonic time basis
tnf = tn(1)+datenum(0,0,0,0,30,0):1/24:tn(end)-datenum(0,0,0,0,30,0);
Tnf = interp1(tn,Tnf,tnf);
nf = interp1(tn,nf,tnf);

clearvars('N','NN')

%---GNS20-PE---%
% read in data
fid = [topdir 'TN415_2023/out/GNS20-PE.csv'];
P = readtable(fid); % [t P T]
P = [datenum(table2array(P(:,1))),table2array(P(:,2:3))]';

% empirical trimming
P=P(:,18000:27027600);

% remove duplicate/nonmonotonic timesteps
ilist=find(diff(P(1,:))<=0);
if ~isempty(ilist)
    P(:,ilist+1)=[];
end

% interpolate to 1 Hz
PP=P(1,1):1/86400:P(1,end);
PP(2,:)=interp1(P(1,:),P(2,:),PP); % P
PP(3,:)=interp1(P(1,:),P(3,:),PP(1,:)); % T

% write to text file, with pressure as Pa
writematrix(PP',[wrtdir 'GNS20-PE_1Hz']) % [t P T]

% decimation loop
tp=[];
Tp=[];
p=[];
i1 = 1;
d2 = floor(P(1,1)) + 1;
while i1<length(P)
    i2 = find(P(1,:)>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segP,~,~] = downsample_uneven(P(1,i1:i2-1),P(2:3,i1:i2-1),1/24);
    tp=cat(2,tp,segt);
    p=cat(2,p,segP(1,:));
    Tp=cat(2,Tp,segP(2,:));
    i1=i2;
    d2=floor(P(1,i2))+1;
end

% scale to hPa, filter tides
p=p/100;
pf=Z_godin(p);
Tpf=Z_godin(Tp);

% remove NaNs from tidal filter
tp(isnan(pf))=[];
Tp(isnan(pf))=[];
Tpf(isnan(pf))=[];
p(isnan(pf))=[];
pf(isnan(pf))=[];

% interpolate onto monotonic time basis
tpf = tp(1)+datenum(0,0,0,0,30,0):1/24:tp(end)-datenum(0,0,0,0,30,0);
Tpf = interp1(tp,Tpf,tpf);
pf = interp1(tp,pf,tpf);

clearvars('P','PP')

%---GNS20-PG---%
% read in data
fid = [topdir 'TN415_2023/out/GNS20-PG.csv'];
G = readtable(fid); % [t P T]
G = [datenum(table2array(G(:,1))),table2array(G(:,2:3))]';

% empirical trimming
G=G(:,70000:27118600);

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
writematrix(GG',[wrtdir 'GNS20-PG_1Hz']) % [t P T]

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

%---GNS20-PI---%
% read in data
fid = [topdir 'TN415_2023/out/GNS20-PI.csv'];
J = readtable(fid); % [t P T]
J = [datenum(table2array(J(:,1))),table2array(J(:,2:3))]';

% empirical trimming
J=J(:,30000:27100000);

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
writematrix(JJ',[wrtdir 'GNS20-PI_1Hz']) % [t P T]

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

%---GNS20-PK---%
% read in data
fid = [topdir 'TN415_2023/out/GNS20-PK.csv'];
K = readtable(fid); % [t P T]
K = [datenum(table2array(K(:,1))),table2array(K(:,2:3))]';

% empirical trimming
K=K(:,240000:18222000);

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
writematrix(KK',[wrtdir 'GNS20-PK_1Hz']) % [t P T]

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

%---GNS20-PL---%
% read in data
fid = [topdir 'TN415_2023/out/GNS20-PL.csv'];
L = readtable(fid); % [t P T]
L = [datenum(table2array(L(:,1))),table2array(L(:,2:3))]';

% empirical trimming
L=L(:,13000:18132000);

% remove duplicate/nonmonotonic timesteps
ilist=find(diff(L(1,:))<=0);
if ~isempty(ilist)
    L(:,ilist+1)=[];
end

% interpolate to 1 Hz
LL=L(1,1):1/86400:L(1,end);
LL(2,:)=interp1(L(1,:),L(2,:),LL); % P
LL(3,:)=interp1(L(1,:),L(3,:),LL(1,:)); % T

% write to text file, with pressure as Pa
writematrix(LL',[wrtdir 'GNS20-PL_1Hz']) % [t P T]

% decimation loop
tl=[];
Tl=[];
l=[];
i1 = 1;
d2 = floor(L(1,1)) + 1;
while i1<length(L)
    i2 = find(L(1,:)>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segL,~,~] = downsample_uneven(L(1,i1:i2-1),L(2:3,i1:i2-1),1/24);
    tl=cat(2,tl,segt);
    l=cat(2,l,segL(1,:));
    Tl=cat(2,Tl,segL(2,:));
    i1=i2;
    d2=floor(L(1,i2))+1;
end

% scale to hPa, filter tides
l=l/100;
lf=Z_godin(l);
Tlf=Z_godin(Tl);

% remove NaNs from tidal filter
tl(isnan(lf))=[];
Tl(isnan(lf))=[];
Tlf(isnan(lf))=[];
l(isnan(lf))=[];
lf(isnan(lf))=[];

% interpolate onto monotonic time basis
tlf = tl(1)+datenum(0,0,0,0,30,0):1/24:tl(end)-datenum(0,0,0,0,30,0);
Tlf = interp1(tl,Tlf,tlf);
lf = interp1(tl,lf,tlf);

clearvars('L','LL')

%-----PLOTTING-----%
figure(2); clf; hold on
bf_plot = detrend(bf);
plot(tbf,bf_plot,'linewidth',1)
text(tbf(end)+10,bf_plot(end-50),{staname{1};[num2str(stadepth(1)) ' m']},'fontsize',14)
gf_plot = detrend(gf);
plot(tgf,gf_plot+20,'linewidth',1)
text(tgf(end)+10,gf_plot(end-50)+20,{staname{9};[num2str(stadepth(9)) ' m']},'fontsize',14)
kf_plot = detrend(kf);
plot(tkf,kf_plot+40,'linewidth',1)
text(tkf(end)+10,kf_plot(end-50)+40,{staname{11};[num2str(stadepth(11)) ' m']},'fontsize',14)
nf_plot = detrend(nf);
plot(tnf,nf_plot+60,'linewidth',1)
text(tnf(end)+10,nf_plot(end-50)+60,{staname{7};[num2str(stadepth(7)) ' m']},'fontsize',14)
mf_plot = detrend(mf);
plot(tmf,mf_plot+80,'linewidth',1)
text(tmf(end)+10,mf_plot(end-50)+80,{staname{6};[num2str(stadepth(6)) ' m']},'fontsize',14)
pf_plot = detrend(pf);
plot(tpf,pf_plot+100,'linewidth',1)
text(tpf(end)+10,pf_plot(end-50)+100,{staname{8};[num2str(stadepth(8)) ' m']},'fontsize',14)
lf_plot = detrend(lf);
plot(tlf,lf_plot+120,'linewidth',1)
text(tlf(end)+10,lf_plot(end-50)+120,{staname{12};[num2str(stadepth(12)) ' m']},'fontsize',14)
df_plot = detrend(df);
plot(tdf,df_plot+140,'linewidth',1)
text(tdf(end)+10,df_plot(end-50)+140,{staname{3};[num2str(stadepth(3)) ' m']},'fontsize',14)
cf_plot = detrend(cf);
plot(tcf,cf_plot+160,'linewidth',1)
text(tcf(end)+10,cf_plot(end-50)+160,{staname{2};[num2str(stadepth(2)) ' m']},'fontsize',14)
hf_plot = detrend(hf);
plot(thf,hf_plot+180,'linewidth',1)
text(thf(end)+10,hf_plot(end-50)+180,{staname{5};[num2str(stadepth(5)) ' m']},'fontsize',14)
jf_plot = detrend(jf);
plot(tjf,jf_plot+200,'linewidth',1)
text(tjf(end)+10,jf_plot(end-50)+200,{staname{10};[num2str(stadepth(10)) ' m']},'fontsize',14)
ff_plot = detrend(ff);
plot(tff,ff_plot+220,'linewidth',1)
text(tff(end)+10,ff_plot(end-50)+220,{staname{4};[num2str(stadepth(4)) ' m']},'fontsize',14)
xlim([datenum(2020,11,01) datenum(2021,12,01)])
datetick('x','keeplimits')
set(gca,'fontsize',14)
ylabel('\DeltaP (hPa)')
title('HOBITSS 2020-2021')
box on; grid on

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print(figdir,'-dpng','-r300')

% combine and save data
t={tb,tc,td,tf,th,tm,tn,tp,tg,tj,tk,tl};
p={b,c,d,f,h,m,n,p,g,j,k,l};
T={Tb,Tc,Td,Tf,Th,Tm,Tn,Tp,Tg,Tj,Tk,Tl};
tf={tbf,tcf,tdf,tff,thf,tmf,tnf,tpf,tgf,tjf,tkf,tlf};
pf={bf,cf,df,ff,hf,mf,nf,pf,gf,jf,kf,lf};
Tf={Tbf,Tcf,Tdf,Tff,Thf,Tmf,Tnf,Tpf,Tgf,Tjf,Tkf,Tlf};
save(svdir,'t','p','T','tf','pf','Tf','staname','stalat','stalon','stadepth')

% write to text files
for i=1:length(t)
    % decimated
    writematrix([t{i}',p{i}',T{i}'],[wrtdir staname{i} '_1hr'])
    % decimated, detided
    writematrix([tf{i}',pf{i}',Tf{i}'],[wrtdir staname{i} '_1hr_detided'])
end