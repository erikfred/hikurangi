% import_H1.m
%
% Read HOBITSS I (2014-2015) data from text/csv files
%

clear; close all

%% HOBITSS I (2014-2015)

topdir='/Volumes/Gorgoroth/apg_database/original/2014-2015_HOBITSS-I/';
wrtdir='/Volumes/Gorgoroth/apg_database/processed/2014-2015_HOBITSS-I/';
figdir='../figures/exploratory/HOBITSS_I/hobitss1_stack';
svdir='../processed_data/HOBITSS_I';

% station info
staname={'EBPR-1','EBPR-2','EBPR-3','SBPR-1','SBPR-2','SBPR-3','SBPR-4','TXBPR-1','TXBPR-2','TXBPR-5'};
stalat=[-38.74617,-38.72967,-38.69393,-38.72100,-38.84742,-38.89325,-38.90785,-38.75643,-38.71345,-38.94778];
stalon=[178.68047,178.61890,178.65133,178.89317,178.87525,178.75520,178.98265,178.99733,178.56863,178.57220];
stadepth=[982.6,1007,1025,2447,2110,1354,3460.1,3532,773,1240.1];

%---EBPR-1---%
% The EBPR data are saved by month
datadir = dir([topdir 'EBPR_P/EBPR1']);
A = [];
for i=3:length(datadir)
    % read in data
    fid = fopen([datadir(i).folder '/' datadir(i).name],'r');
    fspec = '%f %f %f %f %f %f %f'; % [YY MM DD hh mm ss P]
    atemp = textscan(fid,fspec);
    fclose(fid);
    A = cat(2,A,[datenum(2000+atemp{1},atemp{2},atemp{3},atemp{4},atemp{5},atemp{6}),atemp{7}]'); % [t P]
end

% empirical trimming
A = A(:,467212:35315916);

% remove duplicate/nonmonotonic timesteps
ilist=find(diff(A(1,:))<=0);
if ~isempty(ilist)
    A(:,ilist+1)=[];
end

% interpolate to 1 Hz
AA=A(1,1):1/86400:A(1,end);
AA(2,:)=interp1(A(1,:),A(2,:),AA); % P

% write to text file, with pressure as Pa
writematrix([AA(1,:);AA(2,:)*1000]',[wrtdir 'EBPR-1_1Hz']) % [t P]

% decimation loop
ta=[];
a=[];
i1 = 1;
d2 = floor(A(1,1))+1;
while i1<length(A)
    i2 = find(A(1,:)>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segA,~,~] = downsample_uneven(A(1,i1:i2-1),A(2,i1:i2-1),1/24);
    if length(segt)>24
        keyboard
    end
    ta=cat(2,ta,segt);
    a=cat(2,a,segA);
    i1=i2;
    d2=floor(A(1,i2))+1;
end

% convert to hPa, tidal filter
a=a*10; % [hPa]
af=Z_godin(a);

% remove NaNs from tidal filter
ta(isnan(af))=[];
a(isnan(af))=[];
af(isnan(af))=[];

% interpolate onto monotonic time basis
taf = ta(1)+datenum(0,0,0,0,30,0):1/24:ta(end)-datenum(0,0,0,0,30,0);
af = interp1(ta,af,taf);
Ta = []; % no T data
Taf = []; % no T data

clearvars('A','AA')

%---EBPR-2---%
% The EBPR data are saved by month
datadir = dir([topdir 'EBPR_P/EBPR2']);
B = [];
for i=3:length(datadir)
    % read in data
    fid = fopen([datadir(i).folder '/' datadir(i).name],'r');
    fspec = '%f %f %f %f %f %f %f'; % [YY MM DD hh mm ss P]
    btemp = textscan(fid,fspec);
    fclose(fid);
    B = cat(2,B,[datenum(2000+btemp{1},btemp{2},btemp{3},btemp{4},btemp{5},btemp{6}),btemp{7}]'); % [t P]
end

% empirical trimming
B = B(:,316733:35299107);

% remove duplicate/nonmonotonic timesteps
ilist=find(diff(B(1,:))<=0);
if ~isempty(ilist)
    B(:,ilist+1)=[];
end

% interpolate to 1 Hz
BB=B(1,1):1/86400:B(1,end);
BB(2,:)=interp1(B(1,:),B(2,:),BB); % P

% write to text file, with pressure as Pa
writematrix([BB(1,:);BB(2,:)*1000]',[wrtdir 'EBPR-2_1Hz']) % [t P]

% decimation loop
tb=[];
b=[];
i1 = 1;
d2 = floor(B(1,1))+1;
while i1<length(B)
    i2 = find(B(1,:)>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segB,~,~] = downsample_uneven(B(1,i1:i2-1),B(2,i1:i2-1),1/24);
    if length(segt)>24
        keyboard
    end
    tb=cat(2,tb,segt);
    b=cat(2,b,segB);
    i1=i2;
    d2=floor(B(1,i2))+1;
end

% convert to hPa and tidal filter
b=b*10; % [hPa]
bf=Z_godin(b);

% remove NaNs from tidal filter
tb(isnan(bf))=[];
b(isnan(bf))=[];
bf(isnan(bf))=[];

% interpolate onto monotonic time basis
tbf = tb(1)+datenum(0,0,0,0,30,0):1/24:tb(end)-datenum(0,0,0,0,30,0);
bf = interp1(tb,bf,tbf);
Tb = []; % no T data
Tbf = []; % no T data

clearvars('B','BB')

%---EBPR-3---%
% The EBPR data are saved by month
datadir = dir([topdir 'EBPR_P/EBPR3']);
C = [];
for i=3:length(datadir)
    % read in data
    fid = fopen([datadir(i).folder '/' datadir(i).name],'r');
    fspec = '%f %f %f %f %f %f %f'; % [YY MM DD hh mm ss P]
    ctemp = textscan(fid,fspec);
    fclose(fid);
    C = cat(2,C,[datenum(2000+ctemp{1},ctemp{2},ctemp{3},ctemp{4},ctemp{5},ctemp{6}),ctemp{7}]'); % [t P]
end

% empirical trimming
C = C(:,327014:35375663);

% remove duplicate/nonmonotonic timesteps
ilist=find(diff(C(1,:))<=0);
if ~isempty(ilist)
    C(:,ilist+1)=[];
end

% interpolate to 1 Hz
CC=C(1,1):1/86400:C(1,end);
CC(2,:)=interp1(C(1,:),C(2,:),CC); % P

% write to text file, with pressure as Pa
writematrix([CC(1,:);CC(2,:)*1000]',[wrtdir 'EBPR-3_1Hz']) % [t P]

% decimation loop
tc=[];
c=[];
i1 = 1;
d2 = floor(C(1,1))+1;
while i1<length(C)
    i2 = find(C(1,:)>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segC,~,~] = downsample_uneven(C(1,i1:i2-1),C(2,i1:i2-1),1/24);
    if length(segt)>24
        keyboard
    end
    tc=cat(2,tc,segt);
    c=cat(2,c,segC);
    i1=i2;
    d2=floor(C(1,i2))+1;
end

% convert to hPa and tidal filter
c=c*10; % [hPa]
cf=Z_godin(c);

% remove NaNs from tidal filter
tc(isnan(cf))=[];
c(isnan(cf))=[];
cf(isnan(cf))=[];

% interpolate onto monotonic time basis
tcf = tc(1)+datenum(0,0,0,0,30,0):1/24:tc(end)-datenum(0,0,0,0,30,0);
cf = interp1(tc,cf,tcf);
Tc = []; % no T data
Tcf = []; % no T data

clearvars('C','CC')

%---SBPR-1---%
% read in data
fid = fopen([topdir 'SBPR_data/SBPR-1_hobitss2014.txt'],'r');
fspec = '%f %f %f'; % [P T t]
D = textscan(fid,fspec);
fclose(fid);
D = [D{3},D{2},D{1}]'; % [t P T]

% time is as yymmddHHMMSS
temp.base = D(1,:);
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
D(1,:) = datenum(2000+temp.yr,temp.mnth,temp.dy,temp.hr,temp.min,temp.sec);

% empirical trimming
D = D(:,2784:17645688);

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
writematrix([DD(1,:);DD(2,:)*100;DD(3,:)]',[wrtdir 'SBPR-1_1Hz']) % [t P T]

% decimation loop
td=[];
d=[];
Td=[];
i1 = 1;
d2 = floor(D(1,1))+1;
while i1<length(D)
    i2 = find(D(1,:)>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segD,~,~] = downsample_uneven(D(1,i1:i2-1),D(2:3,i1:i2-1),1/24);
    if length(segt)>24
        keyboard
    end
    td=cat(2,td,segt);
    d=cat(2,d,segD(2,:));
    Td=cat(2,Td,segD(1,:));
    i1=i2;
    d2=floor(D(1,i2))+1;
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

clearvars('temp','D','DD')

%---SBPR-2---%
% read in data
fid = fopen([topdir 'SBPR_data/SBPR-2_hobitss2014.txt'],'r');
fspec = '%f %f %f'; % [P T t]
E = textscan(fid,fspec);
fclose(fid);
E = [E{3},E{2},E{1}]'; % [t P T]

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
E = E(:,3670:17556400);

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
writematrix([EE(1,:);EE(2,:)*100;EE(3,:)]',[wrtdir 'SBPR-2_1Hz']) % [t P T]

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
    e=cat(2,e,segE(2,:));
    Te=cat(2,Te,segE(1,:));
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

%---SBPR-3---%
% read in data
fid = fopen([topdir 'SBPR_data/SBPR-3_hobitss2014.txt'],'r');
fspec = '%f %f %f'; % [P T t]
F = textscan(fid,fspec);
fclose(fid);
F = [F{3},F{2},F{1}]'; % [t P T]

% time is as yymmddHHMMSS
temp.base = F(1,:);
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
F(1,:) = datenum(2000+temp.yr,temp.mnth,temp.dy,temp.hr,temp.min,temp.sec);

% empirical trimming
F = F(:,6195:17543253);

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
writematrix([FF(1,:);FF(2,:)*100;FF(3,:)]',[wrtdir 'SBPR-3_1Hz']) % [t P T]

% decimation loop
tf=[];
f=[];
Tf=[];
i1 = 1;
d2 = floor(F(1,1))+1;
while i1<length(F)
    i2 = find(F(1,:)>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segF,~,~] = downsample_uneven(F(1,i1:i2-1),F(2:3,i1:i2-1),1/24);
    if length(segt)>24
        keyboard
    end
    tf=cat(2,tf,segt);
    f=cat(2,f,segF(2,:));
    Tf=cat(2,Tf,segF(1,:));
    i1=i2;
    d2=floor(F(1,i2))+1;
end

% tidal filter
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
ff = interp1(tf,ff,tff);
Tff = interp1(tf,Tff,tff);

clearvars('temp','F','FF')

%---SBPR-4---%
% read in data
fid = fopen([topdir 'SBPR_data/SBPR-4_hobitss2014.txt'],'r');
fspec = '%f %f %f'; % [P T t]
G = textscan(fid,fspec);
fclose(fid);
G = [G{3},G{2},G{1}]'; % [t P T]

% time is as yymmddHHMMSS
temp.base = G(1,:);
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
G(1,:) = datenum(2000+temp.yr,temp.mnth,temp.dy,temp.hr,temp.min,temp.sec);

% empirical trimming
G = G(:,5479:1530724);

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
writematrix([GG(1,:);GG(2,:)*100;GG(3,:)]',[wrtdir 'SBPR-4_1Hz']) % [t P T]

% decimation loop
tg=[];
g=[];
Tg=[];
i1 = 1;
d2 = floor(G(1,1))+1;
while i1<length(G)
    i2 = find(G(1,:)>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segG,~,~] = downsample_uneven(G(1,i1:i2-1),G(2:3,i1:i2-1),1/24);
    if length(segt)>24
        keyboard
    end
    tg=cat(2,tg,segt);
    g=cat(2,g,segG(2,:));
    Tg=cat(2,Tg,segG(1,:));
    i1=i2;
    d2=floor(G(1,i2))+1;
end

% tidal filter
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
gf = interp1(tg,gf,tgf);
Tgf = interp1(tg,Tgf,tgf);

clearvars('temp','G','GG')

%---TXBPR-1---%
% read in data
fid = fopen([topdir 'TXBPR1_1min_Pa.txt'],'r');
fspec = '%f %f'; % [t P]
H = textscan(fid,fspec);
fclose(fid);
H = [datenum(H{1})+datenum(2014,01,01),H{2}]';

% no trimming needed (empirically determined)

% remove duplicate/nonmonotonic timesteps
ilist=find(diff(H(1,:))<=0);
if ~isempty(ilist)
    H(:,ilist+1)=[];
end

% interpolate to 1 Hz
HH=H(1,1):1/86400:H(1,end);
HH(2,:)=interp1(H(1,:),H(2,:),HH); % P

% write to text file, with pressure as Pa
writematrix(HH(1,:)',[wrtdir 'TXBPR-1_1Hz']) % [t P]

% decimation loop
th=[];
h=[];
i1 = 1;
d2 = floor(H(1,1))+1;
while i1<length(H)
    i2 = find(H(1,:)>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segH,~,~] = downsample_uneven(H(1,i1:i2-1),H(2,i1:i2-1),1/24);
    if length(segt)>24
        keyboard
    end
    th=cat(2,th,segt);
    h=cat(2,h,segH);
    i1=i2;
    d2=floor(H(1,i2))+1;
end


% convert to hPa, filter tides
h=h/100;
hf=Z_godin(h);

% remove NaNs from tidal filter
th(isnan(hf))=[];
h(isnan(hf))=[];
hf(isnan(hf))=[];

% interpolate onto monotonic time basis
thf = th(1)+datenum(0,0,0,0,30,0):1/24:th(end)-datenum(0,0,0,0,30,0);
hf = interp1(th,hf,thf);
Th = []; % no T data
Thf = []; % no T data

clearvars('H','HH')

%---TXBPR-2---%
% read in data
fid = fopen([topdir 'TXBPR2_1min_Pa.txt'],'r');
fspec = '%f %f'; % [t P]
J = textscan(fid,fspec);
fclose(fid);
J = [datenum(J{1})+datenum(2014,01,01),J{2}]';

% empirical trimming (instrument failure?)
J=J(:,1:344968);

% remove duplicate/nonmonotonic timesteps
ilist=find(diff(J(1,:))<=0);
if ~isempty(ilist)
    J(:,ilist+1)=[];
end

% interpolate to 1 Hz
JJ=J(1,1):1/86400:J(1,end);
JJ(2,:)=interp1(J(1,:),J(2,:),JJ); % P

% write to text file, with pressure as Pa
writematrix(JJ(1,:)',[wrtdir 'TXBPR-2_1Hz']) % [t P]

% decimation loop
tj=[];
j=[];
i1 = 1;
d2 = floor(J(1,1))+1;
while i1<length(J)
    i2 = find(J(1,:)>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segJ,~,~] = downsample_uneven(J(1,i1:i2-1),J(2,i1:i2-1),1/24);
    if length(segt)>24
        keyboard
    end
    tj=cat(2,tj,segt);
    j=cat(2,j,segJ);
    i1=i2;
    d2=floor(J(1,i2))+1;
end

% convert to hPa, filter tides
j=j/100;
jf=Z_godin(j);

% remove NaNs from tidal filter
tj(isnan(jf))=[];
j(isnan(jf))=[];
jf(isnan(jf))=[];

% interpolate onto monotonic time basis
tjf = tj(1)+datenum(0,0,0,0,30,0):1/24:tj(end)-datenum(0,0,0,0,30,0);
jf = interp1(tj,jf,tjf);
Tj = []; % no T data
Tjf = []; % no T data

clearvars('J','JJ')

%---TXBPR-5---%
% read in data
fid = fopen([topdir 'TXBPR5_1min_Pa.txt'],'r');
fspec = '%f %f'; % [t P]
K = textscan(fid,fspec);
fclose(fid);
K = [datenum(K{1})+datenum(2014,01,01),K{2}]';

% empirical trimming (instrument failure?)
K=K(:,1:404840);

% remove duplicate/nonmonotonic timesteps
ilist=find(diff(K(1,:))<=0);
if ~isempty(ilist)
    K(:,ilist+1)=[];
end

% interpolate to 1 Hz
KK=K(1,1):1/86400:K(1,end);
KK(2,:)=interp1(K(1,:),K(2,:),KK); % P

% write to text file, with pressure as Pa
writematrix(KK(1,:)',[wrtdir 'TXBPR-5_1Hz']) % [t P]

% decimation loop
tk=[];
k=[];
i1 = 1;
d2 = floor(K(1,1))+1;
while i1<length(K)
    i2 = find(K(1,:)>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segK,~,~] = downsample_uneven(K(1,i1:i2-1),K(2,i1:i2-1),1/24);
    if length(segt)>24
        keyboard
    end
    tk=cat(2,tk,segt);
    k=cat(2,k,segK);
    i1=i2;
    d2=floor(K(1,i2))+1;
end

% convert to hPa, tidal filter
k=k/100;
kf=Z_godin(k);

% remove NaNs from tidal filter
tk(isnan(kf))=[];
k(isnan(kf))=[];
kf(isnan(kf))=[];

% interpolate onto monotonic time basis
tkf = tk(1)+datenum(0,0,0,0,30,0):1/24:tk(end)-datenum(0,0,0,0,30,0);
kf = interp1(tk,kf,tkf);
Tk = []; % no T data
Tkf = []; % no T data

clearvars('K','KK')

%-----PLOTTING-----%
figure(2); clf; hold on
hf_plot = detrend(hf);
plot(thf,hf_plot,'linewidth',1)
text(thf(end)+10,hf_plot(end-100),{staname{8};[num2str(stadepth(8)) ' m']},'fontsize',14)
gf_plot = detrend(gf);
plot(tgf,gf_plot+20,'linewidth',1)
text(tgf(end)+10,gf_plot(end-100)+20,{staname{7};[num2str(stadepth(7)) ' m']},'fontsize',14)
df_plot = detrend(df);
plot(tdf,df_plot+40,'linewidth',1)
text(tdf(end)+10,df_plot(end-100)+40,{staname{4};[num2str(stadepth(4)) ' m']},'fontsize',14)
ef_plot = detrend(ef);
plot(tef,ef_plot+60,'linewidth',1)
text(tef(end)+10,ef_plot(end-100)+60,{staname{5};[num2str(stadepth(5)) ' m']},'fontsize',14)
ff_plot = detrend(ff);
plot(tff,ff_plot+80,'linewidth',1)
text(tff(end)+10,ff_plot(end-100)+80,{staname{6};[num2str(stadepth(6)) ' m']},'fontsize',14)
kf_plot = detrend(kf);
plot(tkf,kf_plot+100,'linewidth',1)
text(tkf(end)+10,kf_plot(end-100)+100,{staname{10};[num2str(stadepth(10)) ' m']},'fontsize',14)
cf_plot = detrend(cf)*10;
plot(tcf,cf_plot+120,'linewidth',1)
text(tcf(end)+10,cf_plot(end-100)+120,{staname{3};[num2str(stadepth(3)) ' m']},'fontsize',14)
bf_plot = detrend(bf)*10;
plot(tbf,bf_plot+140,'linewidth',1)
text(tbf(end)+10,bf_plot(end-100)+140,{staname{2};[num2str(stadepth(2)) ' m']},'fontsize',14)
af_plot = detrend(af)*10;
plot(taf,af_plot+160,'linewidth',1)
text(taf(end)+10,af_plot(end-100)+160,{staname{1};[num2str(stadepth(1)) ' m']},'fontsize',14)
jf_plot = detrend(jf);
plot(tjf,jf_plot+180,'linewidth',1)
text(tjf(end)+10,jf_plot(end-100)+180,{staname{9};[num2str(stadepth(9)) ' m']},'fontsize',14)
xlim([datenum(2014,05,01) datenum(2015,09,01)])
datetick('x','keeplimits')
set(gca,'fontsize',14)
ylabel('\DeltaP (hPa)')
title('HOBITSS 2014-2015')
box on; grid on

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print(figdir,'-dpng','-r300')

% combine and save data
t={ta,tb,tc,td,te,tf,tg,th,tj,tk};
p={a,b,c,d,e,f,g,h,j,k};
tf={taf,tbf,tcf,tdf,tef,tff,tgf,thf,tjf,tkf};
pf={af,bf,cf,df,ef,ff,gf,hf,jf,kf};
T={Ta,Tb,Tc,Td,Te,Tf,Tg,Th,Tj,Tk};
Tf={Taf,Tbf,Tcf,Tdf,Tef,Tff,Tgf,Thf,Tjf,Tkf};
save(svdir,'t','p','tf','pf','T','Tf','staname','stalat','stalon','stadepth')

% write to text files
for i=1:length(t)
    % decimated
    writematrix([t{i}',p{i}',T{i}'],[wrtdir staname{i} '_1hr'])
    % decimated, detided
    writematrix([tf{i}',pf{i}',Tf{i}'],[wrtdir staname{i} '_1hr_detided'])
end