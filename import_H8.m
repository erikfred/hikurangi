% import_H8.m
%
% Begin exploring the Hikurangi datasets
%

clear; close all

%% HOBITSS VIII (2021-2022)

topdir='/Volumes/Gorgoroth/apg_database/original/2021-2022_HOBITSS-VIII/';
wrtdir='/Volumes/Gorgoroth/apg_database/processed/2021-2022_HOBITSS-VIII/';
figdir='../figures/exploratory/HOBITSS_VIII/hobitss8_stack';
svdir='../processed_data/HOBITSS_VIII';

% station info
staname={'GNS21-PA','GNS21-PB','GNS21-PC','GNS21-PD','GNS21-PE','GNS21-PF','GNS21-PG',...
    'GNS21-PI','GNS21-PJ','TU21-PA','TU21-PB','TU21-PC','TU21-PD','TU21-PE'};
stalat=[-39.146,-39.070,-39.00645,-38.958,-38.94795,-38.981,-38.864,...
    -38.688,-38.722,-39.66678,-39.66775,-38.89774,-38.21355,-38.19602];
stalon=[178.670,178.530,178.477,178.370,178.57155,178.800,178.670,...
    178.760,178.890,178.47045,178.54308,178.74923,179.04713,179.13581];
stadepth=[2265,1442,1406,1265,1252,1817,1159,...
    992,2461,1512,2640,1366,1194,1904];

%---GNS21-PA---%
% read in data
fid = [topdir 'GNS21-PA.csv'];
A = readtable(fid); % [t P T]
A = [datenum(table2array(A(:,1))),table2array(A(:,2:3))]';

% empirical trimming
A = A(:,23926:29974220);

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
writematrix(AA',[wrtdir 'GNS21-PA_1Hz']) % [t P T]

% decimation loop
ta=[];
Ta=[];
a=[];
i1 = 1;
d2 = floor(A(1,1)) + 1;
while i1<length(A)
    i2 = find(A(1,:)>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segA,~,~] = downsample_uneven(A(1,i1:i2-1),A(2:3,i1:i2-1),1/24);
    ta=cat(2,ta,segt);
    a=cat(2,a,segA(1,:));
    Ta=cat(2,Ta,segA(2,:));
    i1=i2;
    d2=floor(A(1,i2))+1;
end

% scale to hPa, tidal filter
a=a/100;
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

% this yields ~340 days of data
clearvars('A','AA')

%---GNS21-PB---%
% read in data
fid = [topdir 'GNS21-PB.csv'];
B = readtable(fid);
B = [datenum(table2array(B(:,1))),table2array(B(:,2:3))]';

% empirical trimming
B = B(:,22648:30000443);

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
writematrix(BB',[wrtdir 'GNS21-PB_1Hz']) % [t P T]

% decimation loop
tb=[];
Tb=[];
b=[];
i1 = 1;
d2 = floor(B(1,1)) + 1;
while i1<length(B)
    i2 = find(B(1,:)>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segB,~,~] = downsample_uneven(B(1,i1:i2-1),B(2:3,i1:i2-1),1/24);
    tb=cat(2,tb,segt);
    b=cat(2,b,segB(1,:));
    Tb=cat(2,Tb,segB(2,:));
    i1=i2;
    d2=floor(B(1,i2))+1;
end

% scale to hPa, tidal filter
b=b/100;
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

% this yields ~340 days of data
clearvars('B','BB')

%---GNS21-PC---%
% read in data
fid = [topdir 'GNS21-PC.csv'];
C = readtable(fid);
C = [datenum(table2array(C(:,1))),table2array(C(:,2:3))]';

% empirical trimming
C = C(:,16000:30019000);

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
writematrix(CC',[wrtdir 'GNS21-PC_1Hz']) % [t P T]

% anomalous period causes offset in P (but not T)
C(2,4332780:4333350) = NaN;
C(2,4333351:end) = C(2,4333351:end) + 600; % empirically determined
C(2,4332779:4333351) = linspace(C(2,4332779),C(2,4333351),4333351-4332779+1); % brief enough for linear interpolation
C(2:3,8864934:8928771) = NaN; % anomalous period (too long for linear interpolation)
C(2,8928772:end) = C(2,8928772:end) + 69;

% decimation loop
tc=[];
Tc=[];
c=[];
i1 = 1;
d2 = floor(C(1,1)) + 1;
while i1<length(C)
    i2 = find(C(1,:)>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segC,~,~] = downsample_uneven(C(1,i1:i2-1),C(2:3,i1:i2-1),1/24);
    tc=cat(2,tc,segt);
    c=cat(2,c,segC(1,:));
    Tc=cat(2,Tc,segC(2,:));
    i1=i2;
    d2=floor(C(1,i2))+1;
end

% scale to hPa, tidal filter
c=c/100;
cf=Z_godin(c);
Tcf=Z_godin(Tc);

% remove NaNs from tidal filter, and from anomaly in the middle
tc(isnan(cf))=[];
Tc(isnan(cf))=[];
Tcf(isnan(cf))=[];
c(isnan(cf))=[];
cf(isnan(cf))=[];

% interpolate onto monotonic time basis
tcf = tc(1)+datenum(0,0,0,0,30,0):1/24:tc(end)-datenum(0,0,0,0,30,0);
cf = interp1(tc,cf,tcf,'spline'); % using spline bc of NaNs in middle
Tcf = interp1(tc,Tcf,tcf,'spline');
% c = interp1(tc,c,tcf,'spline'); % using spline bc of NaNs in middle

% this yields ~335 days of data
clearvars('C','CC')

%---GNS21-PD---%
% read in data
fid = [topdir 'GNS21-PD.csv'];
D = readtable(fid);
D = [datenum(table2array(D(:,1))),table2array(D(:,2:3))]';

% empirical trimming
D = D(:,11000:30099450);

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
writematrix(DD',[wrtdir 'GNS21-PD_1Hz']) % [t P T]

% unusual behavior leads to offset
D(2,4363500:4370000) = NaN;
D(2,4370001:end) = D(2,4370001:end) + 400;
D(2,4404000:4406000) = NaN;
D(2,4406001:end) = D(2,4406001:end) - 50;
D(2,4422000:4426000) = NaN;
D(2,4426001:end) = D(2,4426001:end) - 100;
D(2,4430000:4435000) = NaN;
D(2,4435001:end) = D(2,4435001:end) - 250;
D(2,4453300:4453700) = NaN;
D(2,4453701:end) = D(2,4453701:end) + 250;
D(2,4514250:4516300) = NaN;
D(2,4516301:end) = D(2,4516301:end) - 220;
D(2,4542500:4545000) = NaN;
D(2,4545001:end) = D(2,4545001:end) + 220;
D(2,7612000:7615500) = NaN;
D(2,7615501:end) = D(2,7615501:end) + 180;

%%%%%---- manual offset correction come from the below model ----%%%%%
% generate model that fits data on BOTH sides of gap
% % try properly fitting and inverting (may cause jumpy behavior as tradeoff
% % to flatness)
% yf1=12.4206012/24;
% yf2=12/24;
% yf3=12.65834751/24;
% yf4=23.93447213/24;
% yf5=6.210300601/24;
% yf6=25.81933871/24;
% yf7=4.140200401/24;
% yf8=8.177140247/24;
% yf9=6/24;
% yf10=6.269173724/24;
% yf11=12.62600509/24;
% yf12=4/24;
% pinv8=[test(3700000:4250000)'; test(4500000:5000000)'];
% tinv8=[ttemp3(3700000:4250000)'; ttemp3(4500000:5000000)']-ttemp3(3700000);
% % construct matrix for inversion
% Gs8=[sin(2*pi/yf1*tinv8),cos(2*pi/yf1*tinv8),...
%     sin(2*pi/yf2*tinv8),cos(2*pi/yf2*tinv8),...
%     sin(2*pi/yf3*tinv8),cos(2*pi/yf3*tinv8),...
%     sin(2*pi/yf4*tinv8),cos(2*pi/yf4*tinv8),...
%     sin(2*pi/yf5*tinv8),cos(2*pi/yf5*tinv8),...
%     sin(2*pi/yf6*tinv8),cos(2*pi/yf6*tinv8),...
%     sin(2*pi/yf7*tinv8),cos(2*pi/yf7*tinv8),...
%     sin(2*pi/yf8*tinv8),cos(2*pi/yf8*tinv8),...
%     sin(2*pi/yf9*tinv8),cos(2*pi/yf9*tinv8),...
%     sin(2*pi/yf10*tinv8),cos(2*pi/yf10*tinv8),...
%     sin(2*pi/yf11*tinv8),cos(2*pi/yf11*tinv8),...
%     sin(2*pi/yf12*tinv8),cos(2*pi/yf12*tinv8),...
%     tinv8,[ones(size(tinv));zeros(size(tinv3))],[zeros(size(tinv));ones(size(tinv3))]];
% % rcond(Gs'*Gs)
% ms8=inv(Gs8'*Gs8)*Gs8'*pinv8;
% pfit8=Gs8*ms8;
% % extrapolate
% tinv9=ttemp3(3700000:5000000)'-ttemp3(3700000);
% Gs9=[sin(2*pi/yf1*tinv9),cos(2*pi/yf1*tinv9),...
%     sin(2*pi/yf2*tinv9),cos(2*pi/yf2*tinv9),...
%     sin(2*pi/yf3*tinv9),cos(2*pi/yf3*tinv9),...
%     sin(2*pi/yf4*tinv9),cos(2*pi/yf4*tinv9),...
%     sin(2*pi/yf5*tinv9),cos(2*pi/yf5*tinv9),...
%     sin(2*pi/yf6*tinv9),cos(2*pi/yf6*tinv9),...
%     sin(2*pi/yf7*tinv9),cos(2*pi/yf7*tinv9),...
%     sin(2*pi/yf8*tinv9),cos(2*pi/yf8*tinv9),...
%     sin(2*pi/yf9*tinv9),cos(2*pi/yf9*tinv9),...
%     sin(2*pi/yf10*tinv9),cos(2*pi/yf10*tinv9),...
%     sin(2*pi/yf11*tinv9),cos(2*pi/yf11*tinv9),...
%     sin(2*pi/yf12*tinv9),cos(2*pi/yf12*tinv9),...
%     tinv9,ones(size(tinv9))];
% pext9=Gs9*ms8(1:end-1);

% decimation loop
td=[];
Td=[];
d=[];
i1 = 1;
d2 = floor(D(1,1)) + 1;
while i1<length(D)
    i2 = find(D(1,:)>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segD,~,~] = downsample_uneven(D(1,i1:i2-1),D(2:3,i1:i2-1),1/24);
    td=cat(2,td,segt);
    d=cat(2,d,segD(1,:));
    Td=cat(2,Td,segD(2,:));
    i1=i2;
    d2=floor(D(1,i2))+1;
end

% scale to hPa, tidal filter
d=d/100;
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
df = interp1(td,df,tdf,'spline');
Tdf = interp1(td,Tdf,tdf);

% this yields ~342 days of data
clearvars('D','DD')

%---GNS21-PE---%
% read in data
fid = [topdir 'GNS21-PE.csv'];
E = readtable(fid);
E = [datenum(table2array(E(:,1))),table2array(E(:,2:3))]';

% empirical trimming
E = E(:,20000:21920000);
% with additional filtering, may be able to extend to:
% E = E(:,20000:29976000);
% SSE deformation seen on GPS just after existing truncation :(

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
writematrix(EE',[wrtdir 'GNS21-PE_1Hz']) % [t P T]

%---------- cleaning up high amplitude noise
B=readtable([wrtdir 'GNS21-PB_1Hz.txt']);
B=table2array(B)';
[~,iA,iB]=intersect(round(E(1,:),6),round(B(1,:),6));
EB=E(2,iA)-B(2,iB);

y1=-1940300;
y2=-1939400;
fcond1=EB>y1 & EB<y2;
fcond1(1:2.1e7)=true;
EB(~fcond1)=NaN; % inject NaNs where things are bad
EB_mf=medfilt1(EB,3600,'truncate'); % smear those NaNs to eliminate nearby bad points
fcond2=isnan(EB_mf); % get NaN indices
iE=iA(fcond2); % same indices in original basis
E(2,iE)=NaN; % NaN-ing original data array
iN=isnan(E(2,:)); %
iNd=diff(iN);
% find continuous segments of NaNs
% assume that E==B for these segments, + linear adjustment to fit endpoints

%----------

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

% scale to hPa, tidal filter
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
ef = interp1(te,ef,tef);
Tef = interp1(te,Tef,tef);

clearvars('E','EE')

%---GNS21-PF---%
% read in data
fid = [topdir 'GNS21-PF.csv'];
F = readtable(fid);
F = [datenum(table2array(F(:,1))),table2array(F(:,2:3))]';

% empirical trimming
F = F(:,14731:29979909);

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
writematrix(FF',[wrtdir 'GNS21-PF_1Hz']) % [t P T]

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

% scale to hPa, tidal filter
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
ff = interp1(tf,ff,tff);
Tff = interp1(tf,Tff,tff);

% this yields ~340 days of data
clearvars('F','FF')

%---GNS21-PG---%
% read in data
fid = [topdir 'GNS21-PG.csv'];
G = readtable(fid);
G = [datenum(table2array(G(:,1))),table2array(G(:,2:3))]';

% empirical trimming
% first ~22 days of pressure data are complete nonsense
G = G(:,1891137:28094835);

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
writematrix(GG',[wrtdir 'GNS21-PG_1Hz']) % [t P T]

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

% scale to hPa, tidal filter
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
gf = interp1(tg,gf,tgf);
Tgf = interp1(tg,Tgf,tgf);

% this yields ~320 days of data
clearvars('G','GG')

%---GNS21-PI---%
% read in data
fid = [topdir 'GNS21-PI.csv'];
Z = readtable(fid);
Z = [datenum(table2array(Z(:,1))),table2array(Z(:,2:3))]';

% empirical trimming
Z = Z(:,14718:32877925);
Z(2:3,27003466:27224753) = NaN; % anomalous period

% remove duplicate/nonmonotonic timesteps
ilist=find(diff(Z(1,:))<=0);
if ~isempty(ilist)
    Z(:,ilist+1)=[];
end

% interpolate to 1 Hz
ZZ=Z(1,1):1/86400:Z(1,end);
ZZ(2,:)=interp1(Z(1,:),Z(2,:),ZZ); % P
ZZ(3,:)=interp1(Z(1,:),Z(3,:),ZZ(1,:)); % T

% write to text file, with pressure as Pa
writematrix(ZZ',[wrtdir 'GNS21-PI_1Hz']) % [t P T]

% decimation loop
tz=[];
Tz=[];
z=[];
i1 = 1;
d2 = floor(Z(1,1)) + 1;
while i1<length(Z)
    i2 = find(Z(1,:)>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segZ,~,~] = downsample_uneven(Z(1,i1:i2-1),Z(2:3,i1:i2-1),1/24);
    tz=cat(2,tz,segt);
    z=cat(2,z,segZ(1,:));
    Tz=cat(2,Tz,segZ(2,:));
    i1=i2;
    d2=floor(Z(1,i2))+1;
end

% scale to hPa, tidal filter
z=z/100;
zf=Z_godin(z);
Tzf=Z_godin(Tz);

% remove NaNs from tidal filter
tz(isnan(zf))=[];
Tz(isnan(zf))=[];
Tzf(isnan(zf))=[];
z(isnan(zf))=[];
zf(isnan(zf))=[];

% interpolate onto monotonic time basis
tzf = tz(1)+datenum(0,0,0,0,30,0):1/24:tz(end)-datenum(0,0,0,0,30,0);
zf = interp1(tz,zf,tzf,'spline'); % again using spline bc of gap in middle
Tzf = interp1(tz,Tzf,tzf,'spline');

% this yields ~375 days of data
clearvars('Z','ZZ')

%---GNS21-PJ---%
% read in data
fid = [topdir 'GNS21-PJ.csv'];
J = readtable(fid);
J = [datenum(table2array(J(:,1))),table2array(J(:,2:3))]';

% empirical trimming
J = J(:,18936:30093788);

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
writematrix(JJ',[wrtdir 'GNS21-PJ_1Hz']) % [t P T]

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

% scale to hPa, tidal filter
j=j/100; %[hPa]
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
jf = interp1(tj,jf,tjf);
Tjf = interp1(tj,Tjf,tjf);

% this yields ~341 days of data
clearvars('J','JJ')

%---TU21-PA---%
% read in data
fid = [topdir 'TU21-PA_131477_RAW.dat'];
K = readtable(fid);
K = [table2array(K(:,3)),table2array(K(:,1:2))]'; % [t P T]

% time is as yymmddHHMMSS
temp.base = K(1,:);
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
K(1,:) = datenum(2000+temp.yr,temp.mnth,temp.dy,temp.hr,temp.min,temp.sec);

% empirical trimming
K = K(:,259218:15093619);

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
writematrix([KK(1,:);KK(2,:)*100;KK(3,:)]',[wrtdir 'TU21-PA_1Hz']) % [t P T]

% decimation loop
tk=[];
k=[];
Tk=[];
i1 = 1;
d2 = floor(K(1,1))+1;
while i1<length(K)
    i2 = find(K(1,:)>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segK,~,~] = downsample_uneven(K(1,i1:i2-1),K(2:3,i1:i2-1),1/24);
    if length(segt)>24
        keyboard
    elseif length(segt)~=length(segK)
        keyboard
    end
    tk=cat(2,tk,segt);
    k=cat(2,k,segK(1,:));
    Tk=cat(2,Tk,segK(2,:));
    i1=i2;
    d2=floor(K(1,i2))+1;
end

% tidal filter
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
kf = interp1(tk,kf,tkf);
Tkf = interp1(tk,Tkf,tkf);

clearvars('temp','K','KK')

%---TU21-PB---%
% read in data
fid = [topdir 'TU21-PB_131479_RAW.dat'];
L = readtable(fid);
L = [table2array(L(:,3)),table2array(L(:,1:2))]'; % [t P T]

% time is as yymmddHHMMSS
temp.base = L(1,:);
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
L(1,:) = datenum(2000+temp.yr,temp.mnth,temp.dy,temp.hr,temp.min,temp.sec);

% empirical trimming
L = L(:,259206:15093612);

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
writematrix([LL(1,:);LL(2,:)*100;LL(3,:)]',[wrtdir 'TU21-PB_1Hz']) % [t P T]

% decimation loop
tl=[];
l=[];
Tl=[];
i1 = 1;
d2 = floor(L(1,1))+1;
while i1<length(L)
    i2 = find(L(1,:)>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segL,~,~] = downsample_uneven(L(1,i1:i2-1),L(2:3,i1:i2-1),1/24);
    if length(segt)>24
        keyboard
    elseif length(segt)~=length(segL)
        keyboard
    end
    tl=cat(2,tl,segt);
    l=cat(2,l,segL(1,:));
    Tl=cat(2,Tl,segL(2,:));
    i1=i2;
    d2=floor(L(1,i2))+1;
end

% tidal filter
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
lf = interp1(tl,lf,tlf);
Tlf = interp1(tl,Tlf,tlf);

clearvars('temp','L','LL')

%---TU21-PC---%
% read in data
fid = [topdir 'TU21-PC_126652_RAW.dat'];
M = readtable(fid);
M = [table2array(M(:,3)),table2array(M(:,1:2))]'; % [t P T]

% time is as yymmddHHMMSS
temp.base = M(1,:);
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
M(1,:) = datenum(2000+temp.yr,temp.mnth,temp.dy,temp.hr,temp.min,temp.sec);

% empirical trimming
M = M(:,230363:15092328);

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
writematrix([MM(1,:);MM(2,:)*100;MM(3,:)]',[wrtdir 'TU21-PC_1Hz']) % [t P T]

% decimation loop
tm=[];
m=[];
Tm=[];
i1 = 1;
d2 = floor(M(1,1))+1;
while i1<length(M)
    i2 = find(M(1,:)>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segM,~,~] = downsample_uneven(M(1,i1:i2-1),M(2:3,i1:i2-1),1/24);
    if length(segt)>24
        keyboard
    elseif length(segt)~=length(segM)
        keyboard
    end
    tm=cat(2,tm,segt);
    m=cat(2,m,segM(1,:));
    Tm=cat(2,Tm,segM(2,:));
    i1=i2;
    d2=floor(M(1,i2))+1;
end

% tidal filter
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
mf = interp1(tm,mf,tmf);
Tmf = interp1(tm,Tmf,tmf);

clearvars('temp','M','MM')

%---TU21-PD---%
% read in data
fid = fopen([topdir 'TU21-PD_125142_RAW.dat']);
fspec = '%s %s %f %f'; % [date time P T]
N = textscan(fid,fspec);
fclose(fid);
t_str=cat(2,cat(1,N{1}{1:end}),repmat(' ',length(N{1}),1),cat(1,N{2}{1:end}));
tN = datenum(t_str,'yyyy/mm/dd HH:MM:SS');
N = [tN,N{3},N{4}]'; % [t P T]

% empirical trimming
N = N(:,798239:30948059);

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
writematrix([NN(1,:);NN(2,:)*100;NN(3,:)]',[wrtdir 'TU21-PD_1Hz']) % [t P T]

% decimation loop
tn=[];
n=[];
Tn=[];
i1 = 1;
d2 = floor(N(1,1))+1;
while i1<length(N)
    i2 = find(N(1,:)>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segN,~,~] = downsample_uneven(N(1,i1:i2-1),N(2:3,i1:i2-1),1/24);
    if length(segt)>24
        keyboard
    elseif length(segt)~=length(segN)
        keyboard
    end
    tn=cat(2,tn,segt);
    n=cat(2,n,segN(1,:));
    Tn=cat(2,Tn,segN(2,:));
    i1=i2;
    d2=floor(N(1,i2))+1;
end

% tidal filter
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
nf = interp1(tn,nf,tnf);
Tnf = interp1(tn,Tnf,tnf);

clearvars('tN','N','NN')

%---TU21-PE---%
% read in data
fid = fopen([topdir 'TU21-PE_126646_RAW.dat']);
fspec = '%s %s %f %f'; % [date time P T]
O = textscan(fid,fspec);
fclose(fid);
t_str=cat(2,cat(1,O{1}{1:end}),repmat(' ',length(O{1}),1),cat(1,O{2}{1:end}));
tO = datenum(t_str,'yyyy/mm/dd HH:MM:SS');
O = [tO,O{3},O{4}]'; % [t P T]

% empirical trimming
O = O(:,808770:30927830);

% remove duplicate/nonmonotonic timesteps
ilist=find(diff(O(1,:))<=0);
if ~isempty(ilist)
    O(:,ilist+1)=[];
end

% interpolate to 1 Hz
OO=O(1,1):1/86400:O(1,end);
OO(2,:)=interp1(O(1,:),O(2,:),OO); % P
OO(3,:)=interp1(O(1,:),O(3,:),OO(1,:)); % T

% write to text file, with pressure as Pa
writematrix([OO(1,:);OO(2,:)*100;OO(3,:)]',[wrtdir 'TU21-PE_1Hz']) % [t P T]

% decimation loop
to=[];
o=[];
To=[];
i1 = 1;
d2 = floor(O(1,1))+1;
while i1<length(O)
    i2 = find(O(1,:)>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segO,~,~] = downsample_uneven(O(1,i1:i2-1),O(2:3,i1:i2-1),1/24);
    if length(segt)>24
        keyboard
    elseif length(segt)~=length(segO)
        keyboard
    end
    to=cat(2,to,segt);
    o=cat(2,o,segO(1,:));
    To=cat(2,To,segO(2,:));
    i1=i2;
    d2=floor(O(1,i2))+1;
end

% tidal filter
of=Z_godin(o);
Tof=Z_godin(To);

% remove NaNs from tidal filter
to(isnan(of))=[];
To(isnan(of))=[];
Tof(isnan(of))=[];
o(isnan(of))=[];
of(isnan(of))=[];

% interpolate onto monotonic time basis
tof = to(1)+datenum(0,0,0,0,30,0):1/24:to(end)-datenum(0,0,0,0,30,0);
of = interp1(to,of,tof);
Tof = interp1(to,Tof,tof);

clearvars('tO','O','OO')

%-----PLOTTING-----%
figure(3); clf; hold on
lf_plot = detrend(lf);
plot(tlf,lf_plot-20,'linewidth',1)
text(tlf(end)+10,lf_plot(end-100)-20,{staname{11};[num2str(stadepth(11)) ' m']},'fontsize',14)
jf_plot = detrend(jf);
plot(tjf,jf_plot,'linewidth',1)
text(tjf(end)+10,jf_plot(end-100),{staname{9};[num2str(stadepth(9)) ' m']},'fontsize',14)
af_plot = detrend(af);
plot(taf,af_plot+20,'linewidth',1)
text(taf(end)+10,af_plot(end-100)+20,{staname{1};[num2str(stadepth(1)) ' m']},'fontsize',14)
of_plot = detrend(of);
plot(tof,of_plot+40,'linewidth',1)
text(tof(end)+10,of_plot(end-100)+40,{staname{14};[num2str(stadepth(14)) ' m']},'fontsize',14)
ff_plot = detrend(ff);
plot(tff,ff_plot+60,'linewidth',1)
text(tff(end)+10,ff_plot(end-100)+60,{staname{6};[num2str(stadepth(6)) ' m']},'fontsize',14)
kf_plot = detrend(kf);
plot(tkf,kf_plot+80,'linewidth',1)
text(tkf(end)+10,kf_plot(end-100)+80,{staname{10};[num2str(stadepth(10)) ' m']},'fontsize',14)
bf_plot = detrend(bf);
plot(tbf,bf_plot+100,'linewidth',1)
text(tbf(end)+10,bf_plot(end-100)+100,{staname{2};[num2str(stadepth(2)) ' m']},'fontsize',14)
cf_plot = detrend(cf);
plot(tcf,cf_plot+120,'linewidth',1)
text(tcf(end)+10,cf_plot(end-100)+120,{staname{3};[num2str(stadepth(3)) ' m']},'fontsize',14)
mf_plot = detrend(mf);
plot(tmf,mf_plot+140,'linewidth',1)
text(tmf(end)+10,mf_plot(end-100)+140,{staname{12};[num2str(stadepth(12)) ' m']},'fontsize',14)
df_plot = detrend(df);
plot(tdf,df_plot+160,'linewidth',1)
text(tdf(end)+10,df_plot(end-100)+160,{staname{4};[num2str(stadepth(4)) ' m']},'fontsize',14)
ef_plot = detrend(ef);
plot(tef,ef_plot+180,'linewidth',1)
text(tef(end)+10,ef_plot(end-100)+180,{staname{5};[num2str(stadepth(5)) ' m']},'fontsize',14)
nf_plot = detrend(nf);
plot(tnf,nf_plot+200,'linewidth',1)
text(tnf(end)+10,nf_plot(end-100)+200,{staname{13};[num2str(stadepth(13)) ' m']},'fontsize',14)
gf_plot = detrend(gf);
plot(tgf,gf_plot+220,'linewidth',1)
text(tgf(end)+10,gf_plot(end-100)+220,{staname{7};[num2str(stadepth(7)) ' m']},'fontsize',14)
zf_plot = detrend(zf);
plot(tzf,zf_plot+240,'linewidth',1)
text(tzf(end)+10,zf_plot(end-100)+240,{staname{8};[num2str(stadepth(8)) ' m']},'fontsize',14)
xlim([datenum(2021,10,01) datenum(2023,02,01)])
datetick('x','keeplimits')
set(gca,'fontsize',14)
ylabel('\DeltaP (hPa)')
box on; grid on
title('HOBITSS 2021-2022')

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print([figdir 'hobitss8_stack'],'-dpng','-r300')

% combine and save data
t={ta,tb,tc,td,te,tf,tg,tz,tj,tk,tl,tm,tn,to};
p={a,b,c,d,e,f,g,z,j,k,l,m,n,o};
T={Ta,Tb,Tc,Td,Te,Tf,Tg,Tz,Tj,Tk,Tl,Tm,Tn,To};
tf={taf,tbf,tcf,tdf,tef,tff,tgf,tzf,tjf,tkf,tlf,tmf,tnf,tof};
pf={af,bf,cf,df,ef,ff,gf,zf,jf,kf,lf,mf,nf,of};
Tf={Taf,Tbf,Tcf,Tdf,Tef,Tff,Tgf,Tzf,Tjf,Tkf,Tlf,Tmf,Tnf,Tof};
save(svdir,'t','p','T','tf','pf','Tf','staname','stalat','stalon','stadepth')

% write to text files
for i=1:length(t)
    % decimated
    writematrix([t{i}',p{i}',T{i}'],[wrtdir staname{i} '_1hr'])
    % decimated, detided
    writematrix([tf{i}',pf{i}',Tf{i}'],[wrtdir staname{i} '_1hr_detided'])
end