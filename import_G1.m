% import_G1.m
%
% Read GONDOR I (2022-2023) APG data from text/csv files
%

clear; close all

%% GONDOR I (2022-2023)

topdir='/Volumes/Gorgoroth/apg_database/original/2022-2023_GONDOR-I/';
wrtdir='/Volumes/Gorgoroth/apg_database/processed/2022-2023_GONDOR-I/';
figdir='../figures/exploratory/GONDOR_I/gondor1_stack';
svdir='../processed_data/GONDOR_I';

% station info
staname={'GNS22-PF2','GNS22-PH','GNS22-PI','GNS22-PJ','GNS22-PK','GNS22-PL',...
    'GNS22-PO','KU22-PA','KU22-PB','KU22-PC','KU22-PD','KU22-PE','POBS08','POBS13'};
stalat=[-38.239982,-38.739712,-38.594789,-38.827763,-38.947053,-38.955079,...
    -39.156341,-38.19342,-38.21206,-38.865517,-39.53722,-39.51798,-38.262884,-39.195760];
stalon=[179.1107,178.6855,178.831001,178.492245,178.838084,178.3662,...
    178.970941,179.03999,179.06711,178.699167,178.35901,178.47652,179.160057,178.490635];
stadepth=[2009,999,958,707,2203,1243,...
    3430,892,1365,1124,790,1610,2807,1481];

%---GNS22-PF2---%
% read in data
fid = [topdir 'GNS22-PF2.csv'];
F = readtable(fid); % [t P T]
F = [datenum(table2array(F(:,1))),table2array(F(:,2:3))]';

% empirical trimming
F = F(:,16000:23446500);

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
writematrix(FF',[wrtdir 'GNS22-PF2_1Hz']) % [t P T]

% sampled at 1 Hz, so just need to ensure even sampling
tfs=F(1,1):1/86400:F(1,end);
fs=interp1(F(1,:),F(2,:),tfs);
Tfs=interp1(F(1,:),F(3,:),tfs);

% decimation loop to 1 sample/hour
tf=[];
f=[];
Tf=[];
i1 = 1;
d2 = floor(tfs(1))+1;
while i1<length(tfs)
    i2 = find(tfs>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segF,~,~] = downsample_uneven(tfs(i1:i2-1),[fs(i1:i2-1);Tfs(i1:i2-1)],1/24);
    if length(segt)>24
        keyboard
    end
    tf=cat(2,tf,segt);
    f=cat(2,f,segF(1,:));
    Tf=cat(2,Tf,segF(2,:));
    i1=i2;
    d2=floor(tfs(i2))+1;
end

% tidal filter, convert to hPa
fs=fs/100;
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

clearvars('F','FF')

%---GNS22-PH---%
% read in data
fid = [topdir 'GNS22-PH.csv'];
H = readtable(fid); % [t P T]
H = [datenum(table2array(H(:,1))),table2array(H(:,2:3))]';

% empirical trimming
H = H(:,22000:23560140);

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
writematrix(HH',[wrtdir 'GNS22-PH_1Hz']) % [t P T]

% sampled at 1 Hz, so just need to ensure even sampling
ths=H(1,1):1/86400:H(1,end);
hs=interp1(H(1,:),H(2,:),ths);
Ths=interp1(H(1,:),H(3,:),ths);

% decimation loop to 1 sample/hour
th=[];
h=[];
Th=[];
i1 = 1;
d2 = floor(ths(1))+1;
while i1<length(ths)
    i2 = find(ths>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segH,~,~] = downsample_uneven(ths(i1:i2-1),[hs(i1:i2-1);Ths(i1:i2-1)],1/24);
    if length(segt)>24
        keyboard
    end
    th=cat(2,th,segt);
    h=cat(2,h,segH(1,:));
    Th=cat(2,Th,segH(2,:));
    i1=i2;
    d2=floor(ths(i2))+1;
end

% tidal filter, convert to hPa
hs=hs/100;
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
hf = interp1(th,hf,thf);
Thf = interp1(th,Thf,thf);

clearvars('H','HH')

%---GNS22-PI---%
% read in data
fid = [topdir 'GNS22-PI.csv'];
Z = readtable(fid); % [t P T]
Z = [datenum(table2array(Z(:,1))),table2array(Z(:,2:3))]';

% empirical trimming
Z = Z(:,11500:7278200);

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
writematrix(ZZ',[wrtdir 'GNS22-PI_1Hz']) % [t P T]

% sampled at 1 Hz, so just need to ensure even sampling
tzs=Z(1,1):1/86400:Z(1,end);
zs=interp1(Z(1,:),Z(2,:),tzs);
Tzs=interp1(Z(1,:),Z(3,:),tzs);

% decimation loop to 1 sample/hour
tz=[];
z=[];
Tz=[];
i1 = 1;
d2 = floor(tzs(1))+1;
while i1<length(tzs)
    i2 = find(tzs>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segZ,~,~] = downsample_uneven(tzs(i1:i2-1),[zs(i1:i2-1);Tzs(i1:i2-1)],1/24);
    if length(segt)>24
        keyboard
    end
    tz=cat(2,tz,segt);
    z=cat(2,z,segZ(1,:));
    Tz=cat(2,Tz,segZ(2,:));
    i1=i2;
    d2=floor(tzs(i2))+1;
end

% tidal filter, convert to hPa
zs=zs/100;
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
zf = interp1(tz,zf,tzf);
Tzf = interp1(tz,Tzf,tzf);

clearvars('Z','ZZ')

%---GNS22-PJ---%
% read in data
fid = [topdir 'GNS22-PJ.csv'];
J = readtable(fid); % [t P T]
J = [datenum(table2array(J(:,1))),table2array(J(:,2:3))]';

% empirical trimming
J = J(:,12000:7351200);

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
writematrix(JJ',[wrtdir 'GNS22-PJ_1Hz']) % [t P T]

% sampled at 1 Hz, so just need to ensure even sampling
tjs=J(1,1):1/86400:J(1,end);
js=interp1(J(1,:),J(2,:),tjs);
Tjs=interp1(J(1,:),J(3,:),tjs);

% decimation loop to 1 sample/hour
tj=[];
j=[];
Tj=[];
i1 = 1;
d2 = floor(tjs(1))+1;
while i1<length(tjs)
    i2 = find(tjs>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segJ,~,~] = downsample_uneven(tjs(i1:i2-1),[js(i1:i2-1);Tjs(i1:i2-1)],1/24);
    if length(segt)>24
        keyboard
    end
    tj=cat(2,tj,segt);
    j=cat(2,j,segJ(1,:));
    Tj=cat(2,Tj,segJ(2,:));
    i1=i2;
    d2=floor(tjs(i2))+1;
end

% tidal filter, convert to hPa
js=js/100;
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
jf = interp1(tj,jf,tjf);
Tjf = interp1(tj,Tjf,tjf);

clearvars('J','JJ')

%---GNS22-PK---%
% read in data
fid = [topdir 'GNS22-PK.csv'];
K = readtable(fid); % [t P T]
K = [datenum(table2array(K(:,1))),table2array(K(:,2:3))]';

% empirical trimming
K = K(:,15000:7520000);

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
writematrix(KK',[wrtdir 'GNS22-PK_1Hz']) % [t P T]

% sampled at 1 Hz, so just need to ensure even sampling
tks=K(1,1):1/86400:K(1,end);
ks=interp1(K(1,:),K(2,:),tks);
Tks=interp1(K(1,:),K(3,:),tks);

% decimation loop to 1 sample/hour
tk=[];
k=[];
Tk=[];
i1 = 1;
d2 = floor(tks(1))+1;
while i1<length(tks)
    i2 = find(tks>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segK,~,~] = downsample_uneven(tks(i1:i2-1),[ks(i1:i2-1);Tks(i1:i2-1)],1/24);
    if length(segt)>24
        keyboard
    end
    tk=cat(2,tk,segt);
    k=cat(2,k,segK(1,:));
    Tk=cat(2,Tk,segK(2,:));
    i1=i2;
    d2=floor(tks(i2))+1;
end

% tidal filter, convert to hPa
ks=ks/100;
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
kf = interp1(tk,kf,tkf);
Tkf = interp1(tk,Tkf,tkf);

clearvars('K','KK')

%---GNS22-PL---%
% read in data
fid = [topdir 'GNS22-PL.csv'];
L = readtable(fid); % [t P T]
L = [datenum(table2array(L(:,1))),table2array(L(:,2:3))]';

% empirical trimming
L = L(:,39000:30720000);

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
writematrix(LL',[wrtdir 'GNS22-PL_1Hz']) % [t P T]

%----------OFFSET CORRECTIONS
% updated (determined from difference with KU22-PC)
C=readtable([wrtdir 'KU22-PC_1Hz']);
C=table2array(C)';
% % offset locations in 'L' indexing
% xline(2664151);xline(2713951);xline(2750401);xline(2837651);xline(4131351)
% xline(5054151);xline(8690351);xline(25315151);xline(26618641)
% % same locations in 'C' indexing (equivalent to shared 'iA/iB' indexing)
% xline(2645000);xline(2694800);xline(2701800);xline(2731250);xline(2818500)
% xline(4112200);xline(5035000);xline(8671200);xline(25296000);xline(26599490)
% % same location in time space (equivalent for all)
% xline(C(1,2645000));xline(C(1,2694800));xline(C(1,2701800));xline(C(1,2731250))
% xline(C(1,2818500));xline(C(1,4112200));xline(C(1,5035000));xline(C(1,8671200))
% xline(C(1,25296000));xline(C(1,26599490))

L(2,2664151:2668061)=NaN;
L(2,2668061+1:end)=L(2,2668061+1:end)-35700;
% formula is [assumed L] = [C at same interval] - [net offset for C] + [net offset for L]
L(2,2664151:2668061)=C(2,2645000:2648910)-linspace(C(2,2645000),C(2,2648910),2668061-2664151+1)+...
    linspace(L(2,2664151-1),L(2,2668061+1),2668061-2664151+1);
L(2,2713951:2714851)=NaN;
L(2,2714851+1:end)=L(2,2714851+1:end)-1565;
L(2,2713951:2714851)=C(2,2694800:2695700)-linspace(C(2,2694800),C(2,2695700),2714851-2713951+1)+...
    linspace(L(2,2713951-1),L(2,2714851+1),2714851-2713951+1);
L(2,2720951:2721951)=NaN;
L(2,2721951+1:end)=L(2,2721951+1:end)-350;
L(2,2720951:2721951)=C(2,2701800:2702800)-linspace(C(2,2701800),C(2,2702800),2721951-2720951+1)+...
    linspace(L(2,2720951-1),L(2,2721951+1),2721951-2720951+1);
L(2,2750401:2753851)=NaN;
L(2,2753851+1:end)=L(2,2753851+1:end)+1280;
L(2,2750401:2753851)=C(2,2731250:2734700)-linspace(C(2,2731250),C(2,2734700),2753851-2750401+1)+...
    linspace(L(2,2750401-1),L(2,2753851+1),2753851-2750401+1);
L(2,2837651:2837691)=NaN;
L(2,2837691+1:end)=L(2,2837691+1:end)+275;
L(2,2837651:2837691)=C(2,2818500:2818540)-linspace(C(2,2818500),C(2,2818540),2837691-2837651+1)+...
    linspace(L(2,2837651-1),L(2,2837691+1),2837691-2837651+1);
L(2,4131351:4132251)=NaN;
L(2,4132251+1:end)=L(2,4132251+1:end)+300;
L(2,4131351:4132251)=C(2,4112200:4113100)-linspace(C(2,4112200),C(2,4113100),4132251-4131351+1)+...
    linspace(L(2,4131351-1),L(2,4132251+1),4132251-4131351+1);
L(2,5054151:5169151)=NaN;
L(2,5169151+1:end)=L(2,5169151+1:end)-1375;
L(2,5054151:5169151)=C(2,5035000:5150000)-linspace(C(2,5035000),C(2,5150000),5169151-5054151+1)+...
    linspace(L(2,5054151-1),L(2,5169151+1),5169151-5054151+1);
L(2,8690351:8690451)=NaN;
L(2,8690451+1:end)=L(2,8690451+1:end)+350;
L(2,8690351:8690451)=C(2,8671200:8671300)-linspace(C(2,8671200),C(2,8671300),8690451-8690351+1)+...
    linspace(L(2,8690351-1),L(2,8690451+1),8690451-8690351+1);
L(2,25315151:25317151)=NaN;
L(2,25317151+1:end)=L(2,25317151+1:end)-10675;
L(2,25315151:25317151)=C(2,25296000:25298000)-linspace(C(2,25296000),C(2,25298000),25317151-25315151+1)+...
    linspace(L(2,25315151-1),L(2,25317151+1),25317151-25315151+1);
L(2,26618641:26621151)=NaN;
L(2,26621151+1:end)=L(2,26621151+1:end)+500;
L(2,26618641:26621151)=C(2,26599490:26602000)-linspace(C(2,26599490),C(2,26602000),26621151-26618641+1)+...
    linspace(L(2,26618641-1),L(2,26621151+1),26621151-26618641+1);
%----------END CORRECTIONS

% sampled at 1 Hz, so just need to ensure even sampling
tls=L(1,1):1/86400:L(1,end);
ls=interp1(L(1,:),L(2,:),tls);
Tls=interp1(L(1,:),L(3,:),tls);

% decimation loop to 1 sample/hour
tl=[];
l=[];
Tl=[];
i1 = 1;
d2 = floor(tls(1))+1;
while i1<length(tls)
    i2 = find(tls>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segL,~,~] = downsample_uneven(tls(i1:i2-1),[ls(i1:i2-1);Tls(i1:i2-1)],1/24);
    if length(segt)>24
        keyboard
    end
    tl=cat(2,tl,segt);
    l=cat(2,l,segL(1,:));
    Tl=cat(2,Tl,segL(2,:));
    i1=i2;
    d2=floor(tls(i2))+1;
end

% tidal filter, convert to hPa
ls=ls/100;
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
lf = interp1(tl,lf,tlf);
Tlf = interp1(tl,Tlf,tlf);

clearvars('L','LL')

%---GNS22-PO---%
% read in data
fid = [topdir 'GNS22-PO.csv'];
O = readtable(fid); % [t P T]
O = [datenum(table2array(O(:,1))),table2array(O(:,2:3))]';

% empirical trimming
O = O(:,30000:7598000);

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
writematrix(OO',[wrtdir 'GNS22-PO_1Hz']) % [t P T]

% sampled at 1 Hz, so just need to ensure even sampling
tos=O(1,1):1/86400:O(1,end);
os=interp1(O(1,:),O(2,:),tos);
Tos=interp1(O(1,:),O(3,:),tos);

% decimation loop to 1 sample/hour
to=[];
o=[];
To=[];
i1 = 1;
d2 = floor(tos(1))+1;
while i1<length(tos)
    i2 = find(tos>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segO,~,~] = downsample_uneven(tos(i1:i2-1),[os(i1:i2-1);Tos(i1:i2-1)],1/24);
    if length(segt)>24
        keyboard
    end
    to=cat(2,to,segt);
    o=cat(2,o,segO(1,:));
    To=cat(2,To,segO(2,:));
    i1=i2;
    d2=floor(tos(i2))+1;
end

% tidal filter, convert to hPa
os=os/100;
o=o/100;
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

clearvars('O','OO')

%---KU22-PA---%
% read in data
fid = fopen([topdir 'KU22PA_125143.dat'],'r');
A = textscan(fid,'%s %s %f %f'); % [date time P T]
fclose(fid);
t_str=cat(2,cat(1,A{1}{:}),repmat(' ',length(A{1}),1),cat(1,A{2}{:}));
tA = datenum(t_str,'yyyy/mm/dd HH:MM:SS');
A = [tA';A{3}';A{4}']; % [t P T]

% empirical trimming
A = A(:,559400:32639400);

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
writematrix([AA(1,:);AA(2,:)*100;AA(3,:)]',[wrtdir 'KU22-PA_1Hz']) % [t P T]

% sampled at 1 Hz, so just need to ensure even sampling
tas=A(1,1):1/86400:A(1,end);
as=interp1(A(1,:),A(2,:),tas);
Tas=interp1(A(1,:),A(3,:),tas);

% decimation loop to 1 sample/hour
ta=[];
a=[];
Ta=[];
i1 = 1;
d2 = floor(tas(1))+1;
while i1<length(tas)
    i2 = find(tas>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segA,~,~] = downsample_uneven(tas(i1:i2-1),[as(i1:i2-1);Tas(i1:i2-1)],1/24);
    if length(segt)>24
        keyboard
    end
    ta=cat(2,ta,segt);
    a=cat(2,a,segA(1,:));
    Ta=cat(2,Ta,segA(2,:));
    i1=i2;
    d2=floor(tas(i2))+1;
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

clearvars('tA','A','AA')

%---KU22-PB---%
% read in data
fid = fopen([topdir 'KU22PB_126647.dat'],'r');
B = textscan(fid,'%s %s %f %f'); % [date time P T]
fclose(fid);
t_str=cat(2,cat(1,B{1}{:}),repmat(' ',length(B{1}),1),cat(1,B{2}{:}));
tB = datenum(t_str,'yyyy/mm/dd HH:MM:SS');
B = [tB';B{3}';B{4}']; % [t P T]

% empirical trimming
B = B(:,555000:end);

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
writematrix([BB(1,:);BB(2,:)*100;BB(3,:)]',[wrtdir 'KU22-PB_1Hz']) % [t P T]

% sampled at 1 Hz, so just need to ensure even sampling
tbs=B(1,1):1/86400:B(1,end);
bs=interp1(B(1,:),B(2,:),tbs);
Tbs=interp1(B(1,:),B(3,:),tbs);

% decimation loop to 1 sample/hour
tb=[];
b=[];
Tb=[];
i1 = 1;
d2 = floor(tbs(1))+1;
while i1<length(tbs)
    i2 = find(tbs>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segB,~,~] = downsample_uneven(tbs(i1:i2-1),[bs(i1:i2-1);Tbs(i1:i2-1)],1/24);
    if length(segt)>24
        keyboard
    end
    tb=cat(2,tb,segt);
    b=cat(2,b,segB(1,:));
    Tb=cat(2,Tb,segB(2,:));
    i1=i2;
    d2=floor(tbs(i2))+1;
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

% interpolate onto monotonic time basis on the hour
tbf = tb(1)+datenum(0,0,0,0,30,0):1/24:tb(end)-datenum(0,0,0,0,30,0);
bf = interp1(tb,bf,tbf);
Tbf = interp1(tb,Tbf,tbf);

clearvars('tB','B','BB')

%---KU22-PC---%
% read in data
fid = fopen([topdir 'KU22PC_126648.dat'],'r');
C = textscan(fid,'%s %s %f %f'); % [date time P T]
fclose(fid);
C{1}(end)=[]; C{2}(end)=[]; % error message as instrument was turned off
t_str=cat(2,cat(1,C{1}{:}),repmat(' ',length(C{1}),1),cat(1,C{2}{:}));
tC = datenum(t_str,'yyyy/mm/dd HH:MM:SS');
C = [tC';C{3}';C{4}']; % [t P T]

% empirical trimming
C = C(:,760000:31466800);

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
writematrix([CC(1,:);CC(2,:)*100;CC(3,:)]',[wrtdir 'KU22-PC_1Hz']) % [t P T]

% sampled at 1 Hz, so just need to ensure even sampling
tcs=C(1,1):1/86400:C(1,end);
cs=interp1(C(1,:),C(2,:),tcs);
Tcs=interp1(C(1,:),C(3,:),tcs);

% decimation loop to 1 sample/hour
tc=[];
c=[];
Tc=[];
i1 = 1;
d2 = floor(tcs(1))+1;
while i1<length(tcs)
    i2 = find(tcs>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segC,~,~] = downsample_uneven(tcs(i1:i2-1),[cs(i1:i2-1);Tcs(i1:i2-1)],1/24);
    if length(segt)>24
        keyboard
    end
    tc=cat(2,tc,segt);
    c=cat(2,c,segC(1,:));
    Tc=cat(2,Tc,segC(2,:));
    i1=i2;
    d2=floor(tcs(i2))+1;
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

% interpolate onto monotonic time basis on the hour
tcf = tc(1)+datenum(0,0,0,0,30,0):1/24:tc(end)-datenum(0,0,0,0,30,0);
cf = interp1(tc,cf,tcf);
Tcf = interp1(tc,Tcf,tcf);

clearvars('tC','C','CC')

%---KU22-PD---%
% read in data
fid = fopen([topdir 'KU22PD_126897.dat'],'r');
D = textscan(fid,'%f %f %f'); % [P T t]
fclose(fid);

% time is as yymmddHHMMSS
temp.base = D{3};
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
D{3} = datenum(2000+temp.yr,temp.mnth,temp.dy,temp.hr,temp.min,temp.sec);

D = [D{3},D{1},D{2}]'; % [t P T]

% empirical trimming
D = D(:,323000:16496500);

% remove duplicate/nonmonotonic timesteps
ilist=find(diff(D(1,:))<=0);
if ~isempty(ilist)
    D(:,ilist+1)=[];
end

% interpolate to 1 Hz
DD=D(1,1):1/86400:D(1,end);
DD(2,:)=interp1(D(1,:),D(2,:),DD); % P
DD(3,:)=interp1(D(1,:),D(3,:),DD(1,:)); % T

% look for mislabeled timestamps
ibad=find(diff(D(1,:))<0);
stepcheck=(D(1,ibad+1)-D(1,ibad))*86400<55 & (D(1,ibad+2)-D(1,ibad+1))*86400>55;
ibad=ibad(stepcheck); % ensures I only make fixes in right spots
D(1,ibad+1)=(D(1,ibad+2)+D(1,ibad))/2;

% write to text file, with pressure as Pa
writematrix([DD(1,:);DD(2,:)*100;DD(3,:)]',[wrtdir 'KU22-PD_1Hz']) % [t P T]

% sampled at 2 Hz, so need to interpolate to 1 Hz
tds=D(1,1):1/86400:D(1,end);
ds=interp1(D(1,:),D(2,:),tds);
Tds=interp1(D(1,:),D(3,:),tds);

% decimation loop to 1 sample/hour
td=[];
d=[];
Td=[];
i1 = 1;
d2 = floor(tds(1))+1;
while i1<length(tds)
    i2 = find(tds>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segD,~,~] = downsample_uneven(tds(i1:i2-1),[ds(i1:i2-1);Tds(i1:i2-1)],1/24);
    if length(segt)>24
        keyboard
    end
    td=cat(2,td,segt);
    d=cat(2,d,segD(1,:));
    Td=cat(2,Td,segD(2,:));
    i1=i2;
    d2=floor(tds(i2))+1;
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

% interpolate onto monotonic time basis on the hour
tdf = td(1)+datenum(0,0,0,0,30,0):1/24:td(end)-datenum(0,0,0,0,30,0);
df = interp1(td,df,tdf);
Tdf = interp1(td,Tdf,tdf);

clearvars('tD','D','DD')

%---KU22-PE---%
% read in data
fid = fopen([topdir 'KU22PE_126645.dat'],'r');
E = textscan(fid,'%s %s %f %f'); % [date time P T]
fclose(fid);
t_str=cat(2,cat(1,E{1}{:}),repmat(' ',length(E{1}),1),cat(1,E{2}{:}));
tE = datenum(t_str,'yyyy/mm/dd HH:MM:SS');
E = [tE';E{3}';E{4}']; % [t P T]

% empirical trimming
E = E(:,400000:32814000);

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
writematrix([EE(1,:);EE(2,:)*100;EE(3,:)]',[wrtdir 'KU22-PE_1Hz']) % [t P T]

% sampled at 1 Hz, so just need to ensure even sampling
tes=E(1,1):1/86400:E(1,end);
es=interp1(E(1,:),E(2,:),tes);
Tes=interp1(E(1,:),E(3,:),tes);

% decimation loop to 1 sample/hour
te=[];
e=[];
Te=[];
i1 = 1;
d2 = floor(tes(1))+1;
while i1<length(tes)
    i2 = find(tes>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segE,~,~] = downsample_uneven(tes(i1:i2-1),[es(i1:i2-1);Tes(i1:i2-1)],1/24);
    if length(segt)>24
        keyboard
    end
    te=cat(2,te,segt);
    e=cat(2,e,segE(1,:));
    Te=cat(2,Te,segE(2,:));
    i1=i2;
    d2=floor(tes(i2))+1;
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

% interpolate onto monotonic time basis on the hour
tef = te(1)+datenum(0,0,0,0,30,0):1/24:te(end)-datenum(0,0,0,0,30,0);
ef = interp1(te,ef,tef);
Tef = interp1(te,Tef,tef);

clearvars('tE','E','EE')

%---POBS-08---%
% read in data
fid = [topdir 'POBS8.mat'];
load(fid,'data')

for i=1:2 % consider both gauges
    if i==1
        Z=[data.ta';data.a1';data.Ta1']; % [t P T]

        % empirical trimming
        Z = Z(:,42681:end);

        % write to text file, with pressure as Pa
        writematrix([Z(1,:);Z(2,:)*100;Z(3,:)]',[wrtdir 'POBS-08_1Hz']) % [t P T]
    else
        continue

        %--- Gauge 2 seems to have some anomalies and should not be used ---%
        % Z=[data.ta';data.a2';data.Ta2']; % [t P T]
        % 
        % % empirical trimming
        % Z = Z(:,42681:end);
        % 
        % % write to text file, with pressure as Pa
        % writematrix([Z(1,:);Z(2,:)*100;Z(3,:)]','/Volumes/Gorgoroth/apg_data/DATABASE/GONDOR-I/POBS-08_G2_1Hz') % [t P T]
    end

    % sampled at 1 Hz, so just need to ensure even sampling
    tzs=Z(1,1):1/86400:Z(1,end);
    zs=interp1(Z(1,:),Z(2,:),tzs);
    Tzs=interp1(Z(1,:),Z(3,:),tzs);

    % decimation loop to 1 sample/hour
    tz=[];
    z=[];
    Tz=[];
    i1 = 1;
    d2 = floor(tzs(1))+1;
    while i1<length(tzs)
        i2 = find(tzs>=d2,1);
        if isempty(i2)
            break
        end
        [segt,segH,~,~] = downsample_uneven(tzs(i1:i2-1),[zs(i1:i2-1);Tzs(i1:i2-1)],1/24);
        if length(segt)>24
            keyboard
        end
        tz=cat(2,tz,segt);
        z=cat(2,z,segH(1,:));
        Tz=cat(2,Tz,segH(2,:));
        i1=i2;
        d2=floor(tzs(i2))+1;
    end

    % tidal filter
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
    zf = interp1(tz,zf,tzf);
    Tzf = interp1(tz,Tzf,tzf);

    if i==1
        tp8=tz; p8=z; Tp8=Tz; % decimated, unfiltered
        tp8f=tzf; p8f=zf; Tp8f=Tzf; % decimated, filtered
        tp8s=tzs; p8s=zs; Tp8s=Tzs; % 1 Hz, unfiltered
    else
        th2=tz; h2=z; Th2=Tz; % decimated, unfiltered
        th2f=tzf; h2f=zf; Th2f=Tzf; % decimated, filtered
        th2s=tzs; h2s=zs; Th2s=Tzs; % 1 Hz, unfiltered
    end
end
clearvars('Z')

%---POBS-13---%
% read in data
fid = [topdir 'POBS13.mat'];
load(fid,'data')

for i=1:2 % consider both gauges
    if i==1
        %--- The two gauges almost look like observations at different locations. Totally inobvious which is right ---%
        %--- Comparison with KU22-PE suggests Gauge 2 is at fault ---%

        Z=[data.ta';data.a1';data.Ta1']; % [t P T]

        % empirical trimming
        Z = Z(:,40625:end-5);

        % write to text file, with pressure as Pa
        writematrix([Z(1,:);Z(2,:)*100;Z(3,:)]',[wrtdir 'POBS-13_1Hz']) % [t P T]
    else
        continue
        %--- The two gauges almost look like observations at different locations. Totally inobvious which is right ---%
        %--- Comparison with KU22-PE suggests Gauge 2 is at fault ---%

        % Z=[data.ta';data.a2';data.Ta2']; % [t P T]
        % 
        % % empirical trimming
        % Z = Z(:,40625:end-5);
        % 
        % % write to text file, with pressure as Pa
        % writematrix([Z(1,:);Z(2,:)*100;Z(3,:)]','/Volumes/Gorgoroth/apg_data/DATABASE/GONDOR-I/POBS-13_G2_1Hz') % [t P T]
    end

    % sampled at 1 Hz, so just need to ensure even sampling
    tzs=Z(1,1):1/86400:Z(1,end);
    zs=interp1(Z(1,:),Z(2,:),tzs);
    Tzs=interp1(Z(1,:),Z(3,:),tzs);

    % decimation loop to 1 sample/hour
    tz=[];
    z=[];
    Tz=[];
    i1 = 1;
    d2 = floor(tzs(1))+1;
    while i1<length(tzs)
        i2 = find(tzs>=d2,1);
        if isempty(i2)
            break
        end
        [segt,segH,~,~] = downsample_uneven(tzs(i1:i2-1),[zs(i1:i2-1);Tzs(i1:i2-1)],1/24);
        if length(segt)>24
            keyboard
        end
        tz=cat(2,tz,segt);
        z=cat(2,z,segH(1,:));
        Tz=cat(2,Tz,segH(2,:));
        i1=i2;
        d2=floor(tzs(i2))+1;
    end

    % tidal filter
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
    zf = interp1(tz,zf,tzf);
    Tzf = interp1(tz,Tzf,tzf);

    if i==1
        tp13=tz; p13=z; Tp13=Tz; % decimated, unfiltered
        tp13f=tzf; p13f=zf; Tp13f=Tzf; % decimated, filtered
        tp13s=tzs; p13s=zs; Tp13s=Tzs; % 1 Hz, unfiltered
    else
        th2=tz; h2=z; Th2=Tz; % decimated, unfiltered
        th2f=tzf; h2f=zf; Th2f=Tzf; % decimated, filtered
        th2s=tzs; h2s=zs; Th2s=Tzs; % 1 Hz, unfiltered
    end
end
clearvars('Z')

%-----PLOTTING-----%
figure(3); clf; hold on
of_plot = detrend(of);
plot(tof,of_plot,'linewidth',1)
text(tof(end)+10,of_plot(end-100),{staname{7};[num2str(stadepth(7)) ' m']},'fontsize',14)
p8f_plot = detrend(p8f);
plot(tp8f,p8f_plot+10,'linewidth',1)
text(tp8f(end)+10,p8f_plot(end-100)+10,{staname{13};[num2str(stadepth(13)) ' m']},'fontsize',14)
kf_plot = detrend(kf);
plot(tkf,kf_plot+20,'linewidth',1)
text(tkf(end)+10,kf_plot(end-100)+20,{staname{5};[num2str(stadepth(5)) ' m']},'fontsize',14)
ff_plot = detrend(ff);
plot(tff,ff_plot+30,'linewidth',1)
text(tff(end)+10,ff_plot(end-100)+30,{staname{1};[num2str(stadepth(1)) ' m']},'fontsize',14)
ef_plot = detrend(ef);
plot(tef,ef_plot+40,'linewidth',1)
text(tef(end)+10,ef_plot(end-100)+40,{staname{12};[num2str(stadepth(12)) ' m']},'fontsize',14)
p13f_plot = detrend(p13f);
plot(tp13f,p13f_plot+50,'linewidth',1)
text(tp13f(end)+10,p13f_plot(end-100)+50,{staname{14};[num2str(stadepth(14)) ' m']},'fontsize',14)
bf_plot = detrend(bf);
plot(tbf,bf_plot+60,'linewidth',1)
text(tbf(end)+10,bf_plot(end-100)+60,{staname{9};[num2str(stadepth(9)) ' m']},'fontsize',14)
lf_plot = detrend(lf);
plot(tlf,lf_plot+70,'linewidth',1)
text(tlf(end)+10,lf_plot(end-100)+70,{staname{6};[num2str(stadepth(6)) ' m']},'fontsize',14)
cf_plot = detrend(cf);
plot(tcf,cf_plot+80,'linewidth',1)
text(tcf(end)+10,cf_plot(end-100)+80,{staname{10};[num2str(stadepth(10)) ' m']},'fontsize',14)
hf_plot = detrend(hf);
plot(thf,hf_plot+90,'linewidth',1)
text(thf(end)+10,hf_plot(end-100)+90,{staname{2};[num2str(stadepth(2)) ' m']},'fontsize',14)
zf_plot = detrend(zf);
plot(tzf,zf_plot+100,'linewidth',1)
text(tzf(end)+10,zf_plot(end-100)+100,{staname{3};[num2str(stadepth(3)) ' m']},'fontsize',14)
af_plot = detrend(af);
plot(taf,af_plot+110,'linewidth',1)
text(taf(end)+10,af_plot(end-100)+110,{staname{8};[num2str(stadepth(8)) ' m']},'fontsize',14)
df_plot = detrend(df);
plot(tdf,df_plot+120,'linewidth',1)
text(tdf(end)+10,df_plot(end-100)+120,{staname{11};[num2str(stadepth(11)) ' m']},'fontsize',14)
jf_plot = detrend(jf);
plot(tjf,jf_plot+130,'linewidth',1)
text(tjf(end)+10,jf_plot(end-100)+130,{staname{4};[num2str(stadepth(4)) ' m']},'fontsize',14)
xlim([datenum(2022,10,01) datenum(2024,01,01)])
datetick('x','keeplimits')
set(gca,'fontsize',14)
ylabel('\DeltaP (hPa)')
box on; grid on
title('GONDOR I 2022-2023')

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print(figdir,'-dpng','-r300')

% combine and save data
t={tf,th,tz,tj,tk,tl,to,ta,tb,tc,td,te,tp8,tp13};
p={f,h,z,j,k,l,o,a,b,c,d,e,p8,p13};
T={Tf,Th,Tz,Tj,Tk,Tl,To,Ta,Tb,Tc,Td,Te,Tp8,Tp13};
tf={tff,thf,tzf,tjf,tkf,tlf,tof,taf,tbf,tcf,tdf,tef,tp8f,tp13f};
pf={ff,hf,zf,jf,kf,lf,of,af,bf,cf,df,ef,p8f,p13f};
Tf={Tff,Thf,Tzf,Tjf,Tkf,Tlf,Tof,Taf,Tbf,Tcf,Tdf,Tef,Tp8f,Tp13f};
save(svdir,'t','p','T','tf','pf','Tf','staname','stalat','stalon','stadepth')

% write to text files
for i=1:length(t)
    % decimated
    writematrix([t{i}',p{i}',T{i}'],[wrtdir staname{i} '_1hr'])
    % decimated, detided
    writematrix([tf{i}',pf{i}',Tf{i}'],[wrtdir staname{i} '_1hr_detided'])
end

%% fixing offsets in PL data

% already incorporated in final form above
% uncomment and use only if residual problems persist

% H=readtable('../apg_data/GONDOR-I/GNS22-PH_1Hz');
% H=table2array(H)';
% L=readtable('../apg_data/GONDOR-I/GNS22-PL_1Hz');
% L=table2array(L)';
% 
% % H and L are closest in depth. Both H and F start too late to overlap with
% % first offset
% 
% [~,ia,ib]=intersect(round(H(1,:),6),round(L(1,:),6));
% HL=H(2,ia)-L(2,ib);
% HL(18010885:18010930)=NaN;
% HL(18010930:end)=HL(18010930:end)+10586;
% il1=ib(18010885);
% il2=ib(18010930);
% L(2,il1:il2)=NaN;
% L(2,il2+1:end)=L(2,il2+1:end)-10586;
% 
% L2(2,il1:il2)=linspace(L2(2,il1-1),L2(2,il2+1),il2-il1+1);
% L2(2,2667910:2668060)=NaN;
% L2(2,2668061:end)=L2(2,2668061:end)-35972;
% L2(2,2667910:2668060)=linspace(L2(2,2667910-1),L2(2,2668060+1),2668060-2667910+1);
% L2(2,2714150:2714210)=NaN;
% L2(2,2714211:end)=L2(2,2714211:end)-1304;
% L2(2,2714150:2714210)=linspace(L2(2,2714150-1),L2(2,2714210+1),2714210-2714150+1);
% L2(2,2750500:2754000)=NaN;
% L2(2,2754001:end)=L2(2,2754001:end)+1168;
% L2(2,2750500:2754000)=linspace(L2(2,2750500-1),L2(2,2754000+1),2754000-2750500+1);
% L2(2,2837658:2837678)=NaN;
% L2(2,2837679:end)=L2(2,2837679:end)+278;
% L2(2,2837658:2837678)=linspace(L2(2,2837658-1),L2(2,2837678+1),2837678-2837658+1);
% L2(2,5056000:5056600)=NaN;
% L2(2,5056601:end)=L2(2,5056601:end)-2483;
% L2(2,5056000:5056600)=linspace(L2(2,5056000-1),L2(2,5056600+1),5056600-5056000+1);
% L2(2,5079500:5082500)=NaN;
% L2(2,5082501:end)=L2(2,5082501:end)+832;
% L2(2,5079500:5081000)=L2(2,5079500-1);
% L2(2,5081000:5082500)=linspace(L2(2,5081000-1),L2(2,5082500+1),5082500-5081000+1);
% L2(2,5163000:5165500)=linspace(L2(2,5163000-1),L2(2,5165500+1),5165500-5163000+1);
% 
% % decimation loop
% tl=[];
% l=[];
% Tl=[];
% i1 = 1;
% d2 = floor(L2(1,1))+1;
% while i1<length(L2)
%     i2 = find(L2(1,:)>=d2,1);
%     if isempty(i2)
%         break
%     end
%     [segt,segF,~,~] = downsample_uneven(L2(1,i1:i2-1),L2(2:3,i1:i2-1),1/24);
%     if length(segt)>24
%         keyboard
%     end
%     tl=cat(2,tl,segt);
%     l=cat(2,l,segF(1,:));
%     Tl=cat(2,Tl,segF(2,:));
%     i1=i2;
%     d2=floor(L2(1,i2))+1;
% end
% 
% % tidal filter, convert to hPa
% l=l/100;
% lf=Z_godin(l);
% Tlf=Z_godin(Tl);
% 
% % remove NaNs from tidal filter
% tl(isnan(lf))=[];
% Tl(isnan(lf))=[];
% Tlf(isnan(lf))=[];
% l(isnan(lf))=[];
% lf(isnan(lf))=[];
% 
% % plot to check things out
% figure(11); clf; hold on
% plot(L2(1,:),L2(2,:),'.')
% plot(tl,lf*100,'linewidth',2)
% xline(L2(1,il1))
% xline(L2(1,2667910))
% xline(L2(1,2714150))
% xline(L2(1,2750500))
% xline(L2(1,2837658))
% xline(L2(1,5056000))
% xline(L2(1,5079500))
% xline(L2(1,5163000))