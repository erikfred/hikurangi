% import_CORK.m
%
% Begin exploring the Hikurangi datasets
%

clear; close all

%% CORK line

topdir='/Volumes/Gorgoroth/apg_database/original/CORK/';
wrtdir='/Volumes/Gorgoroth/apg_database/processed/CORK/';
figdir='../figures/exploratory/CORK/cork_stack';
svdir='../processed_data/CORK';

% station info
staname={'BPR','U1518','U1519'};
stalat=[0,-38.859,-38.727]; % I don't have BPR coordinates
stalon=[0,178.896,178.615]; % I don't have BPR coordinates
stadepth=[0,2631,1000]; % I don't have BPR coordinates

%---BPR2 (incomming plate) Part 1---%
% read in data
fid = fopen([topdir 'incoming_plate_BPR/BPR2_TAN2102.dat'],'r');
fspec = '%f %f %f %f %f'; % [t ? T T-related? P]
A = textscan(fid,fspec,'HeaderLines',7);
fclose(fid);
A = [A{1},A{5},A{3}]'; % [t P T]

%---BPR2 (incomming plate) Part 2---%
% read in data
fid = fopen([topdir 'incoming_plate_BPR/BPR2_TN415_2023.dat'],'r');
fspec = '%f %f %f %f %f'; % [t ? T T-related? P]
B = textscan(fid,fspec,'HeaderLines',7);
fclose(fid);
A = cat(2,A,[B{1},B{5},B{3}]'); % [t P T]

% time is as yymmddHHMMSS
temp.base = A(1,:);
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
A(1,:) = datenum(2000+temp.yr,temp.mnth,temp.dy,temp.hr,temp.min,temp.sec);

% empirical trimming
A = A(:,4070:end);

%----STOPGAP TEMPERATURE FIX----%
% Laura is looking into this at the binary level
tinv=A(1,:)-A(1,1); tinv=tinv/max(tinv); tinv=tinv';
Tinv=A(3,:)';
Pinv=A(2,:)';
Ginv=[ones(size(tinv)),tinv,Tinv];
minv=inv(Ginv'*Ginv)*Ginv'*Pinv;
Pmod=Ginv*minv;
A(2,:)=Pinv'-minv(3)*Tinv';
A(2,1078444:1078488)=linspace(A(2,1078444),A(2,1078488),1078488-1078444+1);
A(2,1080532:1080550)=linspace(A(2,1080532),A(2,1080550),1080550-1080532+1);

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
writematrix([AA(1,:);AA(2,:)*1000;AA(3,:)]',[wrtdir 'LTBPR_1Hz'])

% decimation loop
ta=[];
Ta=[];
a=[];
i1 = 1;
d2 = floor(A(1,1))+1;
while i1<length(A)
    i2 = find(A(1,:)>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segA,~,~] = downsample_uneven(A(1,i1:i2-1),A(2:3,i1:i2-1),1/24);
    if length(segt)>24
        keyboard
    end
    ta=cat(2,ta,segt);
    a=cat(2,a,segA(1,:));
    Ta=cat(2,Ta,segA(2,:));
    i1=i2;
    d2=floor(A(1,i2))+1;
end

% tidal filter, convert to hPa
a=a*10;
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

% check data length, timestamps against expectations
clearvars('temp','A','B','AA')

%---U1518 Part 1---%
% read in data
fid = fopen([topdir 'U1518/U1518_RR1902_clean.dat'],'r');
fspec = '%s %s %f %f %f %f %f %f %f %f %f'; % [t t T T1 P1 T2 P2 T3 P3 T4 P4]
E = textscan(fid,fspec,'HeaderLines',3);
fclose(fid);
t_str=cat(2,cat(1,E{1}{:}),repmat(' ',length(E{1}),1),cat(1,E{2}{:}));
tE = datenum(t_str,'yyyy-mm-dd HH:MM:SS');
E = [tE,E{11},E{10}]'; % [t P1 T1 P2 T2 P3 T3 P4 T4]

%---U1518 Part 2---%
% read in data
fid = fopen([topdir 'U1518/U1518_TAN2102.dat'],'r');
fspec = '%s %s %f %f %f %f %f %f %f %f %f'; % [t t T T1 P1 T2 P2 T3 P3 T4 P4]
F = textscan(fid,fspec);
fclose(fid);
t_str=cat(2,cat(1,F{1}{:}),repmat(' ',length(F{1}),1),cat(1,F{2}{:}));
tF = datenum(t_str,'yyyy-mm-dd HH:MM:SS');
F{11} = [F{11}(1:525835);F{11}(525836:end)-F{11}(525836)+F{11}(525835)+...
    mean(diff(F{11}(525820:525835)))]; % corrects an artificial offset
E = cat(2,E,[tF,F{11},F{10}]'); % [t P T]

%---U1518 Part 3---%
% read in data
fid = fopen([topdir 'U1518/U1518_2023_formatted.dat'],'r');
fspec = '%s %s %f %f %f %f %f %f %f %f %f'; % [t t T T1 P1 T2 P2 T3 P3 T4 P4]
K = textscan(fid,fspec);
fclose(fid);
t_str=cat(2,cat(1,K{1}{:}),repmat(' ',length(K{1}),1),cat(1,K{2}{:}));
tK = datenum(t_str,'yyyy-mm-dd HH:MM:SS');
E = cat(2,E,[tK,K{11}/10+0.15,K{10}]'); % [t P T]

% empirical trimming, convert to Pa
E = E(:,3123:end);
E(2,92703:92706)=linspace(E(2,92703),E(2,92706),4);

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
writematrix([EE(1,:);EE(2,:)*10000;EE(3,:)]',[wrtdir 'U1518_1Hz'])

% decimation loop
te=[];
Te=[];
e=[];
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

% scale to hPa, tidal filter
e=e*100; % [hPa]
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
ef = interp1(te,ef,tef,'spline'); % spline helps with gap at March 2021
Tef = interp1(te,Tef,tef);

clearvars('t_str','tC','C','tD','D','tE','E','EE','tF','F','tK','K')

%---U1519 Part 1---%
% read in data
fid = fopen([topdir 'U1519/U1519_RR1902_clean.dat'],'r');
fspec = '%s %s %f %f %f %f %f %f %f'; % [t t ? T1 P1 T2 P2 T3 P3]
G = textscan(fid,fspec,'HeaderLines',3);
fclose(fid);
t_str=cat(2,cat(1,G{1}{:}),repmat(' ',length(G{1}),1),cat(1,G{2}{:}));
tG = datenum(t_str,'yyyy-mm-dd HH:MM:SS');
G = [tG,G{9},G{8}]'; % [t P T]

%---U1519 Part 2---%
% read in data
fid = fopen([topdir 'U1519/U1519_TAN2102.dat'],'r');
fspec = '%s %s %f %f %f %f %f %f %f'; % [t t ? T1 P1 T2 P2 T3 P3]
H = textscan(fid,fspec);
fclose(fid);
t_str=cat(2,cat(1,H{1}{:}),repmat(' ',length(H{1}),1),cat(1,H{2}{:}));
tH = datenum(t_str,'yyyy-mm-dd HH:MM:SS');
G =  cat(2,G,[tH,H{9},H{8}]'); % [t P T]

%---U1519 Part 3---%
% read in data
fid = fopen([topdir 'U1519/U1519_TN415_2023.dat'],'r');
fspec = '%f %f %f %f %f %f %f %f %f'; % [t ? T T1 P1 T2 P2 T3 P3]
J = textscan(fid,fspec,'HeaderLines',7);
fclose(fid);
% time is as yymmddHHMMSS
temp.base = J{1};
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
J{1} = datenum(2000+temp.yr,temp.mnth,temp.dy,temp.hr,temp.min,temp.sec);
J{9} = J{9}/10; % rescales to match other files
G = cat(2,G,[J{1},J{9},J{8}]'); % [t P T]

% empirical trimming
G = G(:,2874:end);

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
writematrix([GG(1,:);GG(2,:)*10000;GG(3,:)]',[wrtdir 'U1519_1Hz'])

% decimation loop
tg=[];
Tg=[];
g=[];
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
    g=cat(2,g,segG(1,:));
    Tg=cat(2,Tg,segG(2,:));
    i1=i2;
    d2=floor(G(1,i2))+1;
end

% scale to hPa, tidal filter
g=g*100; % [hPa]
gf=Z_godin(g);
Tgf=Z_godin(Tg);

% remove NaNs from tidal filter
tg(isnan(gf))=[];
Tg(isnan(gf))=[];
Tgf(isnan(gf))=[];
g(isnan(gf))=[];
gf(isnan(gf))=[];

% interpolate onto monotonic time basis on the hour
tgf = tg(1)+datenum(0,0,0,0,30,0):1/24:tg(end)-datenum(0,0,0,0,30,0);
gf = interp1(tg,gf,tgf);
Tgf = interp1(tg,Tgf,tgf);

clearvars('t_str','tG','G','GG','tH','H','temp','J')

%-----PLOTTING-----%
figure(3); clf; hold on
% BPR
af_plot=detrend(af);
plot(taf,af_plot,'linewidth',1)
text(taf(end)+10,af_plot(end-100),{staname{1}},'fontsize',14)
% U1518
ef_plot=detrend(ef);
plot(tef,ef_plot+15,'linewidth',1)
text(tef(end)+10,ef_plot(end-100)+15,{staname{2}},'fontsize',14)
% U1519
gf_plot=detrend(gf);
plot(tgf,gf_plot+30,'linewidth',1)
text(tgf(end)+10,gf_plot(end-100)+30,{staname{3}},'fontsize',14)
xlim([datenum(2018,01,01) datenum(2024,01,01)])
datetick('x',10,'keeplimits')
set(gca,'fontsize',14)
ylabel('\DeltaP (cm)')
ylim([-10 40])
box on; grid on
title('CORK wellheads (2018-2023)')

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print(figdir,'-dpng','-r300')

% combine and save data
t={ta,te,tg};
p={a,e,g};
T={Ta,Te,Tg};
tf={taf,tef,tgf};
pf={af,ef,gf};
Tf={Taf,Tef,Tgf};
save(svdir,'t','p','T','tf','pf','Tf','staname','stalat','stalon','stadepth')

% write to text files
for i=1:length(t)
    % decimated
    writematrix([t{i}',p{i}',T{i}'],[wrtdir staname{i} '_1hr'])
    % decimated, detided
    writematrix([tf{i}',pf{i}',Tf{i}'],[wrtdir staname{i} '_1hr_detided'])
end
