% import_G2.m
%
% Read GONDOR II (2023-2024) APG data from text/csv files
%

clear; close all

%% GONDOR II (2023-2024)

% parent directory
topdir='/Volumes/Gorgoroth/apg_database/original/2023-2024_GONDOR-II/';
wrtdir='/Volumes/Gorgoroth/apg_database/processed/2023-2024_GONDOR-II/';
figdir='../figures/exploratory/GONDOR_II/gondor2_stack';
svdir='../processed_data/GONDOR_II';

% station info
mdat=readtable([topdir 'GNS-BPR-surveyed_locations.csv']);
staname=table2cell(mdat(:,1))';
stalat=table2array(mdat(:,3))';
stalon=table2array(mdat(:,2))';
stadepth=table2array(mdat(:,4))';

% remove PH from the list (too short)
ibad=find(strcmp(staname,'GNS23-PH'));
staname(ibad)=[];
stalat(ibad)=[];
stalon(ibad)=[];
stadepth(ibad)=[];

% remove PN from the list (unusable data quality)
ibad=find(strcmp(staname,'GNS23-PN'));
staname(ibad)=[];
stalat(ibad)=[];
stalon(ibad)=[];
stadepth(ibad)=[];

%---GNS22-PK2---%
% read in data
fid = [topdir 'final/GNS22-PK2.csv'];
K = readtable(fid); % [t P T]
K = [datenum(table2array(K(:,1))),table2array(K(:,2:3))]';

% empirical trimming
K = K(:,3249450:51061600);

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
writematrix(KK',[wrtdir 'GNS22-PK2_1Hz']) % [t P T]

% sampled at 1 Hz, so just need to ensure even sampling
tks=K(1,1):1/86400:K(1,end);
ks=interp1(K(1,:),K(2,:),tks);
Tks=interp1(K(1,:),K(3,:),tks);

% clean up anomalous spike
ks(1672675:1675091)=linspace(ks(1672675),ks(1675091),1675091-1672675+1);
Tks(1672675:1675091)=linspace(Tks(1672675),Tks(1675091),1675091-1672675+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% there are a series of irregular behaviors that cause an apparent spike
% in the filtered data. it recovers, but should be noted. approximately
% during the month of November
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%---GNS22-PM---%
% read in data
fid = [topdir 'final/GNS22-PM.csv'];
M = readtable(fid); % [t P T]
M = [datenum(table2array(M(:,1))),table2array(M(:,2:3))]';

% empirical trimming
M = M(:,8366:61519215);

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
writematrix(MM',[wrtdir 'GNS22-PM_1Hz']) % [t P T]

% sampled at 1 Hz, so just need to ensure even sampling
tms=M(1,1):1/86400:M(1,end);
ms=interp1(M(1,:),M(2,:),tms);
Tms=interp1(M(1,:),M(3,:),tms);

% decimation loop to 1 sample/hour
tm=[];
m=[];
Tm=[];
i1 = 1;
d2 = floor(tms(1))+1;
while i1<length(tms)
    i2 = find(tms>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segM,~,~] = downsample_uneven(tms(i1:i2-1),[ms(i1:i2-1);Tms(i1:i2-1)],1/24);
    if length(segt)>24
        keyboard
    end
    tm=cat(2,tm,segt);
    m=cat(2,m,segM(1,:));
    Tm=cat(2,Tm,segM(2,:));
    i1=i2;
    d2=floor(tms(i2))+1;
end

% tidal filter, convert to hPa
ms=ms/100;
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
mf = interp1(tm,mf,tmf);
Tmf = interp1(tm,Tmf,tmf);

clearvars('M','MM')

%---GNS22-PP---%
% read in data
fid = [topdir 'final/GNS22-PP.csv'];
P = readtable(fid); % [t P T]
P = [datenum(table2array(P(:,1))),table2array(P(:,2:3))]';

% empirical trimming
P = P(:,14700:55094000);

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
writematrix(PP',[wrtdir 'GNS22-PP_1Hz']) % [t P T]

% sampled at 1 Hz, so just need to ensure even sampling
tps=P(1,1):1/86400:P(1,end);
ps=interp1(P(1,:),P(2,:),tps);
Tps=interp1(P(1,:),P(3,:),tps);

% decimation loop to 1 sample/hour
tp=[];
p=[];
Tp=[];
i1 = 1;
d2 = floor(tps(1))+1;
while i1<length(tps)
    i2 = find(tps>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segP,~,~] = downsample_uneven(tps(i1:i2-1),[ps(i1:i2-1);Tps(i1:i2-1)],1/24);
    if length(segt)>24
        keyboard
    end
    tp=cat(2,tp,segt);
    p=cat(2,p,segP(1,:));
    Tp=cat(2,Tp,segP(2,:));
    i1=i2;
    d2=floor(tps(i2))+1;
end

% tidal filter, convert to hPa
ps=ps/100;
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
pf = interp1(tp,pf,tpf);
Tpf = interp1(tp,Tpf,tpf);

clearvars('P','PP')

%---GNS23-PF---%
% read in data
fid = [topdir 'final/GNS23-PF.csv'];
F = readtable(fid); % [t P T]
F = [datenum(table2array(F(:,1))),table2array(F(:,2:3))]';

% empirical trimming
F = F(:,8050:28169118);

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
writematrix(FF',[wrtdir 'GNS23-PF_1Hz']) % [t P T]

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

%---GNS23-PH---%
% read in data
fid = [topdir 'final/GNS23-PH.csv'];
% H = readtable(fid); % [t P T]
% H = [datenum(table2array(H(:,1))),table2array(H(:,2:3))]';
% 
% % empirical trimming
% % can probably get more out of this, but will require careful work and
% % probably not worth it, as max duration is only 3 months
% H = H(:,6100:5008350);
% 
% % write to text file, with pressure as Pa
% writematrix([H(1,:);H(2,:);H(3,:)]',[wrtdir 'GNS23-PH_1Hz']) % [t P T]
% 
% % sampled at 1 Hz, so just need to ensure even sampling
% ths=H(1,1):1/86400:H(1,end);
% hs=interp1(H(1,:),H(2,:),ths);
% Ths=interp1(H(1,:),H(3,:),ths);
% 
% % decimation loop to 1 sample/hour
% th=[];
% h=[];
% Th=[];
% i1 = 1;
% d2 = floor(ths(1))+1;
% while i1<length(ths)
%     i2 = find(ths>=d2,1);
%     if isempty(i2)
%         break
%     end
%     [segt,segH,~,~] = downsample_uneven(ths(i1:i2-1),[hs(i1:i2-1);Ths(i1:i2-1)],1/24);
%     if length(segt)>24
%         keyboard
%     end
%     th=cat(2,th,segt);
%     h=cat(2,h,segH(1,:));
%     Th=cat(2,Th,segH(2,:));
%     i1=i2;
%     d2=floor(ths(i2))+1;
% end
% 
% % tidal filter, convert to hPa
% hs=hs/100;
% h=h/100;
% hf=Z_godin(h);
% Thf=Z_godin(Th);
% 
% % remove NaNs from tidal filter
% th(isnan(hf))=[];
% Th(isnan(hf))=[];
% Thf(isnan(hf))=[];
% h(isnan(hf))=[];
% hf(isnan(hf))=[];
% 
% % interpolate onto monotonic time basis
% thf = th(1)+datenum(0,0,0,0,30,0):1/24:th(end)-datenum(0,0,0,0,30,0);
% hf = interp1(th,hf,thf);
% Thf = interp1(th,Thf,thf);
% 
% clearvars('H')

%---GNS23-PI---%
% read in data
fid = [topdir 'final/GNS23-PI.csv'];
H = readtable(fid); % [t P T]
H = [datenum(table2array(H(:,1))),table2array(H(:,2:3))]';

% empirical trimming
H = H(:,7050:28334500);

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
writematrix(HH',[wrtdir 'GNS23-PI_1Hz']) % [t P T]

% sampled at 1 Hz, so just need to ensure even sampling
tzs=H(1,1):1/86400:H(1,end);
zs=interp1(H(1,:),H(2,:),tzs);
Tzs=interp1(H(1,:),H(3,:),tzs);

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

clearvars('H','HH')

%---GNS23-PN---%
%---BAD DATA, NOT WORTH DOING---%
% read in data
fid = [topdir 'final/GNS23-PN.csv'];

%% -----PLOTTING----- %%

figure(3); clf; hold on
kf_plot = detrend(kf);
plot(tkf,kf_plot,'linewidth',1)
text(tkf(end)+10,kf_plot(end-100),{staname{1}},'fontsize',14)
mf_plot = detrend(mf);
plot(tmf,mf_plot+10,'linewidth',1)
text(tmf(end)+10,mf_plot(end-100)+10,{staname{2}},'fontsize',14)
pf_plot = detrend(pf);
plot(tpf,pf_plot+20,'linewidth',1)
text(tpf(end)+10,pf_plot(end-100)+20,{staname{3}},'fontsize',14)
ff_plot = detrend(ff);
plot(tff,ff_plot+30,'linewidth',1)
text(tff(end)+10,ff_plot(end-100)+30,{staname{4}},'fontsize',14)
zf_plot = detrend(zf);
plot(tzf,zf_plot+40,'linewidth',1)
text(tzf(end)+10,zf_plot(end-100)+40,{staname{5}},'fontsize',14)
xlim([datenum(2022,10,01) datenum(2025,01,01)])
datetick('x','keeplimits')
set(gca,'fontsize',14)
ylabel('\DeltaP (hPa)')
box on; grid on
title('GONDOR II 2023-2024')

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print(figdir,'-dpng','-r300')

% combine and save data
t={tk,tm,tp,tf,tz};
p={k,m,p,f,z};
T={Tk,Tm,Tp,Tf,Tz};
tf={tkf,tmf,tpf,tff,tzf};
pf={kf,mf,pf,ff,zf};
Tf={Tkf,Tmf,Tpf,Tff,Tzf};
save(svdir,'t','p','T','tf','pf','Tf','staname','stalat','stalon','stadepth')

% write to text files
for i=1:length(t)
    % decimated
    writematrix([t{i}',p{i}',T{i}'],[wrtdir staname{i} '_1hr'])
    % decimated, detided
    writematrix([tf{i}',pf{i}',Tf{i}'],[wrtdir staname{i} '_1hr_detided'])
end

% combine and save 1 Hz data
ts={tks,tms,tps,tfs,tzs};
ps={ks,ms,ps,fs,zs};
Ts={Tks,Tms,Tps,Tfs,Tzs};
save([svdir '_1Hz'],'ts','ps','Ts','staname','stalat','stalon','stadepth')