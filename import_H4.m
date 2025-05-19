% import_H4.m
%
% Read HOBITSS IV (2017-2018) data from text/csv files
%

clear; close all

%% HOBITSS IV (2017-2018)

topdir='/Volumes/Gorgoroth/apg_database/original/2017-2018_HOBITSS-IV/';
wrtdir='/Volumes/Gorgoroth/apg_database/processed/2017-2018_HOBITSS-IV/';
figdir='../figures/exploratory/HOBITSS_IV/hobitss4_stack';
svdir='../processed_data/HOBITSS_IV';

% station info
staname={'KU17-1','KU17-2','KU17-3','KU17-4','KU17-5','TX17-2','TX17-4'};
stalat=[-38.9117,-38.8501,-38.8932,-38.7108,-38.7251,-39.2969,-39.2327];
stalon=[178.9849,178.8724,178.7548,178.6614,178.8943,178.6532,178.4339];
stadepth=[3479,2121,1352,1038,2468,2349,1625];

%---KU17-1---%
% read in data
fid = fopen([topdir 'KU_OBP/KU17-1_108727.dat'],'r');
fspec = '%s %s %f %f'; % [date time P T]
A = textscan(fid,fspec);
fclose(fid);
t_str=cat(2,cat(1,A{1}{:}),repmat(' ',length(A{1}),1),cat(1,A{2}{:}));
tA = datenum(t_str,'yyyy/mm/dd HH:MM:SS');
A = [tA';A{3}';A{4}']; % [t P T]

% empirical trimming
A = A(:,33364:end);

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
writematrix([AA(1,:);AA(2,:)*100;AA(3,:)]',[wrtdir 'KU17-1_1Hz']) % [t P T]

% decimation loop down from 1 Hz
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

% this yields ~293 days of data -- did instrument fail?
clearvars('t_str','tA','A','AA')

%---KU17-2---%
% read in data
fid = fopen([topdir 'KU_OBP/KU17-2_120398.dat'],'r');
fspec = '%s %s %f %f'; % [date time P T]
B = textscan(fid,fspec);
fclose(fid);
t_str=cat(2,cat(1,B{1}{:}),repmat(' ',length(B{1}),1),cat(1,B{2}{:}));
tB = datenum(t_str,'yyyy/mm/dd HH:MM:SS');
B = [tB';B{3}';B{4}']; % [t P T]

% empirical trimming
B = B(:,233745:40625239);

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
writematrix([BB(1,:);BB(2,:)*100;BB(3,:)]',[wrtdir 'KU17-2_1Hz']) % [t P T]

% decimation loop down from 1 Hz
tb=[];
b=[];
Tb=[];
i1 = 1;
while i1<length(B)
    [~,i2] = min(abs(B(1,:)-(floor(B(1,i1))+1)));
    [segt,segB,~,~] = downsample_uneven(B(1,i1:i2-1),B(2:3,i1:i2-1),1/24);
    tb=cat(2,tb,segt);
    b=cat(2,b,segB(1,:));
    Tb=cat(2,Tb,segB(2,:));
    i1=i2;
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

% this yields ~468 days of data, which matches cruise reports
clearvars('t_str','tB','B','BB')

%---KU17-3---%
% read in data
fid = fopen([topdir 'KU_OBP/KU17-3_131477.dat'],'r');
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
C = C(:,19143:20224947);

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
writematrix([CC(3,:);CC(2,:)*100;CC(1,:)]',[wrtdir 'KU17-3_1Hz']) % [t P T]

% decimation loop down from 0.5 Hz
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

% this yields ~468 days of data, which matches cruise reports
clearvars('temp','C','CC')

%---KU17-4---%
% read in data
fid = fopen([topdir 'KU_OBP/KU17-4_131479.dat'],'r');
fspec = '%f %f %f'; % [P T t]
sizeD = [3 Inf];
D = fscanf(fid,fspec,sizeD);
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
D = D(:,84032:20325014);

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
writematrix([DD(3,:);DD(2,:)*100;DD(1,:)]',[wrtdir 'KU17-4_1Hz']) % [t P T]

% decimation loop down from 0.5 Hz
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

% this yields ~469 days of data, which matches cruise reports
clearvars('temp','D','DD')

%---KU17-5---%
% read in data
fid = fopen([topdir 'KU_OBP/KU17-5_126652.dat'],'r');
fspec = '%f %f %f'; % [P T t]
sizeE = [3 Inf];
E = fscanf(fid,fspec,sizeE);
fclose(fid);

% time is as yymmddHHMMSS
temp.base = E(3,:);
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
E(3,:) = datenum(2000+temp.yr,temp.mnth,temp.dy,temp.hr,temp.min,temp.sec);

% empirical trimming
E = E(:,71705:20451046);

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
writematrix([EE(3,:);EE(2,:)*100;EE(1,:)]',[wrtdir 'KU17-5_1Hz']) % [t P T]

% decimation loop down from 0.5 Hz
te=[];
e=[];
Te=[];
i1 = 1;
d2 = floor(E(3,1))+1;
while i1<length(E)
    i2 = find(E(3,:)>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segE,~,~] = downsample_uneven(E(3,i1:i2-1),E(1:2,i1:i2-1),1/24);
    if length(segt)>24
        keyboard
    end
    te=cat(2,te,segt);
    e=cat(2,e,segE(1,:));
    Te=cat(2,Te,segE(2,:));
    i1=i2;
    d2=floor(E(3,i2))+1;
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

% this yields ~472 days of data, which matches cruise reports
clearvars('temp','E','EE')

%---TX17-2---%
% already as matfile
temp = load([topdir 'data_decimated/TX17-2.mat']);
F = [temp.TXBPR2_P';temp.TXBPR2_Temp';temp.TXBPR2_T']; % [P T t]

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
writematrix(FF',[wrtdir 'TX17-2_1Hz']) % [t P T]

% decimation loop down from 100 Hz
tf=[];
f=[];
Tf=[];
i1 = 1;
while i1<length(F)
    [~,i2] = min(abs(F(3,:)-(F(3,i1)+1)));
    [segt,segF,~,~] = downsample_uneven(F(3,i1:i2-1),F(1:2,i1:i2-1),1/24);
    tf=cat(2,tf,segt);
    f=cat(2,f,segF(1,:));
    Tf=cat(2,Tf,segF(2,:));
    i1=i2;
end

% convert to hPa
f = f/100;

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
tff = tff+datenum(2017,0,0,0,0,0); % [datenum]

% data ends ~1.5 months earlier than expected from cruise reports
clearvars('temp','F','FF')

%---TX17-4---%
% already as matfile
temp = load([topdir 'data_decimated/TX17-4.mat']);
G = [temp.TXBPR4_P';temp.TXBPR4_Temp';temp.TXBPR4_T']; % [P T t]

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
writematrix(GG',[wrtdir 'TX17-4_1Hz']) % [t P T]

% decimation loop down from 100 Hz
tg=[];
g=[];
Tg=[];
i1 = 1;
while i1<length(G)
    [~,i2] = min(abs(G(3,:)-(G(3,i1)+1)));
    [segt,segG,~,~] = downsample_uneven(G(3,i1:i2-1),G(1:2,i1:i2-1),1/24);
    tg=cat(2,tg,segt);
    g=cat(2,g,segG(1,:));
    Tg=cat(2,Tg,segG(2,:));
    i1=i2;
end

% convert to hPa
g = g/100;

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
tgf = tgf+datenum(2017,0,0,0,0,0);

% data ends ~1.5 months earlier than expected from cruise reports
clearvars('temp','G','GG')

%-----PLOTTING-----%
figure(2); clf; hold on
af_plot = detrend(af);
plot(taf,af_plot,'linewidth',1)
text(taf(end)+10,af_plot(end-100),{staname{1};[num2str(stadepth(1)) ' m']},'fontsize',14)
ef_plot = detrend(ef);
plot(tef,ef_plot+20,'linewidth',1)
text(tef(end)+10,ef_plot(end-100)+20,{staname{5};[num2str(stadepth(5)) ' m']},'fontsize',14)
ff_plot = detrend(ff);
plot(tff,ff_plot+40,'linewidth',1)
text(tff(end)+10,ff_plot(end-100)+40,{staname{6};[num2str(stadepth(6)) ' m']},'fontsize',14)
bf_plot = detrend(bf);
plot(tbf,bf_plot+60,'linewidth',1)
text(tbf(end)+10,bf_plot(end-100)+60,{staname{2};[num2str(stadepth(2)) ' m']},'fontsize',14)
gf_plot = detrend(gf);
plot(tgf,gf_plot+80,'linewidth',1)
text(tgf(end)+10,gf_plot(end-100)+80,{staname{7};[num2str(stadepth(7)) ' m']},'fontsize',14)
cf_plot = detrend(cf);
plot(tcf,cf_plot+100,'linewidth',1)
text(tcf(end)+10,cf_plot(end-100)+100,{staname{3};[num2str(stadepth(3)) ' m']},'fontsize',14)
df_plot = detrend(df);
plot(tdf,df_plot+120,'linewidth',1)
text(tdf(end)+10,df_plot(end-100)+120,{staname{4};[num2str(stadepth(4)) ' m']},'fontsize',14)
xlim([datenum(2017,06,01) datenum(2019,01,01)])
datetick('x','keeplimits')
set(gca,'fontsize',14)
ylabel('\DeltaP (hPa)')
title('HOBITSS 2017-2018')
box on; grid on

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print(figdir,'-dpng','-r300')

% combine and save data
t={ta,tb,tc,td,te,tf,tg};
p={a,b,c,d,e,f,g};
T={Ta,Tb,Tc,Td,Te,Tf,Tg};
tf={taf,tbf,tcf,tdf,tef,tff,tgf};
pf={af,bf,cf,df,ef,ff,gf};
Tf={Taf,Tbf,Tcf,Tdf,Tef,Tff,Tgf};
save(svdir,'t','p','T','tf','pf','Tf','staname','stalat','stalon','stadepth')

% write to text files
for i=1:length(t)
    % decimated
    writematrix([t{i}',p{i}',T{i}'],[wrtdir staname{i} '_1hr'])
    % decimated, detided
    writematrix([tf{i}',pf{i}',Tf{i}'],[wrtdir staname{i} '_1hr_detided'])
end