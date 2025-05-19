% import_H2.m
%
% Read HOBITSS II (2015-2016) data from text/csv files
%

clear; close all

%% HOBITSS II (2015-2016)

topdir='/Volumes/Gorgoroth/apg_database/original/2015-2016_HOBITSS-II/';
wrtdir='/Volumes/Gorgoroth/apg_database/processed/2015-2016_HOBITSS-II/';
figdir='../figures/exploratory/HOBITSS_II/hobitss2_stack';
svdir='../processed_data/HOBITSS_II';

% station info
staname={'KU15-1','KU15-3','KU15-4','KU15-5'};
stalat=[0,0,0,0];
stalon=[0,0,0,0];
stadepth=[0,0,0,0];

%---KU15-1---%
% read in data
fid = fopen([topdir 'KU15-1.txt'],'r');
A = textscan(fid,'%f %f %f'); % [P T t]
fclose(fid);

% time is as yymmddHHMMSS
temp.base = A{3};
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
A{3} = datenum(2000+temp.yr,temp.mnth,temp.dy,temp.hr,temp.min,temp.sec);

A = [A{3},A{1},A{2}]'; % [t P T]

% empirical trimming
A = A(:,281340:16030200);

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
writematrix([AA(1,:);AA(2,:)*100;AA(3,:)]',[wrtdir 'KU15-1_1Hz']) % [t P T]

% decimation loop
ta=[];
a=[];
Ta=[];
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

clearvars('temp','A','AA')

%---KU15-2---%
% extremely messy data that only last a few weeks anyway
% % read in data
% fid = fopen([topdir 'KU15-2.txt'],'r');
% B = textscan(fid,'%f %f %f'); % [P T t]
% fclose(fid);
% 
% % time is as yymmddHHMMSS
% temp.base = B{3};
% temp.yr = floor(temp.base/10^10);
% temp.base = temp.base-temp.yr*10^10;
% temp.mnth = floor(temp.base/10^8);
% temp.base = temp.base-temp.mnth*10^8;
% temp.dy = floor(temp.base/10^6);
% temp.base = temp.base-temp.dy*10^6;
% temp.hr = floor(temp.base/10^4);
% temp.base = temp.base-temp.hr*10^4;
% temp.min = floor(temp.base/10^2);
% temp.sec = temp.base-temp.min*10^2;
% B{3} = datenum(2000+temp.yr,temp.mnth,temp.dy,temp.hr,temp.min,temp.sec);
% 
% B = [B{3},B{1},B{2}]'; % [t P T]
% 
% % empirical trimming
% B = B(:,288200:18641804);
% B(:,B(2,:)>7e5)=[];
% 
% % remove duplicate/nonmonotonic timesteps
% ilist=find(diff(B(1,:))<=0);
% if ~isempty(ilist)
%     B(:,ilist+1)=[];
% end
% 
% % interpolate to 1 Hz
% BB=B(1,1):1/86400:B(1,end);
% BB(2,:)=interp1(B(1,:),B(2,:),BB); % P
% BB(3,:)=interp1(B(1,:),B(3,:),BB(1,:)); % T
% 
% % write to text file, with pressure as Pa
% writematrix([BB(1,:);BB(2,:)*100;BB(3,:)]',[wrtdir 'TU13-4_1Hz']) % [t P T]
% 
% % decimation loop
% tb=[];
% b=[];
% Tb=[];
% i1 = 1;
% d2 = floor(B(1,1))+1;
% while i1<length(B)
%     i2 = find(B(1,:)>=d2,1);
%     if isempty(i2)
%         break
%     end
%     [segt,segB,~,~] = downsample_uneven(B(1,i1:i2-1),B(2:3,i1:i2-1),1/24);
%     if length(segt)>24
%         keyboard
%     end
%     tb=cat(2,tb,segt);
%     b=cat(2,b,segB(1,:));
%     Tb=cat(2,Tb,segB(2,:));
%     i1=i2;
%     d2=floor(B(1,i2))+1;
% end
% 
% % tidal filter
% bf=Z_godin(b);
% Tbf=Z_godin(Tb);
% 
% % remove NaNs from tidal filter
% ta(isnan(bf))=[];
% Ta(isnan(bf))=[];
% Taf(isnan(bf))=[];
% b(isnan(bf))=[];
% bf(isnan(bf))=[];
% 
% % interpolate onto monotonic time basis on the hour
% tbf = tb(1)+datenum(0,0,0,0,30,0):1/24:tb(end)-datenum(0,0,0,0,30,0);
% bf = interp1(tb,bf,tbf);
% Tbf = interp1(tb,Tbf,tbf);
% 
% clearvars('temp','B','BB')

%---KU15-3---%
% read in data
fid = fopen([topdir 'KU15-3.txt'],'r');
C = textscan(fid,'%f %f %f'); % [P T t]
fclose(fid);

% time is as yymmddHHMMSS
temp.base = C{3};
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
C{3} = datenum(2000+temp.yr,temp.mnth,temp.dy,temp.hr,temp.min,temp.sec);

C = [C{3},C{1},C{2}]'; % [t P T]

% empirical trimming
C = C(:,323500:31957600);

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
writematrix([CC(1,:);CC(2,:)*100;CC(3,:)]',[wrtdir 'KU15-3_1Hz']) % [t P T]

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

% interpolate onto monotonic time basis on the hour
tcf = tc(1)+datenum(0,0,0,0,30,0):1/24:tc(end)-datenum(0,0,0,0,30,0);
cf = interp1(tc,cf,tcf);
Tcf = interp1(tc,Tcf,tcf);

clearvars('temp','C','CC')

%---KU15-4---%
% read in data
fid = fopen([topdir 'KU15-4.txt'],'r');
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
D = D(:,156450:31957500);

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
writematrix([DD(1,:);DD(2,:)*100;DD(3,:)]',[wrtdir 'KU15-4_1Hz']) % [t P T]

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
    d=cat(2,d,segD(1,:));
    Td=cat(2,Td,segD(2,:));
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

% interpolate onto monotonic time basis on the hour
tdf = td(1)+datenum(0,0,0,0,30,0):1/24:td(end)-datenum(0,0,0,0,30,0);
df = interp1(td,df,tdf);
Tdf = interp1(td,Tdf,tdf);

clearvars('temp','D','DD')

%---KU15-5---%
% read in data
fid = fopen([topdir 'KU15-5.txt'],'r');
F = textscan(fid,'%f %f %f'); % [P T t]
fclose(fid);

% time is as yymmddHHMMSS
temp.base = F{3};
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
F{3} = datenum(2000+temp.yr,temp.mnth,temp.dy,temp.hr,temp.min,temp.sec);

F = [F{3},F{1},F{2}]'; % [t P T]

% empirical trimming
F = F(:,254300:16018650);

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
writematrix([FF(1,:);FF(2,:)*100;FF(3,:)]',[wrtdir 'KU15-5_1Hz']) % [t P T]

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
    f=cat(2,f,segF(1,:));
    Tf=cat(2,Tf,segF(2,:));
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

% interpolate onto monotonic time basis on the hour
tff = tf(1)+datenum(0,0,0,0,30,0):1/24:tf(end)-datenum(0,0,0,0,30,0);
ff = interp1(tf,ff,tff);
Tff = interp1(tf,Tff,tff);

clearvars('temp','F','FF')

%-----PLOTTING-----%
figure(1); clf; hold on
af_plot = detrend(af);
plot(taf,af_plot,'linewidth',1)
text(taf(end)+10,af_plot(end-100),{staname{1};[num2str(stadepth(1)) ' m']},'fontsize',14)
cf_plot = detrend(cf);
plot(tcf,cf_plot+10,'linewidth',1)
text(tcf(end)+10,cf_plot(end-100)+10,{staname{2};[num2str(stadepth(2)) ' m']},'fontsize',14)
df_plot = detrend(df);
plot(tdf,df_plot+20,'linewidth',1)
text(tdf(end)+10,df_plot(end-100)+20,{staname{3};[num2str(stadepth(3)) ' m']},'fontsize',14)
ff_plot = detrend(ff);
plot(tff,ff_plot+30,'linewidth',1)
text(tff(end)+10,ff_plot(end-100)+30,{staname{4};[num2str(stadepth(4)) ' m']},'fontsize',14)
xlim([datenum(2015,06,01) datenum(2016,09,01)])
datetick('x','keeplimits')
set(gca,'fontsize',14)
ylabel('\DeltaP (hPa)')
title('HOBITSS 2015-2016')
box on; grid on

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print(figdir,'-dpng','-r300')

% combine and save data
t={ta,tc,td,tf};
p={a,c,d,f};
T={Ta,Tc,Td,Tf};
tf={taf,tcf,tdf,tff};
pf={af,cf,df,ff};
Tf={Taf,Tcf,Tdf,Tff};
save(svdir,'t','p','T','tf','pf','Tf','staname','stalat','stalon','stadepth')

% write to text files
for i=1:length(t)
    % decimated
    writematrix([t{i}',p{i}',T{i}'],[wrtdir staname{i} '_1hr'])
    % decimated, detided
    writematrix([tf{i}',pf{i}',Tf{i}'],[wrtdir staname{i} '_1hr_detided'])
end