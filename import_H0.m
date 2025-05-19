% import_H0.m
%
% Read 2013-2014 data from text/csv files
%

clear; close all

%% HOBITSS 0 (2013-2014)

topdir='/Volumes/Gorgoroth/apg_database/original/2013-2014_HOBITSS-0/';
wrtdir='/Volumes/Gorgoroth/apg_database/processed/2013-2014_HOBITSS-0/';
figdir='../figures/exploratory/HOBITSS_0/hobitss0_stack';
svdir='../processed_data/HOBITSS_0';

% station info
staname={'TU13-1'};
stalat=[0]; % I don't have geographic info for these stations
stalon=[0];
stadepth=[0];

%---TU13-1---%
% read in data
fid = fopen([topdir 'TU13-1.txt'],'r');
A = textscan(fid,'%f %f %f','TreatAsEmpty','########OBS_STOP########'); % [P T t]
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
A = A(:,288200:18641804);
A(:,A(2,:)>7e5)=[];

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
writematrix([AA(1,:);AA(2,:)*100;AA(3,:)]',[wrtdir 'TU13-1_1Hz']) % [t P T]

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

%---TU13-4---%
% All I have is an empty (0 byte) file. Is this a mistake or did the
% instrument not return data?
% % read in data
% fid = fopen([topdir 'TU13-4.txt'],'r');
% B = textscan(fid,'%f %f %f','TreatAsEmpty','########OBS_STOP########'); % [P T t]
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

% %-----PLOTTING-----%
% figure(1); clf; hold on
% 
% fh=gcf;
% fh.PaperUnits='inches';
% fh.PaperPosition=[0 0 8.5 11];
% print(figdir,'-dpng','-r300')

% combine and save data
t={ta};
p={a};
T={Ta};
tf={taf};
pf={af};
Tf={Taf};
save(svdir,'t','p','T','tf','pf','Tf','staname','stalat','stalon','stadepth')

% write to text files
for i=1:length(t)
    % decimated
    writematrix([t{i}',p{i}',T{i}'],[wrtdir staname{i} '_1hr'])
    % decimated, detided
    writematrix([tf{i}',pf{i}',Tf{i}'],[wrtdir staname{i} '_1hr_detided'])
end