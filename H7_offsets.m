% H7_offsets.m

%
% Try to get to the bottom of the spurious offsets in a few of the HOBITSS
% VII pressure records
%

% load up pertinent data
% Reference location (no apparent issues)
fid = '/Users/erikfred/Downloads/GNS20-PL.csv';
L = readtable(fid); % [t P T]
L = [datenum(table2array(L(:,1))),table2array(L(:,2:3))]';
% Oddities begin in early May
fid = '/Users/erikfred/Downloads/GNS20-PI.csv';
I = readtable(fid); % [t P T]
I = [datenum(table2array(I(:,1))),table2array(I(:,2:3))]';
% Oddities begin in late May
fid = '/Users/erikfred/Downloads/GNS20-PE.csv';
E = readtable(fid); % [t P T]
E = [datenum(table2array(E(:,1))),table2array(E(:,2:3))]';

% % create filter for showcasing odd behavior
% % (can be hard to see behind huge tidal signal)
% ws=1/(L(1,2)-L(1,1)); % sample frequency in 1/d
% wn=2*(1/3)/ws; % cutoff frequency to exclude tides
% [b,a]=butter(2,wn,'low');
% 
% L_filt=L(:,10000:end-10000);
% L_filt(3,:)=filtfilt(b,a,L_filt(2,:));
% 
% I_filt=I(:,10000:end-10000);
% I_filt(3,:)=filtfilt(b,a,I_filt(2,:));
% 
% E_filt=E(:,10000:end-10000);
% E_filt(3,:)=filtfilt(b,a,E_filt(2,:));

%--- decimate and filter as I normally would
L_filt=L(1:2,10000:end-10000); L_filt(2,:)=L_filt(2,:)/100;
I_filt=I(1:2,10000:end-10000); I_filt(2,:)=I_filt(2,:)/100;
E_filt=E(1:2,10000:end-10000); E_filt(2,:)=E_filt(2,:)/100;

% decimation loops
tl=[]; ti=[]; te=[];
l=[]; ii=[]; e=[];
i1 = 1;
while i1<length(L_filt)
    [~,i2] = min(abs(L_filt(1,:)-(floor(L_filt(1,i1))+1)));
    [segt,segL,~,~] = downsample_uneven(L_filt(1,i1:i2-1),L_filt(2,i1:i2-1),1/24);
    tl=cat(2,tl,segt);
    l=cat(2,l,segL);
    i1=i2;
end
i1 = 1;
while i1<length(I_filt)
    [~,i2] = min(abs(I_filt(1,:)-(floor(I_filt(1,i1))+1)));
    [segt,segI,~,~] = downsample_uneven(I_filt(1,i1:i2-1),I_filt(2,i1:i2-1),1/24);
    ti=cat(2,ti,segt);
    ii=cat(2,ii,segI);
    i1=i2;
end
i1 = 1;
while i1<length(E_filt)
    [~,i2] = min(abs(E_filt(1,:)-(floor(E_filt(1,i1))+1)));
    [segt,segE,~,~] = downsample_uneven(E_filt(1,i1:i2-1),E_filt(2,i1:i2-1),1/24);
    te=cat(2,te,segt);
    e=cat(2,e,segE);
    i1=i2;
end

% tidal filters
lf=Z_godin(l); iif=Z_godin(ii); ef=Z_godin(e);
tl(isnan(lf))=[]; ti(isnan(iif))=[]; te(isnan(ef))=[];
lf(isnan(lf))=[]; iif(isnan(iif))=[]; ef(isnan(ef))=[];

% interpolate onto monotonic time basis
tlf = tl(1)+datenum(0,0,0,0,30,0):1/24:tl(end)-datenum(0,0,0,0,30,0);
lf = interp1(tl,lf,tlf);
tif = ti(1)+datenum(0,0,0,0,30,0):1/24:ti(end)-datenum(0,0,0,0,30,0);
iif = interp1(ti,iif,tif);
tef = te(1)+datenum(0,0,0,0,30,0):1/24:te(end)-datenum(0,0,0,0,30,0);
ef = interp1(te,ef,tef);

% generate overlapping time series and differences
[ttemp1,ia1,ib1]=intersect(round(I_filt(1,:),8),round(L_filt(1,:),8));
[ttemp2,ia2,ib2]=intersect(round(E_filt(1,:),8),round(L_filt(1,:),8));
[ttemp3,ia3,ib3]=intersect(round(tef,8),round(tlf,8));

figure(9); clf;
subplot(211); hold on
plot(L_filt(1,:),L_filt(2,:),'linewidth',2)
yyaxis right
plot(I_filt(1,:),I_filt(2,:),'linewidth',2)
xlim([datenum(2021,05,10,0,0,0) datenum(2021,05,10,4,0,0)])
datetick('x','keeplimits')
subplot(212); hold on
plot(ttemp1,I_filt(2,ia1)-L_filt(2,ib1),'k','linewidth',2)
xlim([datenum(2021,05,10,0,0,0) datenum(2021,05,10,4,0,0)])
datetick('x','keeplimits')

figure(10); clf;
subplot(311); hold on
plot(L_filt(1,:),L_filt(2,:)-nanmean(L_filt(2,:)),'linewidth',2)
ylabel('\DeltaP (cm)')
yyaxis right
plot(E_filt(1,:),E_filt(2,:)-nanmean(E_filt(2,:)),'linewidth',2)
xlim([datenum(2021,05,19,12,0,0) datenum(2021,05,21,12,0,0)])
datetick('x','keeplimits')
legend('GNS20-PL','GNS20-PE')
box on; grid on
subplot(312); hold on
plot(ttemp2,detrend(E_filt(2,ia2)-L_filt(2,ib2)),'k','linewidth',2)
xlim([datenum(2021,05,15,0,0,0) datenum(2021,05,30,0,0,0)])
h1=xline(datenum(2021,05,19,12,0,0),'--');
h2=xline(datenum(2021,05,21,12,0,0),'--');
datetick('x','keeplimits')
legend('PE-PL')
ylabel('\DeltaP (cm)')
box on; grid on
subplot(313); hold on
plot(ttemp3,detrend(ef(ia3)-lf(ib3)),'k','linewidth',2)
h3=xline(datenum(2021,05,15),'--');
h4=xline(datenum(2021,05,30),'--');
datetick('x','keeplimits')
legend('PE-PL (filtered)')
ylabel('\DeltaP (cm)')
box on; grid on

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print('../figures/exploratory/HOBITSS_VII/offsets/PE_May20','-dtiff','-r300')

delete(h1); delete(h2); delete(h3); delete(h4);
subplot(311)
xlim([datenum(2021,05,5,12,0,0) datenum(2021,05,7,12,0,0)])
datetick('x','keeplimits')
subplot(312)
xlim([datenum(2021,05,1,0,0,0) datenum(2021,05,15,0,0,0)])
datetick('x','keeplimits')
h1=xline(datenum(2021,05,5,12,0,0),'--');
h2=xline(datenum(2021,05,7,12,0,0),'--');
legend('PE-PL')
subplot(313)
h3=xline(datenum(2021,05,1),'--');
h4=xline(datenum(2021,05,15),'--');
legend('PE-PL (filtered)')

print('../figures/exploratory/HOBITSS_VII/offsets/PE_May06','-dtiff','-r300')

delete(h1); delete(h2); delete(h3); delete(h4);
subplot(311)
xlim([datenum(2021,07,19,12,0,0) datenum(2021,07,20,12,0,0)])
datetick('x','keeplimits')
subplot(312)
xlim([datenum(2021,07,18,0,0,0) datenum(2021,07,26,0,0,0)])
datetick('x','keeplimits')
h1=xline(datenum(2021,07,19,12,0,0),'--');
h2=xline(datenum(2021,07,20,12,0,0),'--');
legend('PE-PL')
subplot(313)
h3=xline(datenum(2021,07,18),'--');
h4=xline(datenum(2021,07,26),'--');
legend('PE-PL (filtered)')

print('../figures/exploratory/HOBITSS_VII/offsets/PE_July20','-dtiff','-r300')

delete(h1); delete(h2); delete(h3); delete(h4);
subplot(311)
xlim([datenum(2021,07,27,00,0,0) datenum(2021,07,28,00,0,0)])
datetick('x','keeplimits')
subplot(312)
xlim([datenum(2021,07,25,0,0,0) datenum(2021,07,31,0,0,0)])
datetick('x','keeplimits')
h1=xline(datenum(2021,07,27,0,0,0),'--');
h2=xline(datenum(2021,07,28,0,0,0),'--');
legend('PE-PL')
subplot(313)
h3=xline(datenum(2021,07,25),'--');
h4=xline(datenum(2021,07,31),'--');
legend('PE-PL (filtered)')

print('../figures/exploratory/HOBITSS_VII/offsets/PE_July27','-dtiff','-r300')

delete(h1); delete(h2); delete(h3); delete(h4);
subplot(311)
xlim([datenum(2021,08,3,00,0,0) datenum(2021,08,4,00,0,0)])
datetick('x','keeplimits')
subplot(312)
xlim([datenum(2021,08,1,0,0,0) datenum(2021,08,8,0,0,0)])
datetick('x','keeplimits')
h1=xline(datenum(2021,08,3,0,0,0),'--');
h2=xline(datenum(2021,08,4,0,0,0),'--');
legend('PE-PL')
subplot(313)
h3=xline(datenum(2021,08,1),'--');
h4=xline(datenum(2021,08,8),'--');
legend('PE-PL (filtered)')

print('../figures/exploratory/HOBITSS_VII/offsets/PE_August3','-dtiff','-r300')

delete(h1); delete(h2); delete(h3); delete(h4);
subplot(311)
xlim([datenum(2021,08,19,12,0,0) datenum(2021,08,20,12,0,0)])
datetick('x','keeplimits')
subplot(312)
xlim([datenum(2021,08,17,0,0,0) datenum(2021,08,22,0,0,0)])
datetick('x',6,'keeplimits')
h1=xline(datenum(2021,08,19,0,0,0),'--');
h2=xline(datenum(2021,08,20,0,0,0),'--');
legend('PE-PL')
subplot(313)
h3=xline(datenum(2021,08,17),'--');
h4=xline(datenum(2021,08,22),'--');
legend('PE-PL (filtered)')

print('../figures/exploratory/HOBITSS_VII/offsets/PE_August19','-dtiff','-r300')
