% AGU_abstract_figure.m
%
% generate a demonstrative figure for my AGU abstract
%

clear; close all

%--- load pressures
load('../processed_data/HOBITSS_VIII_dedrifted_cork.mat');
% unify the time basis
t_min=min(cellfun(@(v)v(1),tf));
t_max=max(cellfun(@(v)v(end),tf));

%--- load and prep GPS
gps=load('../processed_data/GPS.mat','MAHI','PARI','MAKO');
% MAHI
igps=gps.MAHI(1,:)>=t_min & gps.MAHI(1,:)<=t_max;
mahi=gps.MAHI(:,igps);
mg=polyfit(1:250,mahi(2,1:250),1);
mahi(2,:)=mahi(2,:)-(1:length(mahi))*mg(1); mahi(2,:)=mahi(2,:)-mean(mahi(2,:));
% PARI
igps=gps.PARI(1,:)>=t_min & gps.PARI(1,:)<=t_max;
pari=gps.PARI(:,igps);
mg=polyfit(1:250,pari(2,1:250),1);
pari(2,:)=pari(2,:)-(1:length(pari))*mg(1); pari(2,:)=pari(2,:)-mean(pari(2,:));
% MAKO
igps=gps.MAKO(1,:)>=t_min & gps.MAKO(1,:)<=t_max;
mako=gps.MAKO(:,igps);
mg=polyfit(1:250,mako(2,1:250),1);
mako(2,:)=mako(2,:)-(1:length(mako))*mg(1); mako(2,:)=mako(2,:)-mean(mako(2,:));

%--- prep CORK
cork=load('../processed_data/CORK_formation.mat');
% U1518
icork=cork.tf{1}>=t_min & cork.tf{1}<=t_max;
u1518(1,:)=cork.tf{1}(icork);
u1518(2,:)=cork.p2{1}(icork);
mc=polyfit(1:6750,u1518(2,1:6750),1);
u1518(2,:)=u1518(2,:)-(1:length(u1518))*mc(1); u1518(2,:)=u1518(2,:)-mean(u1518(2,:));
% U1519
icork=cork.tf{2}>=t_min & cork.tf{2}<=t_max;
u1519(1,:)=cork.tf{2}(icork);
u1519(2,:)=cork.p2{2}(icork); 
mc=polyfit(1:6750,u1519(2,1:6750),1);
u1519(2,:)=u1519(2,:)-(1:length(u1519))*mc(1); u1519(2,:)=u1519(2,:)-mean(u1519(2,:));

%--- plot time series
figure(86); clf;
subplot(121); hold on
plot(mahi(1,:),mahi(2,:)/10,'k','linewidth',1)
text(mahi(1,end)+10,mahi(2,end)/10,'MAHI','fontsize',14)
plot(pari(1,:),pari(2,:)/10+2,'k','linewidth',1)
text(pari(1,end)+10,pari(2,end)/10+2,'PARI','fontsize',14)
plot(mako(1,:),mako(2,:)/10+4,'k','linewidth',1)
text(mako(1,end)+10,mako(2,end)/10+4,'MAKO','fontsize',14)
plot(u1518(1,:),u1518(2,:)/3+8,'k','linewidth',1)
text(u1518(1,end)+10,u1518(2,end)/3+8,'U1518','fontsize',14)
plot(u1519(1,:),u1519(2,:)/3+10,'k','linewidth',1)
text(u1519(1,end)+10,u1519(2,end)/3+10,'U1519','fontsize',14)
xlim([datenum(2021,10,01) datenum(2023,01,01)])
datetick('x',6,'keeplimits')
set(gca,'fontsize',16)
ylabel('GPS East (arbitrary scale) Formation Pressure (arbitrary scale)')
box on; grid on;
subplot(122); hold on
plot(tfc{1},ccor{1},'r','linewidth',1)
text(tfc{1}(end)+10,0,{'GNS21-PA';'2265 m'},'fontsize',14)
plot(tfc{4},ccor{4}+2.5,'b','linewidth',1)
text(tfc{4}(end)+10,2.5,{'GNS21-PD';'1265 m'},'fontsize',14)
plot(tfc{8},ccor{8}+6,'r','linewidth',1)
text(tfc{8}(end)+10,6,{'GNS21-PJ';'2461 m'},'fontsize',14)
plot(tfc{7},ccor{7}+8.5,'b','linewidth',1)
text(tfc{7}(end)+10,8.5,{'GNS21-PI';'992 m'},'fontsize',14)
plot(tfc{13},ccor{13}+12,'r','linewidth',1)
text(tfc{13}(end)+13,12,{'TU21-PE';'1904 m'},'fontsize',14)
plot(tfc{12},ccor{12}+14,'b','linewidth',1)
text(tfc{12}(end)+10,14,{'TU21-PD';'1194 m'},'fontsize',14)
xlim([datenum(2021,10,01) datenum(2023,01,01)])
datetick('x',6,'keeplimits')
set(gca,'fontsize',16)
ylabel('Pressure (cm)')
box on; grid on;
fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 17 11];
print('../figures/exploratory/HOBITSS_VIII/differences/AGU-results','-dpng','-r300')
print('../figures/exploratory/HOBITSS_VIII/differences/AGU-results','-depsc','-vector')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2022
%--- load pressures
load('../processed_data/HOBITSS_VIII_dedrifted_cork.mat');
% unify the time basis
t_min=min(cellfun(@(v)v(1),tf));
t_max=max(cellfun(@(v)v(end),tf));

%--- load and prep GPS
gps=load('../processed_data/GPS.mat','MAHI','PARI','MAKO');
% MAHI
igps=gps.MAHI(1,:)>=t_min & gps.MAHI(1,:)<=t_max;
mahi=gps.MAHI(:,igps);
mg=polyfit(1:250,mahi(2,1:250),1);
mahi(2,:)=mahi(2,:)-(1:length(mahi))*mg(1); mahi(2,:)=mahi(2,:)-mean(mahi(2,:));
disp(['2022, MAHI, ' num2str(round(mean(mahi(2,end-50:end))-mean(mahi(2,1:51))))])
% PARI
igps=gps.PARI(1,:)>=t_min & gps.PARI(1,:)<=t_max;
pari=gps.PARI(:,igps);
mg=polyfit(1:250,pari(2,1:250),1);
pari(2,:)=pari(2,:)-(1:length(pari))*mg(1); pari(2,:)=pari(2,:)-mean(pari(2,:));
disp(['2022, PARI, ' num2str(round(mean(pari(2,end-50:end))-mean(pari(2,1:51))))])
% MAKO
igps=gps.MAKO(1,:)>=t_min & gps.MAKO(1,:)<=t_max;
mako=gps.MAKO(:,igps);
mg=polyfit(1:250,mako(2,1:250),1);
mako(2,:)=mako(2,:)-(1:length(mako))*mg(1); mako(2,:)=mako(2,:)-mean(mako(2,:));
disp(['2022, MAKO, ' num2str(round(mean(mako(2,end-50:end))-mean(mako(2,1:51))))])

figure(82); clf; hold on
h1=plot(mahi(1,:),mahi(2,:),'k','linewidth',1);
text(mahi(1,end)+10,mahi(2,end),'MAHI','fontsize',14)
plot(pari(1,:),pari(2,:)+20,'k','linewidth',1)
text(pari(1,end)+10,pari(2,end)+20,'PARI','fontsize',14)
plot(mako(1,:),mako(2,:)+40,'k','linewidth',1)
text(mako(1,end)+10,mako(2,end)+40,'MAKO','fontsize',14)
xlim([datenum(2021,10,01) datenum(2023,01,01)])
datetick('x',6,'keeplimits')
set(gca,'fontsize',16)
ylabel('GPS East (mm)')
box on; grid on;

%% 2019
t_min=t_min-3*365;
t_max=t_max-3*365;
%--- load and prep GPS
gps=load('../processed_data/GPS.mat','MAHI','PARI','MAKO');
% MAHI
igps=gps.MAHI(1,:)>=t_min & gps.MAHI(1,:)<=t_max;
mahi=gps.MAHI(:,igps);
mg=polyfit(1:125,mahi(2,1:125),1);
mahi(2,:)=mahi(2,:)-(1:length(mahi))*mg(1); mahi(2,:)=mahi(2,:)-mean(mahi(2,:));
disp(['2019, MAHI, ' num2str(round(mean(mahi(2,end-50:end))-mean(mahi(2,1:51))))])
% PARI
igps=gps.PARI(1,:)>=t_min & gps.PARI(1,:)<=t_max;
pari=gps.PARI(:,igps);
mg=polyfit(1:125,pari(2,1:125),1);
pari(2,:)=pari(2,:)-(1:length(pari))*mg(1); pari(2,:)=pari(2,:)-mean(pari(2,:));
disp(['2019, PARI, ' num2str(round(mean(pari(2,end-50:end))-mean(pari(2,1:51))))])
% MAKO
igps=gps.MAKO(1,:)>=t_min & gps.MAKO(1,:)<=t_max;
mako=gps.MAKO(:,igps);
mg=polyfit(1:125,mako(2,1:125),1);
mako(2,:)=mako(2,:)-(1:length(mako))*mg(1); mako(2,:)=mako(2,:)-mean(mako(2,:));
disp(['2019, MAKO, ' num2str(round(mean(mako(2,end-50:end))-mean(mako(2,1:51))))])

h2=plot(mahi(1,:)+3*365,mahi(2,:)+80,'b','linewidth',1);
text(mahi(1,end)+3*365+10,mahi(2,end)+80,'MAHI','fontsize',14)
plot(pari(1,:)+3*365,pari(2,:)+100,'b','linewidth',1)
text(pari(1,end)+3*365+10,pari(2,end)+100,'PARI','fontsize',14)
plot(mako(1,:)+3*365,mako(2,:)+120,'b','linewidth',1)
text(mako(1,end)+3*365+10,mako(2,end)+120,'MAKO','fontsize',14)

%% 2014
t_min=t_min-4.5*365;
t_max=t_max-4.5*365;
%--- load and prep GPS
gps=load('../processed_data/GPS.mat','MAHI','PARI','MAKO');
% MAHI
igps=gps.MAHI(1,:)>=t_min & gps.MAHI(1,:)<=t_max;
mahi=gps.MAHI(:,igps);
mg=polyfit(1:125,mahi(2,1:125),1);
mahi(2,:)=mahi(2,:)-(1:length(mahi))*mg(1); mahi(2,:)=mahi(2,:)-mean(mahi(2,:));
disp(['2014, MAHI, ' num2str(round(mean(mahi(2,end-50:end))-mean(mahi(2,1:51))))])
% PARI
igps=gps.PARI(1,:)>=t_min & gps.PARI(1,:)<=t_max;
pari=gps.PARI(:,igps);
mg=polyfit(1:125,pari(2,1:125),1);
pari(2,:)=pari(2,:)-(1:length(pari))*mg(1); pari(2,:)=pari(2,:)-mean(pari(2,:));
disp(['2014, PARI, ' num2str(round(mean(pari(2,end-50:end))-mean(pari(2,1:51))))])
% MAKO
igps=gps.MAKO(1,:)>=t_min & gps.MAKO(1,:)<=t_max;
mako=gps.MAKO(:,igps);
mg=polyfit(1:125,mako(2,1:125),1);
mako(2,:)=mako(2,:)-(1:length(mako))*mg(1); mako(2,:)=mako(2,:)-mean(mako(2,:));
disp(['2014, MAKO, ' num2str(round(mean(mako(2,end-50:end))-mean(mako(2,1:51))))])

h3=plot(mahi(1,:)+7.5*365,mahi(2,:)+160,'r','linewidth',1);
text(mahi(1,end)+7.5*365+10,mahi(2,end)+160,'MAHI','fontsize',14)
plot(pari(1,:)+7.5*365,pari(2,:)+180,'r','linewidth',1)
text(pari(1,end)+7.5*365+10,pari(2,end)+180,'PARI','fontsize',14)
plot(mako(1,:)+7.5*365,mako(2,:)+200,'r','linewidth',1)
text(mako(1,end)+7.5*365+10,mako(2,end)+200,'MAKO','fontsize',14)
legend([h3,h2,h1],'2014 SSE','2019 SSE','2022 SSE','location','northwest')

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print('../figures/exploratory/HOBITSS_VIII/differences/AGU-gps-comp','-dpng','-r300')
print('../figures/exploratory/HOBITSS_VIII/differences/AGU-gps-comp','-depsc','-vector')
