%% HOBITSS VIII

clear
load('../processed_data/HOBITSS_VIII.mat');

% sort by latitude
[lat,id]=sort(stalat);
name=staname(id);
lon=stalon(id);
depth=stadepth(id);
tin=tf(id);
pin=pf(id);

%----- CEOF + exponential
% unify the time basis
t_min=max(cellfun(@(v)v(1),tin));
t_max=min(cellfun(@(v)v(end),tin));

for kk=1:length(pin)
    cond=(tin{kk}>=t_min+0.01) & (tin{kk}<=t_max+0.01);
    ttemp{kk}=tin{kk}(cond);
    ptemp{kk}=pin{kk}(cond); ptemp{kk}=detrend(ptemp{kk});
    dP(:,kk)=ptemp{kk}/std(ptemp{kk});
end
% calculate CEOFs from combined dataset
[~,loadings,pcs,~]=ceof(dP,1,true); % only need 1st CEOF
pred1=real(pcs(:,1) * conj(loadings(:,1)).');
% apply best CEOF "seasonal" correction
scaled_pred1=mat2cell(pred1(:,1:length(pin)).*cellfun(@std,ptemp),length(dP),ones(1,length(pin)));

figure(99); clf
subplot(121); hold on
for i=1:7
plot(ttemp{i},ptemp{i}-scaled_pred1{i}'+2*(i-1),'linewidth',1)
text(ttemp{i}(end)+10,mean(ptemp{i}-scaled_pred1{i}'+2*(i-1)),name{i})
end
ylim([-2 14])
datetick('x',6)
box on; grid on
set(gca,'fontsize',14)
subplot(122); hold on
for i=8:13
plot(ttemp{i},ptemp{i}-scaled_pred1{i}'+2*(i-8),'linewidth',1)
text(ttemp{i}(end)+10,mean(ptemp{i}-scaled_pred1{i}'+2*(i-8)),name{i})
end
ylim([-2 14])
datetick('x',6)
box on; grid on
set(gca,'fontsize',14)

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print('../figures/exploratory/HOBITSS_VIII/differences/ceof/test','-dpng','-r300')

%% HOBITSS VII

clear
load('../processed_data/HOBITSS_VII.mat');

% sort by latitude
[lat,id]=sort(stalat);
name=staname(id);
lon=stalon(id);
depth=stadepth(id);
tin=tf(id);
pin=pf(id);

%----- CEOF + exponential
% unify the time basis
t_min=max(cellfun(@(v)v(1),tin));
t_max=min(cellfun(@(v)v(end),tin));

for kk=1:length(pin)
    cond=(tin{kk}>=t_min+0.01) & (tin{kk}<=t_max+0.01);
    ttemp{kk}=tin{kk}(cond);
    ptemp{kk}=pin{kk}(cond); ptemp{kk}=detrend(ptemp{kk});
    dP(:,kk)=ptemp{kk}/std(ptemp{kk});
end
% calculate CEOFs from combined dataset
[~,loadings,pcs,~]=ceof(dP,1,true); % only need 1st CEOF
pred1=real(pcs(:,1) * conj(loadings(:,1)).');
% apply best CEOF "seasonal" correction
scaled_pred1=mat2cell(pred1(:,1:length(pin)).*cellfun(@std,ptemp),length(dP),ones(1,length(pin)));

figure(99); clf
subplot(121); hold on
for i=1:6
plot(ttemp{i},ptemp{i}-scaled_pred1{i}'+2*(i-1),'linewidth',1)
text(ttemp{i}(end)+10,mean(ptemp{i}-scaled_pred1{i}'+2*(i-1)),name{i})
end
ylim([-2 14])
datetick('x',6)
box on; grid on
set(gca,'fontsize',14)
subplot(122); hold on
for i=7:12
plot(ttemp{i},ptemp{i}-scaled_pred1{i}'+2*(i-7),'linewidth',1)
text(ttemp{i}(end)+10,mean(ptemp{i}-scaled_pred1{i}'+2*(i-7)),name{i})
end
ylim([-2 14])
datetick('x',6)
box on; grid on
set(gca,'fontsize',14)

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print('../figures/exploratory/HOBITSS_VII/differences/ceof/test','-dpng','-r300')