% method_comp_plots.m
%
% compares and plots depth matched vs. network average vs. CEOF corrections
%

clear; close all

load('../processed_data/HOBITSS_VIII_dedrifted_cork.mat')

%---------- setup and on-the-fly calculations ----------%
% % ensure arrays are shaped as columns
% for ll=1:length(pf)
%     if size(tf{ll},1)<size(tf{ll},2)
%         tf{ll}=tf{ll}';
%     end
%     if size(pf{ll},1)<size(pf{ll},2)
%         pf{ll}=pf{ll}';
%     end
% end

% approximate CEOFs (no exp/lin cleanup)
% unify the time basis
t_min=max(cellfun(@(v)v(1),tf));
t_max=min(cellfun(@(v)v(end),tf));

for kk=1:length(pf)
    cond=(tf{kk}>=t_min+0.01) & (tf{kk}<=t_max+0.01);
    tc_out{kk}=tf{kk}(cond);
    ptemp{kk}=pf{kk}(cond); ptemp{kk}=detrend(ptemp{kk});
    dP(:,kk)=ptemp{kk}/std(ptemp{kk});
end
% calculate CEOFs from combined dataset
[~,loadings,pcs,~]=ceof(dP,1,true); % only need 1st CEOF
pred1=real(pcs(:,1) * conj(loadings(:,1)).');
scaled_pred1=mat2cell(pred1.*cellfun(@std,ptemp),length(dP),ones(1,length(pf)));
ceof_temp=cellfun(@minus,ptemp,scaled_pred1,'UniformOutput',false);

for iii=1:length(ceof_temp)
    for jjj=1:length(ceof_temp)
        cor_grid(iii,jjj)=xcorr(dP(:,iii),dP(:,jjj),0,'coeff');
    end
end

% empirically-determined best depth matches
sta_in= {'GNS21-PA','GNS21-PB','GNS21-PC','GNS21-PD','GNS21-PF','GNS21-PG',...
    'GNS21-PI','GNS21-PJ','TU21-PA','TU21-PB','TU21-PC','TU21-PD','TU21-PE',...
    'BPR','U1518','U1519'};
sta_ref1={'TU21-PE','TU21-PA','TU21-PA','TU21-PD','TU21-PE','TU21-PD',...
    'TU21-PD','TU21-PA','TU21-PE','GNS21-PJ','TU21-PA','TU21-PC','GNS21-PB',...
    'TU21-PB','TU21-PB','TU21-PD'};
% (references that have little/no signal in CEOF corrected data)
sta_ref2={'TU21-PA','TU21-PA','TU21-PA','TU21-PA','TU21-PA','TU21-PA',...
    'GNS21-PB','TU21-PB','TU21-PC','TU21-PA','TU21-PA','GNS21-PB','GNS21-PB',...
    'TU21-PB','TU21-PB','TU21-PA'};
% (references that have opposite sense of signal in CEOF corrected data)
sta_ref3={'GNS21-PD','GNS21-PB','GNS21-PC','GNS21-PA','GNS21-PF','GNS21-PG',...
    'GNS21-PJ','GNS21-PI','TU21-PA','TU21-PB','TU21-PC','TU21-PE','TU21-PD',...
    'BPR','U1519','U1518'};

for m=1:length(scor)
    i_match=find(strcmp(staname{m},sta_in));
    ista_ref=find(strcmp(sta_ref1{i_match},staname));
    [t1{m},ia,ib]=intersect(round(tf{m},6),round(tf{ista_ref},6));

    p1{m}=scor{m}(ia);
    p2{m}=scor{ista_ref}(ib);

    Gp=p2{m};
    mp=inv(Gp'*Gp)*Gp'*p1{m};
    p2{m}=mp*p2{m};
    pmatch{m}=p1{m}-p2{m};
end

% calculate the network average
for m=1:length(scor)
    cond=(tf{m}>=t_min+0.01) & (tf{m}<=t_max+0.01);
    ta{m}=tf{m}(cond);
    p3{m}=scor{m}(cond);
    p_avg(:,m)=scor{m}(cond)/std(scor{m}(cond));
end
p_avg=mean(p_avg,2);

for m=1:length(scor)
    pa{m}=p_avg;

    Gp=pa{m};
    mp=inv(Gp'*Gp)*Gp'*p3{m};
    pa{m}=mp*pa{m};
    pavg{m}=p3{m}-pa{m};
end
%---------- end setup ----------%

% plot proxies over pressures
n=ceil(length(pf)/2);
figure(99); clf
subplot(121); hold on
for i=1:n
    plot(tf{i},scor{i}+10*(i-1),'k','linewidth',2)
    plot(tfc{i},scaled_pred1{i}+10*(i-1),'color',[255 198 30]/255,'linewidth',1)
    plot(t1{i},p2{i}+10*(i-1),'color',[0 154 222]/255,'linewidth',1)
    plot(ta{i},pa{i}+10*(i-1),'color',[255 31 91]/255,'linewidth',1)
    text(tf{i}(end)+10,10*(i-1),{staname{i}; [' (' sta_ref1{i} ')']},'fontsize',14)
end
ylim([-10 n*10])
datetick('x',6)
box on; grid on
set(gca,'fontsize',16)
subplot(122); hold on
for i=n+1:length(pf)
    plot(tf{i},scor{i}+10*(i-n-1),'k','linewidth',2)
    plot(tfc{i},scaled_pred1{i}+10*(i-n-1),'color',[255 198 30]/255,'linewidth',1)
    plot(t1{i},p2{i}+10*(i-n-1),'color',[0 154 222]/255,'linewidth',1)
    plot(ta{i},pa{i}+10*(i-n-1),'color',[255 31 91]/255,'linewidth',1)
    text(tf{i}(end)+10,10*(i-n-1),{staname{i}; [' (' sta_ref1{i} ')']},'fontsize',14)
end
ylim([-10 n*10])
datetick('x',6)
box on; grid on
set(gca,'fontsize',16)

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 17 11];
print('../figures/exploratory/HOBITSS_VIII/differences/proxies/AGU_stack','-dpng','-r300')
print('../figures/exploratory/HOBITSS_VIII/differences/proxies/AGU_stack','-depsc','-painters')

% bar plots of proxy effectiveness
figure(65); clf; hold on
bar(1:length(scor),cellfun(@rms,scor),'facecolor','none','edgecolor','k','linewidth',1)
hb=bar(1:length(scor),[cellfun(@rms,ccor);cellfun(@rms,pmatch);cellfun(@rms,pavg)],1);
set(hb(1),'facecolor',[255 198 30]/255);
set(hb(2),'facecolor',[0 154 222]/255);
set(hb(3),'facecolor',[255 31 91]/255);
box on; grid on
ylabel('RMS (cm)')
set(gca,'xtick',1:length(scor))
set(gca,'xticklabels',staname)
xtickangle(45)
ax=gca;
ax.YAxis.FontSize=16;
ax.XAxis.FontSize=14;

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 5.5];
print('../figures/exploratory/HOBITSS_VIII/differences/proxies/AGU_bars','-dpng','-r300')
print('../figures/exploratory/HOBITSS_VIII/differences/proxies/AGU_bars','-depsc','-painters')

% plot ceof-corrected data
n=ceil(length(pf)/2);
figure(98); clf
subplot(121); hold on
for i=1:n
    plot(tfc{i},ccor{i}+4*(i-1),'color',[255 198 30]/255,'linewidth',1)
    text(tfc{i}(end)+10,4*(i-1),staname{i},'fontsize',14)
end
ylim([-4 n*4])
datetick('x',6)
box on; grid on
xline(datenum(2022,07,01),'k','linewidth',1)
xline(datenum(2022,10,01),'k','linewidth',1)
set(gca,'fontsize',16)
subplot(122); hold on
for i=n+1:length(pf)
    plot(tfc{i},ccor{i}+4*(i-n-1),'color',[255 198 30]/255,'linewidth',1)
    text(tfc{i}(end)+10,4*(i-n-1),staname{i},'fontsize',14)
end
ylim([-4 n*4])
datetick('x',6)
box on; grid on
xline(datenum(2022,07,01),'k','linewidth',1)
xline(datenum(2022,10,01),'k','linewidth',1)
set(gca,'fontsize',16)

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 17 11];
print('../figures/exploratory/HOBITSS_VIII/differences/proxies/AGU_ceof-cor_stack','-dpng','-r300')

% plot depth-match corrected pressures
n=ceil(length(pf)/2);
figure(97); clf
subplot(121); hold on
for i=1:n
    plot(t1{i},pmatch{i}+4*(i-1),'color',[0 154 222]/255,'linewidth',1)
    text(t1{i}(end)+10,4*(i-1),{staname{i}; [' (' sta_ref1{i} ')']},'fontsize',14)
end
ylim([-4 n*4])
xlim([datenum(2021,10,01) datenum(2023,01,01)])
datetick('x',6,'keeplimits')
box on; grid on
xline(datenum(2022,07,01),'k','linewidth',1)
xline(datenum(2022,10,01),'k','linewidth',1)
set(gca,'fontsize',16)
subplot(122); hold on
for i=n+1:length(pf)
    plot(t1{i},pmatch{i}+4*(i-n-1),'color',[0 154 222]/255,'linewidth',1)
    text(t1{i}(end)+10,4*(i-n-1),{staname{i}; [' (' sta_ref1{i} ')']},'fontsize',14)
end
ylim([-4 n*4])
xlim([datenum(2021,10,01) datenum(2023,01,01)])
datetick('x',6,'keeplimits')
box on; grid on
xline(datenum(2022,07,01),'k','linewidth',1)
xline(datenum(2022,10,01),'k','linewidth',1)
set(gca,'fontsize',16)

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 17 11];
print('../figures/exploratory/HOBITSS_VIII/differences/proxies/AGU_pmatch-cor_stack_ALT','-dpng','-r300')

% plot network average corrected pressures
n=ceil(length(pf)/2);
figure(99); clf
subplot(121); hold on
for i=1:n
    plot(ta{i},pavg{i}+4*(i-1),'color',[255 31 91]/255,'linewidth',1)
    text(tfc{i}(end)+10,4*(i-1),staname{i},'fontsize',14)
end
ylim([-4 n*4])
datetick('x',6)
box on; grid on
xline(datenum(2022,07,01),'k','linewidth',1)
xline(datenum(2022,10,01),'k','linewidth',1)
set(gca,'fontsize',16)
subplot(122); hold on
for i=n+1:length(pf)
    plot(ta{i},pavg{i}+4*(i-n-1),'color',[255 31 91]/255,'linewidth',1)
    text(tfc{i}(end)+10,4*(i-n-1),staname{i},'fontsize',14)
end
ylim([-4 n*4])
datetick('x',6)
box on; grid on
xline(datenum(2022,07,01),'k','linewidth',1)
xline(datenum(2022,10,01),'k','linewidth',1)
set(gca,'fontsize',16)

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 17 11];
print('../figures/exploratory/HOBITSS_VIII/differences/proxies/AGU_pavg-cor_stack','-dpng','-r300')