function [Pout]=proxy_corrections(tin,pin,station)
%
% Code for generating proxy time series for detectability testing, and
% optionally making plots
%

% unify the time basis
t_min=max(cellfun(@(v)v(1),tin));
t_max=min(cellfun(@(v)v(end),tin));

for kk=1:length(pin)
    cond=(tin{kk}>=t_min+0.01) & (tin{kk}<=t_max+0.01);
    tin{kk}=tin{kk}(cond)';
    pin{kk}=pin{kk}(cond)';
end

%-----original data
Pout.base=pin;

%-----depth matching
% empirically-determined best depth matches
sta_in= {'GNS21-PA','GNS21-PB','GNS21-PC','GNS21-PD','GNS21-PF','GNS21-PG',...
    'GNS21-PI','GNS21-PJ','TU21-PA','TU21-PB','TU21-PC','TU21-PD','TU21-PE',...
    'BPR','U1518','U1519'};
sta_ref1={'TU21-PE','TU21-PA','TU21-PA','TU21-PD','TU21-PE','TU21-PD',...
    'TU21-PD','TU21-PB','','','TU21-PA','','',...
    '','TU21-PB','TU21-PD'};
sta_ref2={'GNS21-PJ','TU21-PC','TU21-PC','GNS21-PI','TU21-PA','GNS21-PC',...
    'GNS21-PC','GNS21-PA','','','TU21-PE','','',...
    '','GNS21-PA','GNS21-PC'};
%-----network average
% calculate the average now for use below
for m=1:length(pin)
    pavg(:,m)=pin{m}/std(pin{m});
end
pavg=mean(pavg,2);

for k=1:length(pin)
    t1=tin{k};
    p1{k}=pin{k}';

    %-----network average
    pa{k}=pavg;

    Gp=pa{k};
    mp=inv(Gp'*Gp)*Gp'*p1{k};
    pa{k}=mp*pa{k};
    Pout.avg{k}=p1{k}-pa{k};

    %-----preferred reference
    i_match=find(strcmp(station{k},sta_in));
    ista_ref=find(strcmp(sta_ref1{i_match},station));
    if isempty(ista_ref)
        p2{k}=[];
        p3{k}=[];
        continue
    end
    p2{k}=pin{ista_ref}';

    Gp=p2{k};
    mp=inv(Gp'*Gp)*Gp'*p1{k};
    p2{k}=mp*p2{k};
    Pout.match1{k}=p1{k}-p2{k};

    %-----secondary reference
    i_match=find(strcmp(station{k},sta_in));
    ista_ref=find(strcmp(sta_ref2{i_match},station));
    if isempty(ista_ref)
        p3{k}=[];
        continue
    end
    p3{k}=pin{ista_ref}';

    Gp=p3{k};
    mp=inv(Gp'*Gp)*Gp'*p1{k};
    p3{k}=mp*p3{k};
    Pout.match2{k}=p1{k}-p3{k};
end

%% plotting to ensure it's all working as expected

n=ceil(length(pin)/2);

%-----preferred reference
figure(30); clf
subplot(121); hold on
for i=1:n
    if isempty(p1{i}) || isempty(p2{i})
        continue
    end
    plot(t1,p1{i}+10*(i-1),'color',[0 0.4470 0.7410],'linewidth',1)
    plot(t1,p2{i}+10*(i-1),'r')
    text(t1(end)+10,mean(p1{i}+10*(i-1)),station{i})
end
ylim([-10 n*10])
datetick('x',6)
box on; grid on
set(gca,'fontsize',14)
subplot(122); hold on
for i=n+1:length(pin)
    if isempty(p1{i}) || isempty(p2{i})
        continue
    end
    plot(t1,p1{i}+10*(i-n-1),'color',[0 0.4470 0.7410],'linewidth',1)
    plot(t1,p2{i}+10*(i-n-1),'r')
    text(t1(end)+10,mean(p1{i}+10*(i-n-1)),station{i})
end
ylim([-10 n*10])
datetick('x',6)
box on; grid on
set(gca,'fontsize',14)

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print('../figures/exploratory/HOBITSS_VIII/differences/proxies/depthmatch1_comparison_stack','-dpng','-r300')

% plot CEOF-corrected pressures
figure(31); clf
subplot(121); hold on
for i=1:n
    if isempty(p1{i}) || isempty(p2{i})
        continue
    end
    plot(t1,Pout.match1{i}+3*(i-1),'k','linewidth',1)
    text(t1(end)+10,mean(Pout.match1{i}+3*(i-1)),station{i})
end
ylim([-2 n*3])
datetick('x',6)
box on; grid on
set(gca,'fontsize',14)
subplot(122); hold on
for i=n+1:length(pin)
    if isempty(p1{i}) || isempty(p2{i})
        continue
    end
    plot(t1,Pout.match1{i}+3*(i-n-1),'k','linewidth',1)
    text(t1(end)+10,mean(Pout.match1{i}+3*(i-n-1)),station{i})
end
ylim([-2 n*3])
datetick('x',6)
box on; grid on
set(gca,'fontsize',14)

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print('../figures/exploratory/HOBITSS_VIII/differences/proxies/depthmatch1_difference_stack','-dpng','-r300')

%-----secondary reference
figure(32); clf
subplot(121); hold on
for i=1:n
    if isempty(p1{i}) || isempty(p3{i})
        continue
    end
    plot(t1,p1{i}+10*(i-1),'color',[0 0.4470 0.7410],'linewidth',1)
    plot(t1,p3{i}+10*(i-1),'r')
    text(t1(end)+10,mean(p1{i}+10*(i-1)),station{i})
end
ylim([-10 n*10])
datetick('x',6)
box on; grid on
set(gca,'fontsize',14)
subplot(122); hold on
for i=n+1:length(pin)
    if isempty(p1{i}) || isempty(p3{i})
        continue
    end
    plot(t1,p1{i}+10*(i-n-1),'color',[0 0.4470 0.7410],'linewidth',1)
    plot(t1,p3{i}+10*(i-n-1),'r')
    text(t1(end)+10,mean(p1{i}+10*(i-n-1)),station{i})
end
ylim([-10 n*10])
datetick('x',6)
box on; grid on
set(gca,'fontsize',14)

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print('../figures/exploratory/HOBITSS_VIII/differences/proxies/depthmatch2_comparison_stack','-dpng','-r300')

% plot CEOF-corrected pressures
figure(33); clf
subplot(121); hold on
for i=1:n
    if isempty(p1{i}) || isempty(p3{i})
        continue
    end
    plot(t1,Pout.match2{i}+3*(i-1),'k','linewidth',1)
    text(t1(end)+10,mean(Pout.match2{i}+3*(i-1)),station{i})
end
ylim([-2 n*3])
datetick('x',6)
box on; grid on
set(gca,'fontsize',14)
subplot(122); hold on
for i=n+1:length(pin)
    if isempty(p1{i}) || isempty(p3{i})
        continue
    end
    plot(t1,Pout.match2{i}+3*(i-n-1),'k','linewidth',1)
    text(t1(end)+10,mean(Pout.match2{i}+3*(i-n-1)),station{i})
end
ylim([-2 n*3])
datetick('x',6)
box on; grid on
set(gca,'fontsize',14)

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print('../figures/exploratory/HOBITSS_VIII/differences/proxies/depthmatch2_difference_stack','-dpng','-r300')

%-----network average
figure(34); clf
subplot(121); hold on
for i=1:n
    if isempty(p2{i}) || isempty(p3{i})
        continue
    end
    plot(t1,p1{i}+10*(i-1),'color',[0 0.4470 0.7410],'linewidth',1)
    plot(t1,pa{i}+10*(i-1),'r')
    text(t1(end)+10,mean(p1{i}+10*(i-1)),station{i})
end
ylim([-10 n*10])
datetick('x',6)
box on; grid on
set(gca,'fontsize',14)
subplot(122); hold on
for i=n+1:length(pin)
    if isempty(p2{i}) || isempty(p3{i})
        continue
    end
    plot(t1,p1{i}+10*(i-n-1),'color',[0 0.4470 0.7410],'linewidth',1)
    plot(t1,pa{i}+10*(i-n-1),'r')
    text(t1(end)+10,mean(p1{i}+10*(i-n-1)),station{i})
end
ylim([-10 n*10])
datetick('x',6)
box on; grid on
set(gca,'fontsize',14)

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print('../figures/exploratory/HOBITSS_VIII/differences/proxies/netavg_comparison_stack','-dpng','-r300')

% plot CEOF-corrected pressures
figure(35); clf
subplot(121); hold on
for i=1:n
    if isempty(p2{i}) || isempty(p3{i})
        continue
    end
    plot(t1,Pout.avg{i}+3*(i-1),'k','linewidth',1)
    text(t1(end)+10,mean(Pout.avg{i}+3*(i-1)),station{i})
end
ylim([-2 n*3])
datetick('x',6)
box on; grid on
set(gca,'fontsize',14)
subplot(122); hold on
for i=n+1:length(pin)
    if isempty(p2{i}) || isempty(p3{i})
        continue
    end
    plot(t1,Pout.avg{i}+3*(i-n-1),'k','linewidth',1)
    text(t1(end)+10,mean(Pout.avg{i}+3*(i-n-1)),station{i})
end
ylim([-2 n*3])
datetick('x',6)
box on; grid on
set(gca,'fontsize',14)

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print('../figures/exploratory/HOBITSS_VIII/differences/proxies/netavg_difference_stack','-dpng','-r300')

end