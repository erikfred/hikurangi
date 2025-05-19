function [poly_out,mp_out,sin_out,ms_out,lin_out,ml_out,ceof_out,mc_out,tc_out] = seasonal_correction(tin,pin,stain,yrstr)
%
% Applies each of two seasonal corrections, one sinusoidal and one
% polynomial
% 
% INPUTS:
%   tin - cell array of time vectors
%   pin - cell array of pressure vectors
%   stain - cell array of station IDs
%   yrstr - HOBITSS numeral
%
% OUTPUTS:
%   poly_out - polynomially-corrected pressures
%   mp_out  - polynomial model as 'm(1)*x^3 + m(2)*x^2 + m(3)*x + m(4) [+ m(5)*exp(-x/m(6))]'
%   sin_out - sinusoidally-corrected pressures
%   ms_out  - sinusoidal model as 'm(1)*sin(2*pi*t/365.25) + m(2)*cos(2*pi*t/365.25) + ...
%                                   m(3)*x + m(4) [+ m(5)*exp(-x/m(6))]'
%   lin_out - linearly-corrected pressures
%   ml_out  - linear model as 'm(1)*x + m(2) [+ m(3)*exp(-x/m(4))]'
%
%               Where   t = (t_in-t_in(1))
%                       x = t/max(t_in-t_in(1)) = (t_in-t_in(1))/max(t_in-t_in(1))
%               And the exponential terms are excluded when no significant
%               exponent is identified
%   ceof_out - ceof-corrected pressures
%   mc_out  - linear model as 'm(1)*x + m(2) [+ m(3)*exp(-x/m(4))]' (for residual)
%   tc_out  - truncated time used for CEOF calculation/correction
%

% save time by identifying where fits have already been found and saved
% if exist(['../processed_data/HOBITSS_' yrstr '_dedrifted.mat'],'file')
%     prior = load(['../processed_data/HOBITSS_' yrstr '_dedrifted.mat']);
%     poly_out = prior.pcor; mp_out = prior.pmod;
%     sin_out = prior.scor; ms_out = prior.smod;
%     lin_out = prior.lcor; ml_out = prior.lmod;
% 
%     l_list = find(~ismember(stain,prior.staname));
% else
%     warning('No previous file found. Confirm to continue.')
%     keyboard
    l_list = 1:length(pin);
% end

% ensure arrays are shaped as columns
for ll=1:length(pin)
    if size(tin{ll},1)<size(tin{ll},2)
        tin{ll}=tin{ll}';
    end
    if size(pin{ll},1)<size(pin{ll},2)
        pin{ll}=pin{ll}';
    end
end

% approximate CEOFs (no exp/lin cleanup)
% unify the time basis
t_min=max(cellfun(@(v)v(1),tin));
t_max=min(cellfun(@(v)v(end),tin));

for kk=1:length(pin)
    cond=(tin{kk}>=t_min+0.01) & (tin{kk}<=t_max+0.01);
    tc_out{kk}=tin{kk}(cond);
    ptemp{kk}=pin{kk}(cond); ptemp{kk}=detrend(ptemp{kk});
    dP(:,kk)=ptemp{kk}/std(ptemp{kk});
end
% calculate CEOFs from combined dataset
[~,loadings,pcs,~]=ceof(dP,1,true); % only need 1st CEOF
pred1=real(pcs(:,1) * conj(loadings(:,1)).');
scaled_pred1=mat2cell(pred1.*cellfun(@std,ptemp),length(dP),ones(1,length(pin)));
ceof_temp=cellfun(@minus,ptemp,scaled_pred1,'UniformOutput',false);

for iii=1:length(ceof_temp)
    for jjj=1:length(ceof_temp)
        cor_grid(iii,jjj)=xcorr(dP(:,iii),dP(:,jjj),0,'coeff');
    end
end

for l=l_list
    t=tin{l};
    t2=tc_out{l};
    p=pin{l};
    if isempty(p)
        continue
    else
        dim=size(t);
        if dim(1)<dim(2) % transpose as needed
            t=t';
            p=p';
        end
    end
    lamda_list=linspace(1/max(t-t(1)),180/max(t-t(1)),1000);
    tinv=(t-t(1))/max(t-t(1)); % better time basis
    tinv2=(t2-t2(1))/max(t2-t2(1));
    
    %----- combined exponential + 3rd degree polynomial
    for jj=1:1000
        %construct matrix for inversion
        Gp=[tinv.^3,tinv.^2,tinv,ones(size(tinv))];
        gexp=exp(-tinv/lamda_list(jj)); gexp(gexp<10^-7)=0;
        Gp=[Gp,gexp];
        % rcond(Gp'*Gp)
        mp{l}(:,jj)=inv(Gp'*Gp)*Gp'*p;
        ppe_fit{l}(:,jj)=Gp*mp{l}(:,jj);
        stds{l}(jj)=std(p-ppe_fit{l}(:,jj));
    end
    [~,imin]=min(stds{l});
    poly_out{l}=p-ppe_fit{l}(:,imin);
    mp_out{l}=[mp{l}(:,imin);lamda_list(imin)];
    % keyboard % intervene as necessary
    if false
        Gp=[tinv.^3,tinv.^2,tinv,ones(size(tinv))];
        mp_out{l}=inv(Gp'*Gp)*Gp'*p;
        poly_out{l}=p-Gp*mp_out{l};
    end

    %----- combined exponential + sinusoid
    yf=365.25;
    for jj=1:1000
        %construct matrix for inversion
        Gs=[sin(2*pi/yf*tinv*max(t-t(1))),cos(2*pi/yf*tinv*max(t-t(1))),tinv,ones(size(tinv))];
        gexp=exp(-tinv/lamda_list(jj)); gexp(gexp<10^-7)=0;
        Gs=[Gs,gexp];
        % rcond(Gs'*Gs)
        ms{l}(:,jj)=inv(Gs'*Gs)*Gs'*p;
        pse_fit{l}(:,jj)=Gs*ms{l}(:,jj);
        stds{l}(jj)=std(p-pse_fit{l}(:,jj));
    end
    [~,imin]=min(stds{l});
    sin_out{l}=p-pse_fit{l}(:,imin);
    ms_out{l}=[ms{l}(:,imin);lamda_list(imin)];
    % keyboard % intervene as necessary
    if false
        Gs=[sin(2*pi/yf*tinv*max(t-t(1))),cos(2*pi/yf*tinv*max(t-t(1))),tinv,ones(size(tinv))];
        ms_out{l}=inv(Gs'*Gs)*Gs'*p;
        sin_out{l}=p-Gs*ms_out{l};
    end

    %----- exponential only
    for jj=1:1000
        %construct matrix for inversion
        Gl=[tinv,ones(size(tinv))];
        gexp=exp(-tinv/lamda_list(jj)); gexp(gexp<10^-7)=0;
        Gl=[Gl,gexp];
        % rcond(Gl'*Gl)
        ml{l}(:,jj)=inv(Gl'*Gl)*Gl'*p;
        ple_fit{l}(:,jj)=Gl*ml{l}(:,jj);
        stds{l}(jj)=std(p-ple_fit{l}(:,jj));
    end
    [~,imin]=min(stds{l});
    lin_out{l}=p-ple_fit{l}(:,imin);
    ml_out{l}=[ml{l}(:,imin);lamda_list(imin)];
    % keyboard % intervene as necessary
    if false
        Gl=[tinv,ones(size(tinv))];
        ml_out{l}=inv(Gl'*Gl)*Gl'*p;
        lin_out{l}=p-Gl*ml_out{l};
    end

    %----- CEOF + exponential
    for jj=1:1000
        %construct matrix for inversion
        Gc=[tinv2,ones(size(tinv2))];
        gexp=exp(-tinv2/lamda_list(jj)); gexp(gexp<10^-7)=0;
        Gc=[Gc,gexp];
        % rcond(Gl'*Gl)
        mc{l}(:,jj)=inv(Gc'*Gc)*Gc'*ceof_temp{l};
        pce_fit{l}(:,jj)=Gc*mc{l}(:,jj);
        stds{l}(jj)=std(ceof_temp{l}-pce_fit{l}(:,jj));
    end
    [~,imin]=min(stds{l});
    ptemp2{l}=ptemp{l}-pce_fit{l}(:,imin);
    mc_out{l}=[mc{l}(:,imin);lamda_list(imin)];
    % figure(2); clf
    % subplot(211); hold on
    % plot(ceof_temp{l},'k','linewidth',1)
    % plot(ceof_temp{l}-pce_fit{l}(:,imin),'linewidth',1)
    % subplot(212); hold on
    % plot(stds{l},'o')
    % if l==3
        keyboard % intervene as necessary
    % end
    if false
        ptemp2{l}=ptemp{l};
        mc_out{l}=[];
    end
end

% finalize CEOFs
for kk=1:length(pin)
    dP2(:,kk)=ptemp2{kk}/std(ptemp2{kk});
end
% calculate CEOFs from combined dataset
[~,loadings,pcs,~]=ceof(dP2,1,true); % only need 1st CEOF
pred2=real(pcs(:,1) * conj(loadings(:,1)).');
scaled_pred2=mat2cell(pred2.*cellfun(@std,ptemp2),length(dP2),ones(1,length(pin)));
ceof_out=cellfun(@minus,ptemp2,scaled_pred2,'UniformOutput',false);

% plot CEOFs over pressures
n=ceil(length(pin)/2);
figure(99); clf
subplot(121); hold on
for i=1:n
    plot(tc_out{i},dP2(:,i)+4*(i-1),'r','linewidth',1)
    plot(tc_out{i},pred2(:,i)+4*(i-1),'k')
    text(tc_out{i}(end)+10,mean(ceof_out{i}+4*(i-1)),stain{i})
end
ylim([-2 n*4])
datetick('x',6)
box on; grid on
set(gca,'fontsize',14)
subplot(122); hold on
for i=n+1:length(pin)
    plot(tc_out{i},dP2(:,i)+4*(i-n-1),'r','linewidth',1)
    plot(tc_out{i},pred2(:,i)+4*(i-n-1),'k')
    text(tc_out{i}(end)+10,mean(ceof_out{i}+4*(i-n-1)),stain{i})
end
ylim([-2 n*4])
datetick('x',6)
box on; grid on
set(gca,'fontsize',14)

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
% print(['../figures/exploratory/HOBITSS_' yrstr '/differences/ceof/ceof_comparison_stack'],'-dpng','-r300')

% plot CEOF-corrected pressures
n=ceil(length(pin)/2);
figure(97); clf
subplot(121); hold on
for i=1:n
    plot(tc_out{i},ceof_out{i}+3*(i-1),'color',[0 0.4470 0.7410],'linewidth',1)
    text(tc_out{i}(end)+10,mean(ceof_out{i}+3*(i-1)),stain{i})
end
ylim([-2 n*3])
datetick('x',6)
box on; grid on
set(gca,'fontsize',14)
subplot(122); hold on
for i=n+1:length(pin)
    plot(tc_out{i},ceof_out{i}+3*(i-n-1),'color',[0 0.4470 0.7410],'linewidth',1)
    text(tc_out{i}(end)+10,mean(ceof_out{i}+3*(i-n-1)),stain{i})
end
ylim([-2 n*3])
datetick('x',6)
box on; grid on
set(gca,'fontsize',14)

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
% print(['../figures/exploratory/HOBITSS_' yrstr '/differences/ceof/ceof_corrected_stack'],'-dpng','-r300')

% plot pressures at full scale
n=ceil(length(pin)/2);
figure(98); clf
subplot(121); hold on
for i=1:n
    plot(tc_out{i},ptemp2{i}+10*(i-1),'color',[0 0.4470 0.7410],'linewidth',1)
    text(tc_out{i}(end)+10,mean(ceof_out{i}+10*(i-1)),stain{i})
end
ylim([-10 n*10])
datetick('x',6)
box on; grid on
set(gca,'fontsize',14)
subplot(122); hold on
for i=n+1:length(pin)
    plot(tc_out{i},ptemp2{i}+10*(i-n-1),'color',[0 0.4470 0.7410],'linewidth',1)
    text(tc_out{i}(end)+10,mean(ceof_out{i}+10*(i-n-1)),stain{i})
end
ylim([-10 n*10])
datetick('x',6)
box on; grid on
set(gca,'fontsize',14)

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];

keyboard