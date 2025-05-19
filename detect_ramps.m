function [detections]=detect_ramps(year_index)
%
% Detect and quantify ramps in HOBITSS data.
%

yr_str={'I','III','IV','Vb','VI','VII','VIII'};
i=year_index;
yrstr=yr_str{i};

load(['../processed_data/HOBITSS_' yrstr '_dedrifted.mat'])

dt=1/24; % sample rate
rampdur=[7,14,30,60,90]*24; % assumed ramp duration [hr]
knowntime=false; % true assumes known event timing

%% Organization

t_min=min([tf{:}]); t_max=max([tf{:}]); t_tot=t_min:1/24:t_max;

% add data from CORK wellheads, if applicable
if i==3
    cork = load('../processed_data/CORK_dedrifted.mat');
    stadepth=cat(2,stadepth,cork.stadepth(2:3));
    stalat=cat(2,stalat,cork.stalat(2:3));
    stalon=cat(2,stalon,cork.stalon(2:3));
    staname=cat(2,staname,cork.staname(2:3));

    % determine temporal overlap
    cond=(cork.tf{2}>=t_min) & (cork.tf{2}<=t_max);
    tf=cat(2,tf,{cork.tf{2}(cond)});
    pf=cat(2,pf,{cork.pf{2}(cond)});
    cond=(cork.tf{3}>=t_min) & (cork.tf{3}<=t_max);
    tf=cat(2,tf,{cork.tf{3}(cond)});
    pf=cat(2,pf,{cork.pf{3}(cond)});
    pmod=[pmod,{[0,0,0,0]},{[0,0,0,0]}]; smod=[smod,{[0,0,0,0]},{[0,0,0,0]}];
    lmod=[lmod,{[0,0]},{[0,0]}]; cmod=[cmod,{[]},{[]}];
elseif i>3 % wellheads + BPR online
    cork = load('../processed_data/CORK_dedrifted.mat');
    stadepth=cat(2,stadepth,cork.stadepth);
    stalat=cat(2,stalat,cork.stalat);
    stalon=cat(2,stalon,cork.stalon);
    staname=cat(2,staname,cork.staname);

    % determine temporal overlap
    cond=(cork.tf{1}>=t_min) & (cork.tf{1}<=t_max);
    tf=cat(2,tf,{cork.tf{1}(cond)});
    pf=cat(2,pf,{cork.pf{1}(cond)});
    cond=(cork.tf{2}>=t_min) & (cork.tf{2}<=t_max);
    tf=cat(2,tf,{cork.tf{2}(cond)});
    pf=cat(2,pf,{cork.pf{2}(cond)});
    cond=(cork.tf{3}>=t_min) & (cork.tf{3}<=t_max);
    tf=cat(2,tf,{cork.tf{3}(cond)});
    pf=cat(2,pf,{cork.pf{3}(cond)});
    pmod=[pmod,{[0,0,0,0]},{[0,0,0,0]},{[0,0,0,0]}];
    smod=[smod,{[0,0,0,0]},{[0,0,0,0]},{[0,0,0,0]}];
    lmod=[lmod,{[0,0]},{[0,0]},{[0,0]}]; cmod=[cmod,{[]},{[]},{[]}];
end

% sort by depth
[depth,id]=sort(stadepth);
name=staname(id);
lon=stalon(id);
lat=stalat(id);
tf=tf(id);
pf=pf(id);
pm=pmod(id);
sm=smod(id);
lm=lmod(id);
cm=cmod(id);

%% detect ramps over a range of parameters

% preallocate detection structure for improved speed
var1={'lin','poly','sin','ceof'};
for i1=1:4
    var2={'base','match','avg'};
    for i2=1:3
        structstr=['detections.' var1{i1} '.' var2{i2}];
        eval([structstr '=struct(''amppred'',NaN(length(name),length(t_tot)),'...
            '''stds'',NaN(length(name),length(t_tot)))'])
    end
end

for sta=1:length(name)
    for onset=30*24:24:length(pf{sta})-30*24
        % cycle over ramp durations
        for i=1:length(rampdur)
            if onset+rampdur(i)>length(pf{sta})-30*24
                disp('Insufficient data after onset to assess')
                continue
            end
            if knowntime
                dt=[dt(1);find(ramptim(onset)+shdur/2==t_hrf{1})];
            end
            
            % apply seasonal corrections
            [~,ppe]=ssn_cor_smplexp(tf,pf,name,'poly');
            [~,pse]=ssn_cor_smplexp(tf,pf,name,'sin');
            [~,ple]=ssn_cor_smplexp(tf,pf,name,'lin');
            [~,pce]=ssn_cor_smplexp(tf,pf,name,'ceof');

            % apply proxy corrections
            p_corr.poly=proxy_corrections(tf,ppe,name,sta);
            p_corr.sin=proxy_corrections(tf,pse,name,sta);
            p_corr.lin=proxy_corrections(tf,ple,name,sta);
            p_corr.ceof=proxy_corrections(tf,pce,name,sta);

                %%%%%%%%%%%
                % once differenced (or otherwise corrected), this snippet
                % steps through onsets, masks out the presumed deformation,
                % and finds the best seasonal fit. The goodness of fit for
                % each theoretical ramp provides a measure of the ramp's
                % presence -- masking out an area without a ramp will not
                % significantly improve the fit vs. no masking. The ramp
                % size is determined from the difference in constant of
                % offset from either side of the masked interval.
                t=ttempf(1200:end)'; p=H8.pf{4}(iaf)-H8.pf{12}(ibf); p=p(1200:end)';
                yf=365.25;
                for jj=1:1000
                    t1=t(1:5*jj); t2=t(5*jj+1000:end);
                    p1=p(1:5*jj); p2=p(5*jj+1000:end);

                    tinv=(t-t(1))/max(t-t(1)); % better time basis
                    tinv1=tinv(1:5*jj); tinv2=tinv(5*jj+1000:end);

                    tinv=[tinv1;tinv2];
                    pinv=[p1;p2];

                    % polynomial inversion
                    Gp=[tinv.^3,tinv.^2,tinv,...
                        [ones(size(tinv1));zeros(size(tinv2))],[zeros(size(tinv1));ones(size(tinv2))]];
                    mp(:,jj)=inv(Gp'*Gp)*Gp'*pinv;
                    ppe_fit(:,jj)=Gp*mp(:,jj);
                    step(jj)=mp(4,jj)-mp(5,jj);
                    stds(jj)=std(pinv-ppe_fit(:,jj));

                    % sinusoidal inversion
                    Gs=[sin(2*pi/yf*tinv*max(t-t(1))),cos(2*pi/yf*tinv*max(t-t(1))),tinv,...
                        [ones(size(tinv1));zeros(size(tinv2))],[zeros(size(tinv1));ones(size(tinv2))]];
                    ms(:,jj)=inv(Gs'*Gs)*Gs'*pinv;
                    pse_fit(:,jj)=Gs*ms(:,jj);
                    step2(jj)=ms(4,jj)-ms(5,jj);
                    stds2(jj)=std(pinv-pse_fit(:,jj));

                    % linear inversion
                    Gl=[tinv,[ones(size(tinv1));zeros(size(tinv2))],[zeros(size(tinv1));ones(size(tinv2))]];
                    ml(:,jj)=inv(Gl'*Gl)*Gl'*pinv;
                    ple_fit(:,jj)=Gl*ml(:,jj);
                    step3(jj)=mp(2,jj)-mp(3,jj);
                    stds3(jj)=std(pinv-ple_fit(:,jj));

                    % use linear version on ceof data? should work, as best
                    % fit will still be determined by properly masked offset
                end
                [~,imin]=min(stds);
                poly_out=pinv-ppe_fit(:,imin);
                mp_out=mp(:,imin);

                [~,imin2]=min(stds2);
                sin_out=p-pse_fit(:,imin2);
                ms_out=ms(:,imin2);

                [~,imin3]=min(stds3);
                lin_out=p-ple_fit(:,imin3);
                ml_out=ml(:,imin3);
                %%%%%%%%%%%%%%%

                % attempt to detect ramp
                %--linear
                [detections.lin.base.amppred(sta,onset),detections.lin.base.stds(sta,onset)]...
                    =findoffset_lin(p_corr.lin.base{sta},dt);
                [detections.lin.match.amppred(sta,onset),detections.lin.match.stds(sta,onset)]...
                    =findoffset_lin(p_corr.lin.match{sta},dt);
                [detections.lin.avg.amppred(sta,onset),detections.lin.avg.stds(sta,onset)]...
                    =findoffset_lin(p_corr.lin.avg{sta},dt);
                %--polynomial
                [detections.poly.base.amppred(sta,onset),detections.poly.base.stds(sta,onset)]...
                    =findoffset_poly(p_corr.poly.base{sta},dt);
                [detections.poly.match.amppred(sta,onset),detections.poly.match.stds(sta,onset)]...
                    =findoffset_poly(p_corr.poly.match{sta},dt);
                [detections.poly.avg.amppred(sta,onset),detections.poly.avg.stds(sta,onset)]...
                    =findoffset_poly(p_corr.poly.avg{sta},dt);
                %--sinusoid
                [detections.sin.base.amppred(sta,onset),detections.sin.base.stds(sta,onset)]...
                    =findoffset_sin(p_corr.sin.base{sta},dt);
                [detections.sin.match.amppred(sta,onset),detections.sin.match.stds(sta,onset)]...
                    =findoffset_sin(p_corr.sin.match{sta},dt);
                [detections.sin.avg.amppred(sta,onset),detections.sin.avg.stds(sta,onset)]...
                    =findoffset_sin(p_corr.sin.avg{sta},dt);
                %--CEOF1
                [detections.ceof.base.amppred(sta,onset),detections.ceof.base.stds(sta,onset)]...
                    =findoffset_ceof(p_corr.ceof.base{sta},dt);
                [detections.ceof.match.amppred(sta,onset),detections.ceof.match.stds(sta,onset)]...
                    =findoffset_ceof(p_corr.ceof.match{sta},dt);
                [detections.ceof.avg.amppred(sta,onset),detections.ceof.avg.stds(sta,onset)]...
                    =findoffset_ceof(p_corr.ceof.avg{sta},dt);
        end
    end
end

% add in 'key' to stations used
%--amplitude
detections.lin.base.staname=name;detections.lin.match.staname=name;detections.lin.avg.staname=name;
detections.poly.base.staname=name;detections.poly.match.staname=name;detections.poly.avg.staname=name;
detections.sin.base.staname=name;detections.sin.match.staname=name;detections.sin.avg.staname=name;
detections.ceof.base.staname=name;detections.ceof.match.staname=name;detections.ceof.avg.staname=name;

keyboard
end
%% ---------INTERNAL FUNCTIONS----------

% generate ramp with specified features
function [pout]=make_ramp(tin,pin,rampamp,rampdur,ramptim)

for i=1:length(tin)
    rmp=linspace(0,rampamp,rampdur*24)';
    [~,imin]=min(abs(tin{i}-ramptim));
    sse_ramp=[zeros(imin-1,1);rmp;rampamp*ones(length(tin{i})-(length(rmp)+imin-1),1)];
    sse_ramp=sse_ramp-mean(sse_ramp);
    pout{i}=pin{i}+sse_ramp;
end

end

% detect and quantify step by assuming form linear->ramp->linear
function [ipt,rampamp,sigflag]=quantify_step(X,A,sname,ind,dt,known_amp,proxy)
X=X-mean(X); % remove any residual average
t=[1:length(X)]';
eventlength=1; % in hours

% detect ramp by usual method
if length(dt)==1
    for jj=1:length(eventlength)
        for ii=24*30+1:24:length(t)-24*30 % don't allow detections in first and last month
            xseg1=t(1:ii);
            xseg3=t(ii+eventlength(jj):end);
            xseg=t([1:ii,ii+eventlength(jj):end]);
            yseg1=X(1:ii);
            yseg3=X(ii+eventlength(jj):end);
            
            % set up inversion
            G=[xseg,[ones(size(yseg1));zeros(size(yseg3))],...
                [zeros(size(yseg1));ones(size(yseg3))]];
            m=inv(G'*G)*G'*[yseg1;yseg3];
            
            % correct timeseries w/ determined slope
            yseg1_corr=yseg1-polyval([m(1),m(2)],xseg1);
            yseg3_corr=yseg3-polyval([m(1),m(3)],xseg3);
            rmslist(ii,jj)=std([yseg1_corr;yseg3_corr]);
        end
    end
    rmslist(rmslist==0)=1000;
    [~,ipt]=min(rmslist);
else
    ipt=dt(2);
    dt=dt(1);
end

% re-calculate CEOFs, with modeled deformation ramp, selecting that with best fit
% correct out sensor drift
% better time basis
temp1=(t-t(1))/max(t-t(1));
if length(sname)>10
    lamda_list=linspace(1/max(t-t(1)),180/max(t-t(1)),1000)*24;
else
    lamda_list=linspace(1/max(t-t(1)),720/max(t-t(1)),1000)*24;
end
eventlength=14*24; % in hours
halflength=eventlength/2;

for k=1:length(A)
    % extract exponential time constant (empirically determined)
    staname2={'LA21';'LA22';'LA23';'LA25';'LA26';'LA28';'LA30';'LA32';'LA33';'LA34';'LA39';...
        'LT01';'LT02';'LT03';'LT04';'LT06';'LT07';'LT09';'LT10';'LT11';'LT12';'LT13';'LT14';'LT16';'LT20'};
%     k_list_man=[1000,1000,157,26,1000,28,1000,786,1000,1000,1000,1000,0.00,...
%         559,33,373,1000,0.00,1000,1000,1000,861,1000,876,144];
%     k_list_amp=[27.67,27.43,-10.39,-2.71,0.00,-4.69,8.52,87.39,30.09,33.34,-91.77,-9.98,0.00,...
%         -15.31,-4.60,-3.53,-11.21,0.00,-9.86,200.60,-11.93,-3.83,-18.80,16.75,-3.94];
    k_list_man=[1,311,53,19,NaN,22,1,162,315,479,1000,...
        1,NaN,318,36,254,288,NaN,1,964,1,130,400,502,146];
    k_list_amp=[0,14.79,-10.49,-2.50,NaN,-4.69,0,72.56,36.51,73.69,-763.67,...
        0,NaN,-12.94,-4.36,-4.80,-5.40,NaN,0,191.80,0,-1.75,-10.31,10.57,-5.15];
    
    ind2=find(strcmp(staname2,sname{k}));
    lamda=lamda_list(k_list_man(ind2));
    amp=k_list_amp(ind2);
    expfit=amp*exp(-temp1/lamda);
    
    pcor{k}=detrend(A{k}-expfit);
    pfit{k}=A{k}-pcor{k};
end
% assemble normalized matrix of P, attempting to model ramp
if size(pcor,2)==10
    proxies=load('../pressure_matfiles/JoanGomberg/SlopesDataCEOF.mat');
else
    proxies=load('../pressure_matfiles/JoanGomberg/ShelfDataCEOF.mat');
end

% faster method that iteratively decreases step size of search
loopcheck1=0;
minguess=-20;
maxguess=20;
guessstep=10;
guessrange=minguess:guessstep:maxguess;
while loopcheck1<=2
    fit_check=[];
    loopcheck1=loopcheck1+1;
    for mm=1:length(guessrange)
        amp_guess=guessrange(mm)+known_amp;
        ramp_guess=zeros(size(X));
        ramp_guess(ipt-halflength:ipt+halflength-1)=[linspace(0,amp_guess,eventlength)]';
        ramp_guess(ipt+halflength:end)=ramp_guess(ipt+halflength:end)+amp_guess;
        ptemp=pcor;
        ptemp{ind}=detrend(ptemp{ind}-ramp_guess);
        for kk=1:length(ptemp)-1
            jj=find(strcmp(proxies.staname,sname{kk}));
            dP(:,kk)=pcor{kk}/std(pcor{kk});
            dT(:,kk)=proxies.datatemp(:,jj)/std(proxies.datatemp(:,jj));
            dS(:,kk)=proxies.datassh(:,jj)/std(proxies.datassh(:,jj));
        end
        % calculate CEOFs from combined dataset
        [~,loadings,pcs,~]=ceof([dP,dT,dS],1,true); % only need 1st CEOF
        pred1=real(pcs(:,1) * conj(loadings(:,1)).');
        % apply CEOF "seasonal" correction
        scaled_pred1=mat2cell(pred1(:,1:length(ptemp)).*cellfun(@std,ptemp),9047,ones(1,length(ptemp)));
        ptemp_cor=cellfun(@minus,ptemp,scaled_pred1,'UniformOutput',false);
        % apply proxy corrections
        ptemp_cor2=proxy_corrections(ptemp,ptemp_cor,sname,'ceof',ind,false);
        
        if strcmp(proxy,'base')
            pcheck=detrend(ptemp_cor2.base{ind});
        elseif strcmp(proxy,'ref')
            pcheck=detrend(ptemp_cor2.ref{ind});
        elseif strcmp(proxy,'match')
            pcheck=detrend(ptemp_cor2.match{ind});
        elseif strcmp(proxy,'avg')
            pcheck=detrend(ptemp_cor2.avg{ind});
        elseif strcmp(proxy,'temp')
            pcheck=detrend(ptemp_cor2.temp{ind});
        elseif strcmp(proxy,'ssh')
            pcheck=detrend(ptemp_cor2.ssh{ind});
        else
            warning('Proxy not properly indicated for CEOF calculation')
        end
        fit_check(mm)=std(pcheck); % will minimize RMS of proxy-corrected pressure
    end
    if loopcheck1~=3
        [~,tempind]=mink(fit_check,2);
        minguess=min(guessrange(tempind));
        maxguess=max(guessrange(tempind));
        guessstep=guessstep/10;
        guessrange=minguess:guessstep:maxguess;
    end
end

% use best fit CEOF calculation
[~,mm]=min(fit_check);
amp_guess=guessrange(mm)+known_amp;
ramp_guess=zeros(size(X));
ramp_guess(ipt-halflength:ipt+halflength-1)=[linspace(0,amp_guess,eventlength)]';
ramp_guess(ipt+halflength:end)=ramp_guess(ipt+halflength:end)+amp_guess;
ptemp=pcor;
ptemp{ind}=detrend(ptemp{ind}-ramp_guess);
for kk=1:length(ptemp)-1
    jj=find(strcmp(proxies.staname,sname{kk}));
    dP(:,kk)=pcor{kk}/std(pcor{kk});
    dT(:,kk)=proxies.datatemp(:,jj)/std(proxies.datatemp(:,jj));
    dS(:,kk)=proxies.datassh(:,jj)/std(proxies.datassh(:,jj));
end
% calculate CEOFs from combined dataset
[~,loadings,pcs,~]=ceof([dP,dT,dS],1,true); % only need 1st CEOF
pred1=real(pcs(:,1) * conj(loadings(:,1)).');
% apply best CEOF "seasonal" correction
scaled_pred1=mat2cell(pred1(:,1:length(ptemp)).*cellfun(@std,ptemp),9047,ones(1,length(ptemp)));
pcor_cor=cellfun(@minus,pcor,scaled_pred1,'UniformOutput',false);
% apply proxy corrections
pcor_cor2=proxy_corrections(pcor,pcor_cor,sname,'ceof',ind,false);

if strcmp(proxy,'base')
    X2=pcor_cor2.base{ind}; X2=X2-mean(X2);
elseif strcmp(proxy,'ref')
    X2=pcor_cor2.ref{ind}; X2=X2-mean(X2);
elseif strcmp(proxy,'match')
    X2=pcor_cor2.match{ind}; X2=X2-mean(X2);
elseif strcmp(proxy,'avg')
    X2=pcor_cor2.avg{ind}; X2=X2-mean(X2);
elseif strcmp(proxy,'temp')
    X2=pcor_cor2.temp{ind}; X2=X2-mean(X2);
elseif strcmp(proxy,'ssh')
    X2=pcor_cor2.ssh{ind}; X2=X2-mean(X2);
else
    warning('Proxy not properly indicated for CEOF calculation')
end

% determine ramp amplitude for newly-corrected pressure
for kk=1:length(ipt)
    xseg=t([1:ipt(kk)-halflength,ipt(kk)+halflength+1:end]);
    yseg1=X2(1:ipt(kk)-halflength);
    yseg3=X2(ipt(kk)+halflength+1:end);
    
    % set up inversion
    G=[xseg,[ones(size(yseg1));zeros(size(yseg3))],...
        [zeros(size(yseg1));ones(size(yseg3))]];
    m=inv(G'*G)*G'*[yseg1;yseg3];
    mdl(:,kk)=G*m;
    
    rampamp(kk)=m(3)-m(2);
    
    % correct timeseries w/ determined slope
    linpred=0:m(1):(length(X2)-1)*m(1);
    Xcorr=X2-linpred';
end

if false
    figure(13); hold on
    plot(t,X2+1*15,'k','linewidth',1)
    plot(t(xseg),G*m+1*15,'k','linewidth',0.5)
    set(gca,'fontsize',16)
    box on; grid minor
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 11 8.5];
end

sigflag=getsig(Xcorr,dt,rampamp,ipt);

end

function [ipt,rampamp,sigflag]=quantify_step2(X,A,sname,ind,dt,known_amp,proxy)
X=X-mean(X); % remove any residual average
t=[1:length(X)]';
eventlength=1; % in hours

% detect ramp by usual method
if length(dt)==1
    for jj=1:length(eventlength)
        for ii=24*30+1:24:length(t)-24*30 % don't allow detections in first and last month
            xseg1=t(1:ii);
            xseg3=t(ii+eventlength(jj):end);
            xseg=t([1:ii,ii+eventlength(jj):end]);
            yseg1=X(1:ii);
            yseg3=X(ii+eventlength(jj):end);
            
            % set up inversion
            G=[xseg,[ones(size(yseg1));zeros(size(yseg3))],...
                [zeros(size(yseg1));ones(size(yseg3))]];
            m=inv(G'*G)*G'*[yseg1;yseg3];
            
            % correct timeseries w/ determined slope
            yseg1_corr=yseg1-polyval([m(1),m(2)],xseg1);
            yseg3_corr=yseg3-polyval([m(1),m(3)],xseg3);
            rmslist(ii,jj)=std([yseg1_corr;yseg3_corr]);
        end
    end
    rmslist(rmslist==0)=1000;
    [~,ipt]=min(rmslist);
else
    ipt=dt(2);
    dt=dt(1);
end

% re-calculate CEOFs, with modeled deformation ramp, selecting that with best fit
% correct out sensor drift
% better time basis
temp1=(t-t(1))/max(t-t(1));
if length(sname)>10
    lamda_list=linspace(1/max(t-t(1)),180/max(t-t(1)),1000)*24;
else
    lamda_list=linspace(1/max(t-t(1)),720/max(t-t(1)),1000)*24;
end
eventlength=14*24; % in hours
halflength=eventlength/2;

% correct out sensor drift
for k=1:length(A)
    % extract exponential time constant (empirically determined)
    staname2={'LA21';'LA22';'LA23';'LA25';'LA26';'LA28';'LA30';'LA32';'LA33';'LA34';'LA39';...
        'LT01';'LT02';'LT03';'LT04';'LT06';'LT07';'LT09';'LT10';'LT11';'LT12';'LT13';'LT14';'LT16';'LT20'};
%     k_list_man=[1000,1000,157,26,1000,28,1000,786,1000,1000,1000,1000,0.00,...
%         559,33,373,1000,0.00,1000,1000,1000,861,1000,876,144];
%     k_list_amp=[27.67,27.43,-10.39,-2.71,0.00,-4.69,8.52,87.39,30.09,33.34,-91.77,-9.98,0.00,...
%         -15.31,-4.60,-3.53,-11.21,0.00,-9.86,200.60,-11.93,-3.83,-18.80,16.75,-3.94];
    k_list_man=[1,311,53,19,NaN,22,1,162,315,479,1000,...
        1,NaN,318,36,254,288,NaN,1,964,1,130,400,502,146];
    k_list_amp=[0,14.79,-10.49,-2.50,NaN,-4.69,0,72.56,36.51,73.69,-763.67,...
        0,NaN,-12.94,-4.36,-4.80,-5.40,NaN,0,191.80,0,-1.75,-10.31,10.57,-5.15];
    
    ind2=find(strcmp(staname2,sname{k}));
    lamda=lamda_list(k_list_man(ind2));
    amp=k_list_amp(ind2);
    expfit=amp*exp(-temp1/lamda);
    
    pcor{k}=detrend(A{k}-expfit);
    pfit{k}=A{k}-pcor{k};
end

% assemble normalized matrix of P, attempting to model ramp

% faster method that iteratively decreases step size of search
loopcheck1=0;
minguess=-20;
maxguess=20;
guessstep=10;
guessrange=minguess:guessstep:maxguess;
while loopcheck1<=2
    fit_check=[];
    loopcheck1=loopcheck1+1;
    for mm=1:length(guessrange)
        amp_guess=guessrange(mm)+known_amp;
        ramp_guess=zeros(size(X));
        ramp_guess(ipt-halflength:ipt+halflength-1)=[linspace(0,amp_guess,eventlength)]';
        ramp_guess(ipt+halflength:end)=ramp_guess(ipt+halflength:end)+amp_guess;
        ptemp=pcor;
        ptemp{ind}=detrend(ptemp{ind}-ramp_guess);
        for kk=1:length(ptemp)-1 % keeps reference station out
            dP(:,kk)=ptemp{kk}/std(ptemp{kk});
        end
        % calculate CEOFs from combined dataset
        [~,loadings,pcs,~]=ceof(dP,1,true); % only need 1st CEOF
        pred1=real(pcs(:,1) * conj(loadings(:,1)).');
        % apply CEOF "seasonal" correction
        scaled_pred1=mat2cell(pred1(:,1:length(ptemp)-1).*cellfun(@std,ptemp(1:end-1)),9047,ones(1,length(ptemp)-1));
        ptemp_cor=cellfun(@minus,ptemp(1:end-1),scaled_pred1,'UniformOutput',false);
        ptemp_cor(length(ptemp))=ptemp(end);
        % apply proxy corrections
        ptemp_cor2=proxy_corrections(ptemp,ptemp_cor,sname,'ceof',ind,false);
        
        if strcmp(proxy,'base')
            pcheck=detrend(ptemp_cor2.base{ind});
        elseif strcmp(proxy,'ref')
            pcheck=detrend(ptemp_cor2.ref{ind});
        elseif strcmp(proxy,'match')
            pcheck=detrend(ptemp_cor2.match{ind});
        elseif strcmp(proxy,'avg')
            pcheck=detrend(ptemp_cor2.avg{ind});
        elseif strcmp(proxy,'temp')
            pcheck=detrend(ptemp_cor2.temp{ind});
        elseif strcmp(proxy,'ssh')
            pcheck=detrend(ptemp_cor2.ssh{ind});
        else
            warning('Proxy not properly indicated for CEOF calculation')
        end
        fit_check(mm)=std(pcheck); % will minimize RMS of proxy-corrected pressure
    end
    if loopcheck1~=3
        [~,tempind]=mink(fit_check,2);
        minguess=min(guessrange(tempind));
        maxguess=max(guessrange(tempind));
        guessstep=guessstep/10;
        guessrange=minguess:guessstep:maxguess;
    end
end

% use best fit CEOF calculation
[~,mm]=min(fit_check);
amp_guess=guessrange(mm)+known_amp;
ramp_guess=zeros(size(X));
ramp_guess(ipt-halflength:ipt+halflength-1)=[linspace(0,amp_guess,eventlength)]';
ramp_guess(ipt+halflength:end)=ramp_guess(ipt+halflength:end)+amp_guess;
ptemp=pcor;
ptemp{ind}=detrend(ptemp{ind}-ramp_guess);
for kk=1:length(ptemp)-1 % keeps reference station out
    dP(:,kk)=ptemp{kk}/std(ptemp{kk});
end
% calculate CEOFs from combined dataset
[~,loadings,pcs,~]=ceof(dP,1,true); % only need 1st CEOF
pred1=real(pcs(:,1) * conj(loadings(:,1)).');
% apply best CEOF "seasonal" correction
scaled_pred1=mat2cell(pred1(:,1:length(ptemp)-1).*cellfun(@std,ptemp(1:end-1)),9047,ones(1,length(ptemp)-1));
pcor_cor=cellfun(@minus,pcor(1:end-1),scaled_pred1,'UniformOutput',false);
pcor_cor(length(pcor))=pcor(end);
% apply proxy corrections
pcor_cor2=proxy_corrections(pcor,pcor_cor,sname,'ceof',ind,false);

if strcmp(proxy,'base')
    X2=pcor_cor2.base{ind}; X2=X2-mean(X2);
elseif strcmp(proxy,'ref')
    X2=pcor_cor2.ref{ind}; X2=X2-mean(X2);
elseif strcmp(proxy,'match')
    X2=pcor_cor2.match{ind}; X2=X2-mean(X2);
elseif strcmp(proxy,'avg')
    X2=pcor_cor2.avg{ind}; X2=X2-mean(X2);
elseif strcmp(proxy,'temp')
    X2=pcor_cor2.temp{ind}; X2=X2-mean(X2);
elseif strcmp(proxy,'ssh')
    X2=pcor_cor2.ssh{ind}; X2=X2-mean(X2);
else
    warning('Proxy not properly indicated for CEOF calculation')
end

% determine ramp amplitude for newly-corrected pressure
for kk=1:length(ipt)
    xseg=t([1:ipt(kk)-halflength,ipt(kk)+halflength+1:end]);
    yseg1=X2(1:ipt(kk)-halflength);
    yseg3=X2(ipt(kk)+halflength+1:end);
    
    % set up inversion
    G=[xseg,[ones(size(yseg1));zeros(size(yseg3))],...
        [zeros(size(yseg1));ones(size(yseg3))]];
    m=inv(G'*G)*G'*[yseg1;yseg3];
    mdl(:,kk)=G*m;
    
    rampamp(kk)=m(3)-m(2);
    
    % correct timeseries w/ determined slope
    linpred=0:m(1):(length(X2)-1)*m(1);
    Xcorr=X2-linpred';
end

if false
    figure(13); clf; hold on
    plot(t,X2,'k','linewidth',1)
    plot(t(xseg),G*m,'k','linewidth',0.5)
end

sigflag=getsig(Xcorr,dt,rampamp,ipt);

end

% CEOF-based step detection
function [ipt,better_ceof,rampamp,sigflag]=findoffset_ceof(X,A,sname,ind,dt,varargin)
X=X-mean(X); % remove any residual average
t=[1:length(X)]';
iref=find(strcmp(sname,'LA21')); % find reference station index
A(iref)=[]; % remove reference station
sname(iref)=[]; % remove reference station
eventlength=1; % in hours

% detect ramp by usual method
if length(dt)==1
    for jj=1:length(eventlength)
        for ii=24*30+1:24:length(t)-24*30 % don't allow detections in first and last month
            xseg1=t(1:ii);
            xseg3=t(ii+eventlength(jj):end);
            xseg=t([1:ii,ii+eventlength(jj):end]);
            yseg1=X(1:ii);
            yseg3=X(ii+eventlength(jj):end);
            
            % set up inversion
            G=[xseg,[ones(size(yseg1));zeros(size(yseg3))],...
                [zeros(size(yseg1));ones(size(yseg3))]];
            m=inv(G'*G)*G'*[yseg1;yseg3];
            
            % correct timeseries w/ determined slope
            yseg1_corr=yseg1-polyval([m(1),m(2)],xseg1);
            yseg3_corr=yseg3-polyval([m(1),m(3)],xseg3);
            rmslist(ii,jj)=std([yseg1_corr;yseg3_corr]);
        end
    end
    rmslist(rmslist==0)=1000;
    [~,ipt]=min(rmslist);
else
    ipt=dt(2);
    dt=dt(1);
end

% re-calculate CEOFs, with modeled deformation ramp, selecting that with best fit
% correct out sensor drift
% better time basis
temp1=(t-t(1))/max(t-t(1));
if length(sname)>10
    lamda_list=linspace(1/max(t-t(1)),180/max(t-t(1)),1000)*24;
else
    lamda_list=linspace(1/max(t-t(1)),720/max(t-t(1)),1000)*24;
end
eventlength=14*24; % in hours
halflength=eventlength/2;

for k=1:length(A)
    % extract exponential time constant (empirically determined)
    staname2={'LA21';'LA22';'LA23';'LA25';'LA26';'LA28';'LA30';'LA32';'LA33';'LA34';'LA39';...
        'LT01';'LT02';'LT03';'LT04';'LT06';'LT07';'LT09';'LT10';'LT11';'LT12';'LT13';'LT14';'LT16';'LT20'};
%     k_list_man=[1000,1000,157,26,1000,28,1000,786,1000,1000,1000,1000,0.00,...
%         559,33,373,1000,0.00,1000,1000,1000,861,1000,876,144];
%     k_list_amp=[27.67,27.43,-10.39,-2.71,0.00,-4.69,8.52,87.39,30.09,33.34,-91.77,-9.98,0.00,...
%         -15.31,-4.60,-3.53,-11.21,0.00,-9.86,200.60,-11.93,-3.83,-18.80,16.75,-3.94];
    k_list_man=[1,311,53,19,NaN,22,1,162,315,479,1000,...
        1,NaN,318,36,254,288,NaN,1,964,1,130,400,502,146];
    k_list_amp=[0,14.79,-10.49,-2.50,NaN,-4.69,0,72.56,36.51,73.69,-763.67,...
        0,NaN,-12.94,-4.36,-4.80,-5.40,NaN,0,191.80,0,-1.75,-10.31,10.57,-5.15];
    
    ind2=find(strcmp(staname2,sname{k}));
    lamda=lamda_list(k_list_man(ind2));
    amp=k_list_amp(ind2);
    expfit=amp*exp(-temp1/lamda);
    
    pcor{k}=detrend(A{k}-expfit);
    pfit{k}=A{k}-pcor{k};
end
% assemble normalized matrix of P, attempting to model ramp
if size(pcor,2)==10
    proxies=load('../pressure_matfiles/JoanGomberg/SlopesDataCEOF.mat');
else
    proxies=load('../pressure_matfiles/JoanGomberg/ShelfDataCEOF.mat');
end

% faster method that iteratively decreases step size of search
loopcheck1=0;
minguess=-20;
maxguess=20;
guessstep=10;
guessrange=minguess:guessstep:maxguess;
while loopcheck1<=2
    fit_check=[];
    loopcheck1=loopcheck1+1;
    for mm=1:length(guessrange)
        amp_guess=guessrange(mm);
        ramp_guess=zeros(size(X));
        ramp_guess(ipt-halflength:ipt+halflength-1)=[linspace(0,amp_guess,eventlength)]';
        ramp_guess(ipt+halflength:end)=ramp_guess(ipt+halflength:end)+amp_guess;
        ptemp=pcor;
        ptemp{ind}=detrend(ptemp{ind}-ramp_guess);
        for kk=1:length(ptemp)-1
            jj=find(strcmp(proxies.staname,sname{kk}));
            dP(:,kk)=pcor{kk}/std(pcor{kk});
            dT(:,kk)=proxies.datatemp(:,jj)/std(proxies.datatemp(:,jj));
            dS(:,kk)=proxies.datassh(:,jj)/std(proxies.datassh(:,jj));
        end
        % calculate CEOFs from combined dataset
        [~,loadings,pcs,~]=ceof([dP,dT,dS],1,true); % only need 1st CEOF
        pred1=real(pcs(:,1) * conj(loadings(:,1)).');
        fit_check(mm)=std(ptemp{ind}-pred1(:,ind)*std(ptemp{ind}));
    end
    if loopcheck1~=3
        [~,tempind]=mink(fit_check,2);
        minguess=min(guessrange(tempind));
        maxguess=max(guessrange(tempind));
        guessstep=guessstep/10;
        guessrange=minguess:guessstep:maxguess;
    end
end

% use best fit CEOF calculation
[~,mm]=min(fit_check);
amp_guess=guessrange(mm);
ramp_guess=zeros(size(X));
ramp_guess(ipt-halflength:ipt+halflength-1)=[linspace(0,amp_guess,eventlength)]';
ramp_guess(ipt+halflength:end)=ramp_guess(ipt+halflength:end)+amp_guess;
ptemp=pcor;
ptemp{ind}=detrend(ptemp{ind}-ramp_guess);
for kk=1:length(ptemp)-1
    jj=find(strcmp(proxies.staname,sname{kk}));
    dP(:,kk)=pcor{kk}/std(pcor{kk});
    dT(:,kk)=proxies.datatemp(:,jj)/std(proxies.datatemp(:,jj));
    dS(:,kk)=proxies.datassh(:,jj)/std(proxies.datassh(:,jj));
end
% calculate CEOFs from combined dataset
[~,loadings,pcs,~]=ceof([dP,dT,dS],1,true); % only need 1st CEOF
pred1=real(pcs(:,1) * conj(loadings(:,1)).');

% apply CEOF corrections to P
pfit2=pfit;
pcor2=pcor;
for ll=1:length(pcor)-1
    pfit2{ll}=pfit{ll}+real(pred1(:,ll))*std(ptemp{ll});
    pcor2{ll}=pcor{ll}-real(pred1(:,ll))*std(ptemp{ll});
end
better_ceof=pcor2;
X2=pcor2{ind}; X2=X2-mean(X2);
% determine ramp amplitude for newly-corrected pressure
for kk=1:length(ipt)
    xseg=t([1:ipt(kk)-halflength,ipt(kk)+halflength+1:end]);
    yseg1=X2(1:ipt(kk)-halflength);
    yseg3=X2(ipt(kk)+halflength+1:end);
    
    % set up inversion
    G=[xseg,[ones(size(yseg1));zeros(size(yseg3))],...
        [zeros(size(yseg1));ones(size(yseg3))]];
    m=inv(G'*G)*G'*[yseg1;yseg3];
    mdl(:,kk)=G*m;
    
    rampamp(kk)=m(3)-m(2);
    
    % correct timeseries w/ determined slope
    linpred=0:m(1):(length(X2)-1)*m(1);
    Xcorr=X2-linpred';
end

if length(varargin)==3
    figure(1); clf; hold on
    plot(X2,'b','linewidth',1)
    plot(xseg,G*m,'b','linewidth',0.5)
%     plot(Xcorr,'r','linewidth',1)
    xline(ipt); xline(ipt+14*24)
    yline(m(3),'k--')
    yline(m(2),'k--')
    legend(num2str(round(rampamp(kk),1)),'location','best')
    set(gca,'xticklabels',[])
    set(gca,'fontsize',16)
    box on; grid minor
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 11 8.5];
    print(['../AACSE_figures/ramp_detectability/single_station_detections/detrend/method_comp/' ...
        varargin{2} '/' varargin{3} 'cm/' varargin{1} '_' num2str(ipt)],'-djpeg','-r100')
    
    figure(13); hold on
    plot(t,X2+6*15,'k','linewidth',1)
    plot(t(xseg),G*m+6*15,'k','linewidth',0.5)
    set(gca,'fontsize',16)
    box on; grid minor
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 11 8.5];
    
%     figure(2); clf
%     subplot(211); hold on
%     subplot(212); hold on
%     subplot(211)
%     h1=plot(pcor{ind}+20,'r','linewidth',1);
%     plot(ramp_guess+20,'r','linewidth',1)
%     h2=plot(detrend(pcor{ind}-ramp_guess),'b','linewidth',1);
%     plot(pred1(:,ind)*std(ptemp{ind}),'k','linewidth',0.5)
%     xline(ipt); xline(ipt+14*24)
%     set(gca,'xticklabels',[])
%     set(gca,'fontsize',16)
%     legend([h1,h2],'with ramp','ramp corrected')
%     box on; grid minor
%     subplot(212)
%     plot(fit_check,'k','linewidth',1)
%     set(gca,'xtick',0:50:450); set(gca,'xticklabels',-20:5:25)
%     ylabel('misfit (cm)')
%     set(gca,'fontsize',16)
%     box on; grid minor
%     fh=gcf;
%     fh.PaperUnits='inches';
%     fh.PaperPosition=[0 0 11 8.5];
%     
%     figure(4); clf
%     subplot(121); hold on
%     subplot(122); hold on
%     for i=1:10
%         subplot(121)
%         plot(ptemp{i}+20*i,'b','linewidth',1)
%         plot(pred1(:,i)*std(ptemp{i})+20*i,'r','linewidth',1)
%         subplot(122)
%         plot(pcor2{i}+20*i,'k','linewidth',1)
%     end
%     subplot(121); box on; grid minor
%     set(gca,'xticklabels',[]); set(gca,'fontsize',16)
%     title('ramp-corrected P & CEOF')
%     subplot(122); box on; grid minor
%     set(gca,'xticklabels',[]); set(gca,'fontsize',16)
%     title('with-ramp difference')
%     fh=gcf;
%     fh.PaperUnits='inches';
%     fh.PaperPosition=[0 0 11 8.5];
end

sigflag=getsig(Xcorr,dt,rampamp,ipt);

end

function [ipt,better_ceof,rampamp,sigflag]=findoffset_ceof2(X,A,sname,ind,dt,varargin)
X=X-mean(X); % remove any residual average
t=[1:length(X)]';
iref=find(strcmp(sname,'LA21')); % find reference station index
A(iref)=[]; % remove reference station
sname(iref)=[]; % remove reference station
eventlength=1; % in hours

% detect ramp by usual method
if length(dt)==1
    for jj=1:length(eventlength)
        for ii=24*30+1:24:length(t)-24*30 % don't allow detections in first and last month
            xseg1=t(1:ii);
            xseg3=t(ii+eventlength(jj):end);
            xseg=t([1:ii,ii+eventlength(jj):end]);
            yseg1=X(1:ii);
            yseg3=X(ii+eventlength(jj):end);
            
            % set up inversion
            G=[xseg,[ones(size(yseg1));zeros(size(yseg3))],...
                [zeros(size(yseg1));ones(size(yseg3))]];
            m=inv(G'*G)*G'*[yseg1;yseg3];
            
            % correct timeseries w/ determined slope
            yseg1_corr=yseg1-polyval([m(1),m(2)],xseg1);
            yseg3_corr=yseg3-polyval([m(1),m(3)],xseg3);
            rmslist(ii,jj)=std([yseg1_corr;yseg3_corr]);
        end
    end
    rmslist(rmslist==0)=1000;
    [~,ipt]=min(rmslist);
else
    ipt=dt(2);
    dt=dt(1);
end

% re-calculate CEOFs, with modeled deformation ramp, selecting that with best fit
% correct out sensor drift
% better time basis
temp1=(t-t(1))/max(t-t(1));
if length(sname)>10
    lamda_list=linspace(1/max(t-t(1)),180/max(t-t(1)),1000)*24;
else
    lamda_list=linspace(1/max(t-t(1)),720/max(t-t(1)),1000)*24;
end
eventlength=14*24; % in hours
halflength=eventlength/2;

for k=1:length(A)
    % extract exponential time constant (empirically determined)
    staname2={'LA21';'LA22';'LA23';'LA25';'LA26';'LA28';'LA30';'LA32';'LA33';'LA34';'LA39';...
        'LT01';'LT02';'LT03';'LT04';'LT06';'LT07';'LT09';'LT10';'LT11';'LT12';'LT13';'LT14';'LT16';'LT20'};
%     k_list_man=[1000,1000,157,26,1000,28,1000,786,1000,1000,1000,1000,0.00,...
%         559,33,373,1000,0.00,1000,1000,1000,861,1000,876,144];
%     k_list_amp=[27.67,27.43,-10.39,-2.71,0.00,-4.69,8.52,87.39,30.09,33.34,-91.77,-9.98,0.00,...
%         -15.31,-4.60,-3.53,-11.21,0.00,-9.86,200.60,-11.93,-3.83,-18.80,16.75,-3.94];
    k_list_man=[1,311,53,19,NaN,22,1,162,315,479,1000,...
        1,NaN,318,36,254,288,NaN,1,964,1,130,400,502,146];
    k_list_amp=[0,14.79,-10.49,-2.50,NaN,-4.69,0,72.56,36.51,73.69,-763.67,...
        0,NaN,-12.94,-4.36,-4.80,-5.40,NaN,0,191.80,0,-1.75,-10.31,10.57,-5.15];
    
    ind2=find(strcmp(staname2,sname{k}));
    lamda=lamda_list(k_list_man(ind2));
    amp=k_list_amp(ind2);
    expfit=amp*exp(-temp1/lamda);
    
    pcor{k}=detrend(A{k}-expfit);
    pfit{k}=A{k}-pcor{k};
end

% assemble normalized matrix of P, attempting to model ramp

% faster method that iteratively decreases step size of search
loopcheck1=0;
minguess=-20;
maxguess=20;
guessstep=10;
guessrange=minguess:guessstep:maxguess;
while loopcheck1<=2
    fit_check=[];
    loopcheck1=loopcheck1+1;
    for mm=1:length(guessrange)
        amp_guess=guessrange(mm);
        ramp_guess=zeros(size(X));
        ramp_guess(ipt-halflength:ipt+halflength-1)=[linspace(0,amp_guess,eventlength)]';
        ramp_guess(ipt+halflength:end)=ramp_guess(ipt+halflength:end)+amp_guess;
        ptemp=pcor;
        ptemp{ind}=detrend(ptemp{ind}-ramp_guess);
        for kk=1:length(ptemp)
            dP(:,kk)=ptemp{kk}/std(ptemp{kk});
        end
        % calculate CEOFs from combined dataset
        [~,loadings,pcs,~]=ceof(dP,1,true); % only need 1st CEOF
        pred1=real(pcs(:,1) * conj(loadings(:,1)).');
        fit_check(mm)=std(ptemp{ind}-pred1(:,ind)*std(ptemp{ind}));
    end
    if loopcheck1~=3
        [~,tempind]=mink(fit_check,2);
        minguess=min(guessrange(tempind));
        maxguess=max(guessrange(tempind));
        guessstep=guessstep/10;
        guessrange=minguess:guessstep:maxguess;
    end
end

% use best fit CEOF calculation
[~,mm]=min(fit_check);
amp_guess=guessrange(mm);
ramp_guess=zeros(size(X));
ramp_guess(ipt-halflength:ipt+halflength-1)=[linspace(0,amp_guess,eventlength)]';
ramp_guess(ipt+halflength:end)=ramp_guess(ipt+halflength:end)+amp_guess;
ptemp=pcor;
ptemp{ind}=detrend(ptemp{ind}-ramp_guess);
for kk=1:length(ptemp)
    dP(:,kk)=ptemp{kk}/std(ptemp{kk});
end
% calculate CEOFs from combined dataset
[~,loadings,pcs,~]=ceof(dP,1,true); % only need 1st CEOF
pred1=real(pcs(:,1) * conj(loadings(:,1)).');

% apply CEOF corrections to P
for ll=1:length(pcor)
    pfit2{ll}=pfit{ll}+real(pred1(:,ll))*std(ptemp{ll});
    pcor2{ll}=pcor{ll}-real(pred1(:,ll))*std(ptemp{ll});
end
better_ceof=pcor2;
X2=pcor2{ind}; X2=X2-mean(X2);
% determine ramp amplitude for newly-corrected pressure
for kk=1:length(ipt)
    xseg=t([1:ipt(kk)-halflength,ipt(kk)+halflength+1:end]);
    yseg1=X2(1:ipt(kk)-halflength);
    yseg3=X2(ipt(kk)+halflength+1:end);
    
    % set up inversion
    G=[xseg,[ones(size(yseg1));zeros(size(yseg3))],...
        [zeros(size(yseg1));ones(size(yseg3))]];
    m=inv(G'*G)*G'*[yseg1;yseg3];
    mdl(:,kk)=G*m;
    
    rampamp(kk)=m(3)-m(2);
    
    % correct timeseries w/ determined slope
    linpred=0:m(1):(length(X2)-1)*m(1);
    Xcorr=X2-linpred';
end

if length(varargin)==3
    figure(13); hold on
    plot(t,X2+6*15,'k','linewidth',1)
    plot(t(xseg),G*m+6*15,'k','linewidth',0.5)
    set(gca,'fontsize',16)
    box on; grid minor
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 11 8.5];
end

if false
    figure(13); clf; hold on
    plot(t,X2,'k','linewidth',1)
    plot(t(xseg),G*m,'k')
end

sigflag=getsig(Xcorr,dt,rampamp,ipt);

end

% polynominal-based ramp detection
function [ipt,mdl,rampamp,sigflag]=findoffset_poly(X,dt,varargin)
X=X-mean(X); % remove any residual average
t=[1:length(X)]';
% detect on assumption of single time step ramp duration
if length(dt)==1
    for ii=24*30+1:24:length(t)-24*30 % don't allow detections in first and last month
        xseg1=t(1:ii);
        xseg3=t(ii+1:end);
        xseg=t([1:ii,ii+1:end]);
        yseg1=X(1:ii);
        yseg3=X(ii+1:end);
        
        % set up inversion
        G=[xseg,[ones(size(yseg1));zeros(size(yseg3))],...
            [zeros(size(yseg1));ones(size(yseg3))]];
        m=inv(G'*G)*G'*[yseg1;yseg3];
        
        % correct timeseries w/ determined slope
        yseg1_corr=yseg1-polyval([m(1),m(2)],xseg1);
        yseg3_corr=yseg3-polyval([m(1),m(3)],xseg3);
        rmslist(ii)=std([yseg1_corr;yseg3_corr]);
    end
    rmslist(rmslist==0)=1000;
    [~,ipt]=min(rmslist);
else
    ipt=dt(2);
    dt=dt(1);
end

% determine amplitude on assumption of [eventlength] ramp duration
eventlength=14*24; % in hours
halflength=eventlength/2;

for kk=1:length(ipt)
    %mask out assumed deformation signal
    x=(t-t(1))/max(t-t(1));
    xseg1=x(1:ipt(kk)-halflength);
    xseg2=x(ipt(kk)+halflength+1:end);
    yseg1=X(1:ipt(kk)-halflength);
    yseg2=X(ipt(kk)+halflength+1:end);
    
    %construct matrix for inversion
    G=[[xseg1.^3;xseg2.^3],[xseg1.^2;xseg2.^2],[xseg1;xseg2],...
        [ones(size(xseg1));zeros(size(xseg2))],[zeros(size(xseg1));ones(size(xseg2))]];
    m=inv(G'*G)*G'*([yseg1;yseg2]);
    mdl(:,kk)=m(1)*x.^3+m(2)*x.^2+m(3)*x; % exclude offsets from model to preserve ramp
    Xcorr=X-mdl(:,kk);
    rampamp(kk)=m(5)-m(4);
    
end

if ~isempty(varargin)
    figure(1); clf; hold on
    plot(t,X,'b','linewidth',1)
    plot(t([1:length(xseg1),length(xseg1)+eventlength:end]),G*m,'b','linewidth',0.5)
%     plot(Xcorr,'r','linewidth',1)
    xline(ipt); xline(ipt+14*24)
    yline(m(5),'k--')
    yline(m(4),'k--')
    legend(num2str(round(rampamp(kk),1)),'location','best')
    set(gca,'xticklabels',[])
    set(gca,'fontsize',16)
    box on; grid minor
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 11 8.5];
    print(['../AACSE_figures/ramp_detectability/single_station_detections/detrend/method_comp/' ...
        varargin{2} '/' varargin{3} 'cm/' varargin{1} '_' num2str(ipt)],'-djpeg','-r100')
    
    figure(11); hold on
    plot(t,X+1*15,'k','linewidth',1)
    plot(t([1:length(xseg1),length(xseg1)+eventlength:end]),G*m+1*15,'k','linewidth',0.5)
    set(gca,'fontsize',16)
    box on; grid minor
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 11 8.5];
end

sigflag=getsig(Xcorr,dt,rampamp,ipt);

end

% sinusoid-based ramp detection
function [ipt,mdl,rampamp,sigflag]=findoffset_sin(X,dt,varargin)
X=X-mean(X); % remove any residual average
t=[1:length(X)]';
% detect on assumption of single time step ramp duration
if length(dt)==1
    for ii=24*30+1:24:length(t)-24*30 % don't allow detections in first and last month
        xseg1=t(1:ii);
        xseg3=t(ii+1:end);
        xseg=t([1:ii,ii+1:end]);
        yseg1=X(1:ii);
        yseg3=X(ii+1:end);
        
        % set up inversion
        G=[xseg,[ones(size(yseg1));zeros(size(yseg3))],...
            [zeros(size(yseg1));ones(size(yseg3))]];
        m=inv(G'*G)*G'*[yseg1;yseg3];
        
        % correct timeseries w/ determined slope
        yseg1_corr=yseg1-polyval([m(1),m(2)],xseg1);
        yseg3_corr=yseg3-polyval([m(1),m(3)],xseg3);
        rmslist(ii)=std([yseg1_corr;yseg3_corr]);
    end
    rmslist(rmslist==0)=1000;
    [~,ipt]=min(rmslist);
else
    ipt=dt(2);
    dt=dt(1);
end

% determine amplitude on assumption of [eventlength] ramp duration
eventlength=14*24; % in hours
halflength=eventlength/2;

for kk=1:length(ipt)
    %mask out assumed deformation signal
    x=(t-t(1))/max(t-t(1));
    xseg1=x(1:ipt(kk)-halflength);
    xseg2=x(ipt(kk)+halflength+1:end);
    yseg1=X(1:ipt(kk)-halflength);
    yseg2=X(ipt(kk)+halflength+1:end);
    yf=365.25*24;
    
    %construct matrix for inversion
    G=[[sin(2*pi/yf*xseg1*max(t-t(1)));sin(2*pi/yf*xseg2*max(t-t(1)))],...
        [cos(2*pi/yf*xseg1*max(t-t(1)));cos(2*pi/yf*xseg2*max(t-t(1)))],...
        [xseg1;xseg2],[ones(size(xseg1));zeros(size(xseg2))],[zeros(size(xseg1));ones(size(xseg2))]];
    m=inv(G'*G)*G'*([yseg1;yseg2]);
    mdl(:,kk)=m(1)*sin(2*pi/yf*x*max(t-t(1)))+m(2)*cos(2*pi/yf*x*max(t-t(1)))+m(3)*x; % exclude offsets from model to preserve ramp
    Xcorr=X-mdl(:,kk);
    rampamp(kk)=m(5)-m(4);
end

if ~isempty(varargin)
    figure(2); clf; hold on
    plot(X,'b')
    plot([1:ipt,ipt+eventlength:length(X)],G*m,'b')
    plot(Xcorr,'r')
    xline(ipt,'k','linewidth',1)
    yline(m(5),'k--')
    yline(m(4),'k--')
    legend(num2str(round(rampamp(kk),1)),'location','best')
    box on
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 11 8.5];
    print(['../AACSE_figures/ramp_detectability/single_station_detections/detrend/method_comp/' ...
        varargin{2} '/' varargin{3} 'cm/' varargin{1} '_' num2str(ipt)],'-djpeg','-r100')
    
    figure(12); hold on
    plot(datenum(2018,06,08,01,00,00)+t/24,X+1*15,'k','linewidth',1)
    plot(datenum(2018,06,08,01,00,00)+t([1:length(xseg1),length(xseg1)+eventlength:end])/24,G*m+1*15,'k','linewidth',0.5)
    text(datenum(2018,06,08,01,00,00)+t(end)/24+5,median(X)+1*15,'Pssh','fontsize',16)
    set(gca,'fontsize',16)
    box on; grid minor
    datetick('x')
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 11 8.5];
end

if false
    figure(13); clf; hold on
    plot(t,X,'k','linewidth',1)
    plot([1:ipt,ipt+eventlength:length(X)],G*m,'k')
end

sigflag=getsig(Xcorr,dt,rampamp,ipt);

end

% modified student's t-test
function sigflag=getsig(X,dt,rampamp,ipt)

Neff1=getNeff(X(1:ipt),dt);
Neff2=getNeff(X(ipt+1:end),dt);
std2=std(X(ipt+1:end));
std1=std(X(1:ipt));
tstat=abs(rampamp)/sqrt(std2*std2/Neff2 + std1*std1/Neff1);
dof=Neff1+Neff2-2;
tcrit=tinv(.975,dof); %tinv is a 1-tail test, so for .05 probability in a 2-tail test, we need to test for .025 in the 1-tail function
%tcrit=0;
sigflag=0;
est=ipt*dt;
sigflag=tstat-tcrit;

if(tstat < tcrit)
  %Do not reject the null hypothesis
  sigflag=0;
end

end

% effective degrees of freedom
function [Neff]=getNeff(X,dt)
iplt=0;

[Ruu_simm, TM_simm] = xcorr(X-mean(X),'coeff'); %Auto-correlation value and lags
Ruu = Ruu_simm(TM_simm>=0); %Ruu = correlations with positive lags (since symmetric)
TM = TM_simm(TM_simm>=0); %TM = positive lags

if(iplt > 0)
npts=length(X);
smp=[0:1:npts-1];
figure(1); clf; hold on;
subplot(2,1,1); hold on;
plot(smp*dt,X,'-k');
title('Input Pressure');

subplot(2,1,2); hold on;
plot(TM*dt,Ruu,'-k');
title('Auto Correlation');
end

%Integrate
ipt=find(Ruu <= 0);
ie=length(Ruu);
if(length(ipt)> 0)
   ie=ipt(1);
end
if ie<2
    keyboard
end
TT_integral  = trapz(TM(1:ie),Ruu(1:ie));

Neff=length(X)/TT_integral;

end

function [pfit,pcor]=ssn_cor_smplexp(tin,pin,name,mdl,mdl_str)
%
% Applies polynomial, sinusoidal, or CEOF1 seasonal correction to input
% data.
%

% scaled time basis
temp1=(tin{1}-tin{1}(1))/max(tin{1}-tin{1}(1));

if strcmp(mdl_str,'poly') % exponential + 3rd degree polynomial
    for i=1:length(pin)
        % correct sensor drift
        if length(mdl{i})==4
            expfit=0;
        elseif length(mdl{i})==6
            expfit=mdl{i}(end-1)*exp(-temp1/mdl{i}(end));
        else
            warning('Unexpected exp fit model. Check pmod.')
            keyboard
        end

        %construct matrix for inversion
        Gp=[temp1.^3,temp1.^2,temp1,ones(size(temp1))];
        mp=inv(Gp'*Gp)*Gp'*(pin{i}-expfit);
        pfit{i}=Gp*mp;
        pcor{i}=pin{i}-pfit{i}-expfit;
    end
elseif strcmp(mdl_str,'sin') % exponential + sinusoid
    for l=1:length(pin)
        % correct sensor drift
        if length(mdl{l})==4
            expfit=0;
        elseif length(mdl{l})==6
            expfit=mdl{l}(end-1)*exp(-temp1/mdl{l}(end));
        else
            warning('Unexpected exp fit model. Check smod.')
            keyboard
        end

        %construct matrix for inversion
        Gs=[sin(2*pi/yf*temp1*max(tin{l}-tin{l}(1))),cos(2*pi/yf*temp1*max(tin{l}-tin{l}(1))),ones(size(temp1)),temp1];
        ms=inv(Gs'*Gs)*Gs'*(pin{l}-expfit);
        pfit{l}=Gs*ms;
        pcor{l}=pin{l}-pfit{l}-expfit;
    end
elseif strcmp(mdl_str,'lin')
    for k=1:length(pin)
        % correct sensor drift
        if length(mdl{k})==2
            expfit=0;
        elseif length(mdl{k})==4
            expfit=mdl{k}(end-1)*exp(-temp1/mdl{k}(end));
        else
            warning('Unexpected exp fit model. Check lmod.')
            keyboard
        end

        %construct matrix for inversion
        Gl=[temp1,ones(size(temp1))];
        ml=inv(Gl'*Gl)*Gl'*(pin{k}-expfit);
        pfit{k}=Gl*ml;
        pcor{k}=pin{k}-pfit{k}-expfit;
    end
elseif strcmp(mdl_str,'ceof')
    for j=1:length(pin)
        % correct sensor drift
        if length(mdl{j})==0
            expfit=0;
        elseif length(mdl{j})==4
            expfit=mdl{j}(end-1)*exp(-temp1/mdl{j}(end));
        else
            warning('Unexpected exp fit model. Check lmod.')
            keyboard
        end

        %construct matrix for inversion
        Gc=[temp1,ones(size(temp1))];
        mc=inv(Gc'*Gc)*Gc'*(pin{j}-expfit);
        pfit{j}=Gl*ml;
        pcor_temp{j}=pin{j}-pfit{j}-expfit;
    end

    % unify the time basis
    t_min=max(cellfun(@(v)v(1),tin));
    t_max=min(cellfun(@(v)v(end),tin));
    
    % once working, I can code in logic for maximizing temporal extent by
    % excluding truncated stations
    for kk=1:length(pcor_temp)
        cond=(tin{kk}>=t_min+0.01) & (tin{kk}<=t_max+0.01);
        ptemp{kk}=pcor_temp{kk}(cond)'; ptemp{kk}=detrend(ptemp{kk});
        dP(:,kk)=ptemp{kk}/std(ptemp{kk});
    end
    % calculate CEOFs from combined dataset
    [~,loadings,pcs,~]=ceof(dP,1,true); % only need 1st CEOF
    pred1=real(pcs(:,1) * conj(loadings(:,1)).');
    scaled_pred1=mat2cell(pred1.*cellfun(@std,ptemp),length(dP),ones(1,length(pcor_temp)));
    pcor=cellfun(@minus,ptemp,scaled_pred1,'UniformOutput',false);
end
end

function make_ceof_plots(tin,dP,dT,dS,pcor,pred1)
for i=1:length(pcor)-1
    plot(tin{1},dP(:,i)*std(pcor{i})+20*(i-1),'b','linewidth',1)
    plot(tin{1},real(pred1(:,i))*std(pcor{i})+20*(i-1),'r','linewidth',1)
end
datetick('x')
xlim([datenum(2018,6,1) datenum(2019,7,1)])
datetick('x',3,'keeplimits')
xtickangle(45)
set(gca,'fontsize',14)
ylabel('P (cm)')
ylim([-20 240])
box on
grid minor

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print('../AACSE_figures/pressure_stacks/ceof/AGU/Pceof_slope','-dtiff','-r300')

figure(9); clf; hold on
for i=1:length(pcor)-1
    plot(tin{1},dT(:,i)*std(pcor{i})+20*(i-1),'k','linewidth',1)
    plot(tin{1},real(pred1(:,i+length(pcor)))*std(pcor{i})+20*(i-1),'r','linewidth',1)
end
datetick('x')
xlim([datenum(2018,6,1) datenum(2019,7,1)])
datetick('x',3,'keeplimits')
xtickangle(45)
set(gca,'fontsize',14)
ylabel('T (cm)')
ylim([-20 240])
box on
grid minor

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print('../AACSE_figures/pressure_stacks/ceof/AGU/Tceof_slope','-dtiff','-r300')

figure(10); clf; hold on
for i=1:length(pcor)-1
    plot(tin{1},dS(:,i)*std(pcor{i})+20*(i-1),'k','linewidth',1)
    plot(tin{1},real(pred1(:,i+2*length(pcor)))*std(pcor{i})+20*(i-1),'r','linewidth',1)
end
datetick('x')
xlim([datenum(2018,6,1) datenum(2019,7,1)])
datetick('x',3,'keeplimits')
xtickangle(45)
set(gca,'fontsize',14)
ylabel('SSH (cm)')
ylim([-20 240])
box on
grid minor

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print('../AACSE_figures/pressure_stacks/ceof/AGU/SSHceof_slope','-dtiff','-r300')
end