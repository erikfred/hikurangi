% CORK_formation_ResidualTides.m
%
% Taking a look at your residual signals
%

clear; close all

%% same stuff as before, just loading up the data
staname={'U1518','U1519'};
stalat=[-38.859,-38.727];
stalon=[178.896,178.615];
stadepth=[0,0];

%---U1518 Part 1---%
fid = fopen('/Volumes/TOSHIBA_EXT/apg_data/CORK_LTBPR/U1518/U1518_RR1902_clean.dat','r');
fspec = '%s %s %f %f %f %f %f %f %f %f %f'; % [t t T T1 P1 T2 P2 T3 P3 T4 P4]
E = textscan(fid,fspec,'HeaderLines',3);
fclose(fid);
t_str=cat(2,cat(1,E{1}{:}),repmat(' ',length(E{1}),1),cat(1,E{2}{:}));
tE = datetime(t_str,'InputFormat','yyyy-MM-dd HH:mm:ss');
E = [E{5},E{4},E{7},E{6},E{9},E{8},E{11},E{10}]; % [P1 T1 P2 T2 P3 T3 P4 T4]

%---U1518 Part 2---%
fid = fopen('/Volumes/TOSHIBA_EXT/apg_data/CORK_LTBPR/U1518/U1518_TAN2102.dat','r');
fspec = '%s %s %f %f %f %f %f %f %f %f %f'; % [t t T T1 P1 T2 P2 T3 P3 T4 P4]
F = textscan(fid,fspec);
fclose(fid);
t_str=cat(2,cat(1,F{1}{:}),repmat(' ',length(F{1}),1),cat(1,F{2}{:}));
tF = datetime(t_str,'InputFormat','yyyy-MM-dd HH:mm:ss');
F{11} = [F{11}(1:525835);F{11}(525836:end)-F{11}(525836)+F{11}(525835)+...
    mean(diff(F{11}(525820:525835)))]; % corrects an artificial offset
tE = [tE;tF];
E = cat(1,E,[F{5},F{4},F{7},F{6},F{9},F{8},F{11},F{10}]); % [P1 T1 P2 T2 P3 T3 P4 T4]

%---U1518 Part 3---%
fid = fopen('/Volumes/TOSHIBA_EXT/apg_data/CORK_LTBPR/U1518/U1518_2023_formatted.dat','r');
fspec = '%s %s %f %f %f %f %f %f %f %f %f'; % [t t T T1 P1 T2 P2 T3 P3 T4 P4]
K = textscan(fid,fspec);
fclose(fid);
t_str=cat(2,cat(1,K{1}{:}),repmat(' ',length(K{1}),1),cat(1,K{2}{:}));
tK = datetime(t_str,'InputFormat','yyyy-MM-dd HH:mm:ss');
tE = [tE;tK];
E = cat(1,E,[K{5}/10,K{4},K{7}/10,K{6},K{9}/10,K{8},K{11}/10+0.15,K{10}]); % [P1 T1 P2 T2 P3 T3 P4 T4]

tE = tE(3123:end);
E = E(3123:end,:);
E(92703:92706,7)=linspace(E(92703,7),E(92706,7),4);

%% okay, now let's take a look

% sub-selection from early on (June 2018)
i1=100000; i2=110000;
p1 = E(i1:i2,1); % this is the deepest formation pressure
p2 = E(i1:i2,7); % this is the wellhead pressure
t = datenum(tE(i1:i2)); % this converts time to a number that we can manipulate

% plot the pressure selections
figure(11); clf; hold on
plot(t,p1,'o-')
yyaxis right
plot(t,p2,'s-')
yyaxis left
legend('formation pressure','wellhead pressure')
datetick('x')
ylabel('Pressure (dbar)')
set(gca,'fontsize',14)
box on; grid on

% look at their cross-correlation
[c,lags]=xcorr(detrend(p1),detrend(p2),'unbiased');
figure(12); clf; hold on
plot(lags,c,'o-')
xlabel('Lag (min)')
ylabel('Correlation Coefficient')
set(gca,'fontsize',14)
xlim([-100 100])
box on; grid on

% do the inversion to model wellhead to formation
tinv = (t-t(1))/max(t-t(1)); % it is computationally easier to use a low-order number system
G = [ones(size(t)),tinv,p2];
m=inv(G'*G)*G'*p1;
p1_model = G*m;
p1_corrected = p1-m(3)*p2;

figure(13); clf; hold on
plot(t,p1-p1_model+mean(p1)-mean(p1-p1_model+mean(p1)),'k-')
legend('residuals')
ylabel('Pressure (dbar)')
datetick('x')
set(gca,'fontsize',14)
box on; grid on

disp(['Misfit RMS = ' num2str(rms(p1-p1_model)) ' dbar'])

% now lets look at a later segment (January 2022)
i1=2000000; i2=2010000;
p1 = E(i1:i2,1); % this is the deepest formation pressure
p2 = E(i1:i2,7); % this is the wellhead pressure
t = datenum(tE(i1:i2)); % this converts time to a number that we can manipulate

% plot the pressure selections
figure(21); clf; hold on
plot(t,p1,'o-')
yyaxis right
plot(t,p2,'s-')
yyaxis left
legend('formation pressure','wellhead pressure')
datetick('x')
ylabel('Pressure (dbar)')
set(gca,'fontsize',14)
box on; grid on

% look at their cross-correlation
[c,lags]=xcorr(detrend(p1),detrend(p2),'unbiased');
figure(22); clf; hold on
plot(lags,c,'o-')
xlabel('Lag (min)')
ylabel('Correlation Coefficient')
set(gca,'fontsize',14)
xlim([-100 100])
box on; grid on

% do the inversion to model wellhead to formation
tinv = (t-t(1))/max(t-t(1)); % it is computationally easier to use a low-order number system
G = [ones(size(t)),tinv,p2];
m=inv(G'*G)*G'*p1;
p1_model = G*m;
p1_corrected = p1-m(3)*p2;

figure(23); clf; hold on
plot(t,p1-p1_model+mean(p1)-mean(p1-p1_model+mean(p1)),'k-')
legend('residuals')
ylabel('Pressure (dbar)')
datetick('x')
set(gca,'fontsize',14)
box on; grid on

disp(['Misfit RMS = ' num2str(rms(p1-p1_model)) ' dbar'])

%-----------------------%

% now we'll just do the above for the entire dataset

tinv=datenum(tE-tE(1)); tinv=tinv/max(tinv);
lamda=1:0.1:20;
for i=1:length(lamda)
    Gl=[ones(size(tinv(30000:end))),tinv(30000:end),E(30000:end,7),exp(-lamda(i)*tinv(30000:end))];
    ml=inv(Gl'*Gl)*Gl'*E(30000:end,1);
    fit_test(i)=rms(E(30000:end,1)-Gl*ml);
end
[~,imin]=min(fit_test);
Gl=[ones(size(tinv(30000:end))),tinv(30000:end),E(30000:end,7),exp(-lamda(imin)*tinv(30000:end))];
ml=inv(Gl'*Gl)*Gl'*E(30000:end,1);
E_cor(:,1)=E(30000:end,1)-Gl*ml;

% But if we assume that our model based on the subset is a good one for the
% entire dataset (minus the exponential part), then we can still apply the
% model correction to the entire thing:
E_cor = E(:,1);
E_cor(:,1) = E(:,1)-(ml(1)+tinv*ml(2)+ml(3)*E(:,7));

% We need to calculate another model (of the same form) for the other
% formation pressures
ml2=inv(Gl'*Gl)*Gl'*E(1000000:2500000,3);
E_cor(:,2) = E(:,3)-(ml2(1)+tinv*ml2(2)+ml2(3)*E(:,7));

ml3=inv(Gl'*Gl)*Gl'*E(1000000:2500000,5);
E_cor(:,3) = E(:,5)-(ml3(1)+tinv*ml3(2)+ml3(3)*E(:,7));

figure(13); clf; hold on
plot(tE,E_cor(:,1),'linewidth',1)
plot(tE,E_cor(:,2),'linewidth',1)
plot(tE,E_cor(:,3),'linewidth',1)
legend('formation 1','formation 2','formation 3')
ylabel('Pressure (dbar)')
datetick('x')
set(gca,'fontsize',14)
box on; grid on

fit((1:260)',test,'a*exp(b*x)+c*x+d','start',[8 -0.04 -0.01 0])

%% final thoughts

% You can see in the plots of the corrected pressure that there are still a
% number of features that could be cleaned up, like the exponential and the
% blips where the different data files were stitched together. I'd be
% happy to chat about these things down the road, but I think this is
% plenty for you to go off of for now.

% Use the above code as a guide to do the same import and processing for
% the 1519 data and reach out if you hit any snags.