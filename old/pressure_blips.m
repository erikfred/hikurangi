% pressure_blips.m
%
% Import troublesome APG data from any year WITHOUT correcting
% blips/offsets. Decimate to [15 minute] and [15 second] sample interval to
% mimic what's done on the DARTs.
%

%----- GNS21-PC
fid = '/Volumes/TOSHIBA_EXT/apg_data/2021_2022_Gisborne/GNS21-PC.csv';
P = readtable(fid);
P = [datenum(table2array(P(:,1))),table2array(P(:,2:3))]';

% trim ends
P = P(:,16000:30019000);

% decimation loop
GNS21PC.t_m=[]; GNS21PC.T_m=[]; GNS21PC.p_m=[];
GNS21PC.t_t=[]; GNS21PC.T_t=[]; GNS21PC.p_t=[];
i1 = 1;
d2 = floor(P(1,1)) + 1;
while i1<length(P)
    i2 = find(P(1,:)>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segC,~,~] = downsample_uneven(P(1,i1:i2-1),P(2:3,i1:i2-1),1/24/4);
    GNS21PC.t_m=cat(2,GNS21PC.t_m,segt);
    GNS21PC.p_m=cat(2,GNS21PC.p_m,segC(1,:));
    GNS21PC.T_m=cat(2,GNS21PC.T_m,segC(2,:));
    [segt,segC,~,~] = downsample_uneven(P(1,i1:i2-1),P(2:3,i1:i2-1),1/24/60/4);
    GNS21PC.t_t=cat(2,GNS21PC.t_t,segt);
    GNS21PC.p_t=cat(2,GNS21PC.p_t,segC(1,:));
    GNS21PC.T_t=cat(2,GNS21PC.T_t,segC(2,:));
    i1=i2;
    d2=floor(P(1,i2))+1;
end

% timestamps of anomalies (empirical)
GNS21PC.t_anom=P(1,8864934);

% plot
figure(44); clf; hold on
plot(GNS21PC.t_t,GNS21PC.p_t,'.-')
plot(GNS21PC.t_m,GNS21PC.p_m,'o-')
xline(GNS21PC.t_anom)

clearvars('P')

%----- GNS21-PD
fid = '/Volumes/TOSHIBA_EXT/apg_data/2021_2022_Gisborne/GNS21-PD.csv';
P = readtable(fid);
P = [datenum(table2array(P(:,1))),table2array(P(:,2:3))]';

% trim ends
P = P(:,11000:30099450);

% decimation loop
GNS21PD.t_m=[]; GNS21PD.T_m=[]; GNS21PD.p_m=[];
GNS21PD.t_t=[]; GNS21PD.T_t=[]; GNS21PD.p_t=[];
i1 = 1;
d2 = floor(P(1,1)) + 1;
while i1<length(P)
    i2 = find(P(1,:)>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segD,~,~] = downsample_uneven(P(1,i1:i2-1),P(2:3,i1:i2-1),1/24/4);
    GNS21PD.t_m=cat(2,GNS21PD.t_m,segt);
    GNS21PD.p_m=cat(2,GNS21PD.p_m,segD(1,:));
    GNS21PD.T_m=cat(2,GNS21PD.T_m,segD(2,:));
    [segt,segD,~,~] = downsample_uneven(P(1,i1:i2-1),P(2:3,i1:i2-1),1/24/60/4);
    GNS21PD.t_t=cat(2,GNS21PD.t_t,segt);
    GNS21PD.p_t=cat(2,GNS21PD.p_t,segD(1,:));
    GNS21PD.T_t=cat(2,GNS21PD.T_t,segD(2,:));
    i1=i2;
    d2=floor(P(1,i2))+1;
end

% timestamps of anomalies (empirical)
GNS21PD.t_anom=[P(1,4363500),P(1,4430000),...
    P(1,4453300),P(1,4542500)];

% plot
figure(44); clf; hold on
plot(GNS21PD.t_t,GNS21PD.p_t,'.-')
plot(GNS21PD.t_m,GNS21PD.p_m,'o-')
for i=1:length(GNS21PD.t_anom)
    xline(GNS21PD.t_anom(i))
end

clearvars('P')

%----- GNS21-PI
fid = '/Volumes/TOSHIBA_EXT/apg_data/2021_2022_Gisborne/GNS21-PI.csv';
P = readtable(fid);
P = [datenum(table2array(P(:,1))),table2array(P(:,2:3))]';

% trim ends
P = P(:,14718:32877925);

% decimation loop
GNS21PI.t_m=[]; GNS21PI.T_m=[]; GNS21PI.p_m=[];
GNS21PI.t_t=[]; GNS21PI.T_t=[]; GNS21PI.p_t=[];
i1 = 1;
d2 = floor(P(1,1)) + 1;
while i1<length(P)
    i2 = find(P(1,:)>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segD,~,~] = downsample_uneven(P(1,i1:i2-1),P(2:3,i1:i2-1),1/24/4);
    GNS21PI.t_m=cat(2,GNS21PI.t_m,segt);
    GNS21PI.p_m=cat(2,GNS21PI.p_m,segD(1,:));
    GNS21PI.T_m=cat(2,GNS21PI.T_m,segD(2,:));
    [segt,segD,~,~] = downsample_uneven(P(1,i1:i2-1),P(2:3,i1:i2-1),1/24/60/4);
    GNS21PI.t_t=cat(2,GNS21PI.t_t,segt);
    GNS21PI.p_t=cat(2,GNS21PI.p_t,segD(1,:));
    GNS21PI.T_t=cat(2,GNS21PI.T_t,segD(2,:));
    i1=i2;
    d2=floor(P(1,i2))+1;
end

% timestamps of anomalies (empirical)
GNS21PI.t_anom=[P(1,27003466)];

% plot
figure(44); clf; hold on
plot(GNS21PI.t_t,GNS21PI.p_t,'.-')
plot(GNS21PI.t_m,GNS21PI.p_m,'o-')
for i=1:length(GNS21PI.t_anom)
    xline(GNS21PI.t_anom(i))
end

clearvars('P')

%----- GNS22-PL
fid = '../apg_data/GONDOR-I/GNS22-PL.csv';
P = readtable(fid); % [t P T]
P = [datenum(table2array(P(:,1))),table2array(P(:,2:3))]';

% trim ends
P = P(:,39000:30720000);

% decimation loop
GNS22PL.t_m=[]; GNS22PL.T_m=[]; GNS22PL.p_m=[];
GNS22PL.t_t=[]; GNS22PL.T_t=[]; GNS22PL.p_t=[];
i1 = 1;
d2 = floor(P(1,1)) + 1;
while i1<length(P)
    i2 = find(P(1,:)>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segD,~,~] = downsample_uneven(P(1,i1:i2-1),P(2:3,i1:i2-1),1/24/4);
    GNS22PL.t_m=cat(2,GNS22PL.t_m,segt);
    GNS22PL.p_m=cat(2,GNS22PL.p_m,segD(1,:));
    GNS22PL.T_m=cat(2,GNS22PL.T_m,segD(2,:));
    [segt,segD,~,~] = downsample_uneven(P(1,i1:i2-1),P(2:3,i1:i2-1),1/24/60/4);
    GNS22PL.t_t=cat(2,GNS22PL.t_t,segt);
    GNS22PL.p_t=cat(2,GNS22PL.p_t,segD(1,:));
    GNS22PL.T_t=cat(2,GNS22PL.T_t,segD(2,:));
    i1=i2;
    d2=floor(P(1,i2))+1;
end

% timestamps of anomalies (empirical)
GNS22PL.t_anom=[P(1,4363500),P(1,4430000),...
    P(1,4453300),P(1,4542500)];

% plot
figure(44); clf; hold on
plot(GNS22PL.t_t,GNS22PL.p_t,'.-')
plot(GNS22PL.t_m,GNS22PL.p_m,'o-')
for i=1:length(GNS22PL.t_anom)
    xline(GNS22PL.t_anom(i))
end

clearvars('P')