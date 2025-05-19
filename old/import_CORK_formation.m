% import_CORK_formation.m
%
% Quick and dirty processing of formation pressures for reference
%

clear; close all

% station info
staname={'U1518','U1519'};
stalat=[-38.859,-38.727]; % BPR location approximated from map view
stalon=[178.896,178.615]; % BPR location approximated from map view
stadepth=[0,0]; % what are the sensor depths?

%---U1518 Part 1---%
% read in data
% fid = fopen('../apg_data/CORK_LTBPR/U1518/U1518_RR1902_clean.dat','r');
fid = fopen('/Volumes/Gorgoroth/apg_data/CORK_LTBPR/U1518/U1518_RR1902_clean.dat','r');
fspec = '%s %s %f %f %f %f %f %f %f %f %f'; % [t t T T1 P1 T2 P2 T3 P3 T4 P4]
E = textscan(fid,fspec,'HeaderLines',3);
fclose(fid);
t_str=cat(2,cat(1,E{1}{:}),repmat(' ',length(E{1}),1),cat(1,E{2}{:}));
tE = datenum(t_str,'yyyy-mm-dd HH:MM:SS');
E = [tE,E{5},E{4},E{7},E{6},E{9},E{8},E{11},E{10}]'; % [t P1 T1 P2 T2 P3 T3 P4 T4]

%---U1518 Part 2---%
% read in data
% fid = fopen('../apg_data/CORK_LTBPR/U1518/U1518_TAN2102.dat','r');
fid = fopen('/Volumes/Gorgoroth/apg_data/CORK_LTBPR/U1518/U1518_TAN2102.dat','r');
fspec = '%s %s %f %f %f %f %f %f %f %f %f'; % [t t T T1 P1 T2 P2 T3 P3 T4 P4]
F = textscan(fid,fspec);
fclose(fid);
t_str=cat(2,cat(1,F{1}{:}),repmat(' ',length(F{1}),1),cat(1,F{2}{:}));
tF = datenum(t_str,'yyyy-mm-dd HH:MM:SS');
F{11} = [F{11}(1:525835);F{11}(525836:end)-F{11}(525836)+F{11}(525835)+...
    mean(diff(F{11}(525820:525835)))]; % corrects an artificial offset
E = cat(2,E,[tF,F{5},F{4},F{7},F{6},F{9},F{8},F{11},F{10}]'); % [t P T]

%---U1518 Part 3---%
% read in data
% fid = fopen('../apg_data/CORK_LTBPR/U1518_2023_formatted.dat','r');
fid = fopen('/Volumes/Gorgoroth/apg_data/CORK_LTBPR/U1518/U1518_2023_formatted.dat','r');
fspec = '%s %s %f %f %f %f %f %f %f %f %f'; % [t t T T1 P1 T2 P2 T3 P3 T4 P4]
K = textscan(fid,fspec);
fclose(fid);
t_str=cat(2,cat(1,K{1}{:}),repmat(' ',length(K{1}),1),cat(1,K{2}{:}));
tK = datenum(t_str,'yyyy-mm-dd HH:MM:SS');
E = cat(2,E,[tK,K{5}/10,K{4},K{7}/10,K{6},K{9}/10,K{8},K{11}/10+0.15,K{10}]'); % [t P T]

% empirical trimming
E = E(:,3123:end);
E(2,92703:92706)=linspace(E(2,92703),E(2,92706),4);

% use wellhead pressure to remove oceanographic effects
tinv=E(1,:)-E(1,1); tinv=tinv/max(tinv);
Gl=[ones(size(tinv(1000000:2500000)))',tinv(1000000:2500000)',E(8,1000000:2500000)'];
ml=inv(Gl'*Gl)*Gl'*E(2,1000000:2500000)';

E_cor = E(1,:);
E_cor(2,:) = E(2,:)-(ml(1)+tinv*ml(2)+ml(3)*E(8,:));

ml=inv(Gl'*Gl)*Gl'*E(4,1000000:2500000)';
E_cor(3,:) = E(4,:)-(ml(1)+tinv*ml(2)+ml(3)*E(8,:));

ml=inv(Gl'*Gl)*Gl'*E(6,1000000:2500000)';
E_cor(4,:) = E(6,:)-(ml(1)+tinv*ml(2)+ml(3)*E(8,:));

% wellhead pressures have their own oddities that need trimming
E_cor(2:end,469337:469383)=NaN;
E_cor(2:end,1562830:1562890)=NaN;
E_cor(2:end,2628028:end)=NaN;

% decimation loop
te=[];
Te=[];
e=[];
e1=[];
e2=[];
e3=[];
i1 = 1;
d2 = floor(E(1,1))+1;
while i1<length(E)
    i2 = find(E(1,:)>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segE,~,~] = downsample_uneven(E(1,i1:i2-1),[E(8:9,i1:i2-1);E_cor(2:4,i1:i2-1)],1/24);
    if length(segt)>24
        keyboard
    end
    te=cat(2,te,segt);
    e=cat(2,e,segE(1,:));
    e1=cat(2,e1,segE(3,:));
    e2=cat(2,e2,segE(4,:));
    e3=cat(2,e3,segE(5,:));
    Te=cat(2,Te,segE(2,:));
    i1=i2;
    d2=floor(E(1,i2))+1;
end

% scale to hPa, tidal filter
e=e*100; % [hPa]
e1=e1*100; e2=e2*100; e3=e3*100;
ef=Z_godin(e);
Tef=Z_godin(Te);

% remove NaNs from tidal filter
te(isnan(ef))=[];
Te(isnan(ef))=[];
Tef(isnan(ef))=[];
e(isnan(ef))=[];
e1(isnan(ef))=[];
e2(isnan(ef))=[];
e3(isnan(ef))=[];
ef(isnan(ef))=[];

% interpolate onto monotonic time basis on the hour
tef = te(1)+datenum(0,0,0,0,30,0):1/24:te(end)-datenum(0,0,0,0,30,0);
ef = interp1(te,ef,tef,'spline'); % spline helps with gap at March 2021
e1 = interp1(te,e1,tef,'spline');
e2 = interp1(te,e2,tef,'spline');
e3 = interp1(te,e3,tef,'spline');
Tef = interp1(te,Tef,tef);

clearvars('t_str','tC','C','tD','D','tE','E','tF','F','tK','K')

%---U1519 Part 1---%
% read in data
fid = fopen('../apg_data/CORK_LTBPR/U1519/U1519_RR1902_clean.dat','r');
% fid = fopen('/Volumes/Gorgoroth/apg_data/CORK_LTBPR/U1519/U1519_RR1902_clean.dat','r');
fspec = '%s %s %f %f %f %f %f %f %f'; % [t t ? T1 P1 T2 P2 T3 P3]
G = textscan(fid,fspec,'HeaderLines',3);
fclose(fid);
t_str=cat(2,cat(1,G{1}{:}),repmat(' ',length(G{1}),1),cat(1,G{2}{:}));
tG = datenum(t_str,'yyyy-mm-dd HH:MM:SS');
G = [tG,G{5},G{4},G{7},G{6},G{9},G{8}]'; % [t P1 T1 P2 T2 P T]

%---U1519 Part 2---%
% read in data
fid = fopen('../apg_data/CORK_LTBPR/U1519/U1519_TAN2102.dat','r');
% fid = fopen('/Volumes/Gorgoroth/apg_data/CORK_LTBPR/U1519/U1519_TAN2102.dat','r');
fspec = '%s %s %f %f %f %f %f %f %f'; % [t t ? T1 P1 T2 P2 T3 P3]
H = textscan(fid,fspec);
fclose(fid);
t_str=cat(2,cat(1,H{1}{:}),repmat(' ',length(H{1}),1),cat(1,H{2}{:}));
tH = datenum(t_str,'yyyy-mm-dd HH:MM:SS');
G =  cat(2,G,[tH,H{5},H{4},H{7},H{6},H{9},H{8}]'); % [t P1 T1 P2 T2 P T]

%---U1519 Part 3---%
% read in data
fid = fopen('../apg_data/CORK_LTBPR/U1519/U1519_TN415_2023.dat','r');
% fid = fopen('/Volumes/Gorgoroth/apg_data/CORK_LTBPR/U1519/U1519_TN415_2023.dat','r');
fspec = '%f %f %f %f %f %f %f %f %f'; % [t ? T T1 P1 T2 P2 T3 P3]
J = textscan(fid,fspec,'HeaderLines',7);
fclose(fid);
% time is as yymmddHHMMSS
temp.base = J{1};
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
J{1} = datenum(2000+temp.yr,temp.mnth,temp.dy,temp.hr,temp.min,temp.sec);
G = cat(2,G,[J{1},J{5}/9.1415-99.62,J{4},J{7}/10,J{6},J{9}/10,J{8}]'); % [t P1 T1 P2 T2 P T]

% empirical trimming
G = G(:,2874:end);

% use wellhead pressure to remove oceanographic effects
tinv=G(1,:)-G(1,1); tinv=tinv/max(tinv);
Gl=[ones(size(tinv(1000000:2500000)))',tinv(1000000:2500000)',G(6,1000000:2500000)'];
ml=inv(Gl'*Gl)*Gl'*G(2,1000000:2500000)';

G_cor = G(1,:);
G_cor(2,:) = G(2,:)-(ml(1)+tinv*ml(2)+ml(3)*G(6,:));

ml=inv(Gl'*Gl)*Gl'*G(4,1000000:2500000)';
G_cor(3,:) = G(4,:)-(ml(1)+tinv*ml(2)+ml(3)*G(6,:));

% wellhead pressures have their own oddities that need trimming
G_cor(2:3,439372:439418)=NaN;
G_cor(2:3,1530180:1530270)=NaN;
G_cor(2:3,2598520:end)=NaN;

% decimation loop
tg=[];
Tg=[];
g=[];
g1=[];
g2=[];
i1 = 1;
d2 = floor(G(1,1))+1;
while i1<length(G)
    i2 = find(G(1,:)>=d2,1);
    if isempty(i2)
        break
    end
    [segt,segG,~,~] = downsample_uneven(G(1,i1:i2-1),[G(6:7,i1:i2-1);G_cor(2:3,i1:i2-1)],1/24);
    if length(segt)>24
        keyboard
    end
    tg=cat(2,tg,segt);
    g=cat(2,g,segG(1,:));
    g1=cat(2,g1,segG(3,:));
    g2=cat(2,g2,segG(4,:));
    Tg=cat(2,Tg,segG(2,:));
    i1=i2;
    d2=floor(G(1,i2))+1;
end

% scale to hPa, tidal filter
g=g*100; % [hPa]
g1=g1*100; g2=g2*100;
gf=Z_godin(g);
Tgf=Z_godin(Tg);

% remove NaNs from tidal filter
tg(isnan(gf))=[];
Tg(isnan(gf))=[];
Tgf(isnan(gf))=[];
g(isnan(gf))=[];
g1(isnan(gf))=[];
g2(isnan(gf))=[];
gf(isnan(gf))=[];

% interpolate onto monotonic time basis on the hour
tgf = tg(1)+datenum(0,0,0,0,30,0):1/24:tg(end)-datenum(0,0,0,0,30,0);
gf = interp1(tg,gf,tgf);
g1 = interp1(tg,g1,tgf,'spline');
g2 = interp1(tg,g2,tgf,'spline');
Tgf = interp1(tg,Tgf,tgf);

clearvars('t_str','tG','G','tH','H','temp','J')

%-----PLOTTING-----%
% make adjustments as needed for targeted figs
figure(2); clf; hold on
plot(tef(tef>=datenum(2020,10,01) & tef<=datenum(2021,10,1)),...
    detrend(e1(tef>=datenum(2020,10,01) & tef<=datenum(2021,10,1))),...
    'linewidth',2,'color',[0,0,128]/255)
plot(tef(tef>=datenum(2020,10,01) & tef<=datenum(2021,10,1)),...
    detrend(e2(tef>=datenum(2020,10,01) & tef<=datenum(2021,10,1)))-5,...
    'linewidth',2,'color',[0,0,255]/255)
plot(tef(tef>=datenum(2020,10,01) & tef<=datenum(2021,10,1)),...
    detrend(e3(tef>=datenum(2020,10,01) & tef<=datenum(2021,10,1)))-10,...
    'linewidth',2,'color',[100,149,237]/255)
plot(tgf(tgf>=datenum(2020,10,01) & tgf<=datenum(2021,10,1)),...
    detrend(g1(tgf>=datenum(2020,10,01) & tgf<=datenum(2021,10,1)))-25,...
    'linewidth',2,'color',[139,0,0]/255)
plot(tgf(tgf>=datenum(2020,10,01) & tgf<=datenum(2021,10,1)),...
    detrend(g2(tgf>=datenum(2020,10,01) & tgf<=datenum(2021,10,1)))-30,...
    'linewidth',2,'color',[255,0,0]/255)
legend('U1518-a','U1518-b','U1518-c','U1519-a','U1519-b','location','northwest')
ylim([-40 10])
datetick('x',3)
ylabel('\DeltaP (cm)')
set(gca,'fontsize',14)
box on; grid on

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print('../figures/exploratory/HOBITSS_VII/differences/CORK-formation_20-21','-dpng','-r300')

% combine and save data
tf={tef,tgf};
pf={ef,gf};
p1={e1,g1};
p2={e2,g2};
p3={e3,[]};
save('../processed_data/CORK_formation','tf','pf','p1','p2','p3','staname','stalat','stalon','stadepth')
