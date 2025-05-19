% import_CORK_formation.m
%
% Read in CORK data from text files and use wellhead pressure to correct
% for oceanographic effects
%

clear; close all

%% station info
staname={'U1518','U1519'};
stalat=[-38.859,-38.727]; % BPR location approximated from map view
stalon=[178.896,178.615]; % BPR location approximated from map view
stadepth=[0,0]; % what are the sensor depths?

%% read in the data from text files
%---U1518 Part 1---%
% read in data
% you'll have to change the directory and file names to wherever you have these data saved
fid = fopen('/Volumes/TOSHIBA_EXT/apg_data/CORK_LTBPR/U1518/U1518_RR1902_clean.dat','r');
fspec = '%s %s %f %f %f %f %f %f %f %f %f'; % [t t T T1 P1 T2 P2 T3 P3 T4 P4]
E = textscan(fid,fspec,'HeaderLines',3);
fclose(fid);
% convert timestamp to MATLAB time format
t_str=cat(2,cat(1,E{1}{:}),repmat(' ',length(E{1}),1),cat(1,E{2}{:}));
tE = datetime(t_str,'InputFormat','yyyy-MM-dd HH:mm:ss');
% merge P and T data together into a single matrix
E = [E{5},E{4},E{7},E{6},E{9},E{8},E{11},E{10}]; % [P1 T1 P2 T2 P3 T3 P4 T4]

%---U1518 Part 2---%
% read in data
% change directory info
fid = fopen('/Volumes/TOSHIBA_EXT/apg_data/CORK_LTBPR/U1518/U1518_TAN2102.dat','r');
fspec = '%s %s %f %f %f %f %f %f %f %f %f'; % [t t T T1 P1 T2 P2 T3 P3 T4 P4]
F = textscan(fid,fspec);
fclose(fid);
% convert timestamp to MATLAB time format
t_str=cat(2,cat(1,F{1}{:}),repmat(' ',length(F{1}),1),cat(1,F{2}{:}));
tF = datetime(t_str,'InputFormat','yyyy-MM-dd HH:mm:ss');
F{11} = [F{11}(1:525835);F{11}(525836:end)-F{11}(525836)+F{11}(525835)+...
    mean(diff(F{11}(525820:525835)))]; % corrects an artificial offset
% append everything onto the previous matrices 
tE = [tE;tF];
E = cat(1,E,[F{5},F{4},F{7},F{6},F{9},F{8},F{11},F{10}]); % [P1 T1 P2 T2 P3 T3 P4 T4]

%---U1518 Part 3---%
% read in data
% change directory info
fid = fopen('/Volumes/TOSHIBA_EXT/apg_data/CORK_LTBPR/U1518/U1518_2023_formatted.dat','r');
fspec = '%s %s %f %f %f %f %f %f %f %f %f'; % [t t T T1 P1 T2 P2 T3 P3 T4 P4]
K = textscan(fid,fspec);
fclose(fid);
% convert timestamp to MATLAB time format
t_str=cat(2,cat(1,K{1}{:}),repmat(' ',length(K{1}),1),cat(1,K{2}{:}));
tK = datetime(t_str,'InputFormat','yyyy-MM-dd HH:mm:ss');
% append everything onto the previous matrices
tE = [tE;tK];
% I had to adjust the units for this one, and account for an offset in the wellhead pressure
E = cat(1,E,[K{5}/10,K{4},K{7}/10,K{6},K{9}/10,K{8},K{11}/10+0.15,K{10}]); % [P1 T1 P2 T2 P3 T3 P4 T4]

% empirical trimming of blips in the data
tE = tE(3123:end);
E = E(3123:end,:);
E(92703:92706,7)=linspace(E(92703,7),E(92706,7),4);

%% oceanographic correction

%----- The Theory -----%
% everything from here through the dashed lines below can be deleted without
% affecting the rest of the code

% sub-select just a bit of the data as an example
i1=100000; i2=110000;
p1 = E(i1:i2,1); % this is the deepest formation pressure
p2 = E(i1:i2,7); % this is the wellhead pressure
t = datenum(tE(i1:i2)); % this converts time to a number that we can manipulate

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

% Like we talked about, the correlation between the two pressure records is visually obvious.
% Now here's the math to calculate a "linear inversion" that maps the wellhead pressure to
% the formation pressure.
% We want to solve p1 = c + b*t + a*p2, where c is a constant that accounts for the different
% base pressure between sensors, b*t is a linear term that accounts for differences in the
% sensor drift between the two, and a*p2 attempts to scale the wellhead pressure to the formation.
% We already have p1, p2, and t from the dataset but need to determine a, b, and c.
% This is a  system of equations with 3 unknowns but the dataset's length of equations, so it is
% an overdetermined system. This means that we can't solve it exactly, but instead we can
% calculate values of a, b, and c that minimize the misfit (or difference) between the two sides
% of the equation.

tinv = (t-t(1))/max(t-t(1)); % it is computationally easier to use a low-order number system
G = [ones(size(t)),tinv,p2];

% Using the above, we can rewrite our equation in matrix form as:
% p1 = G * m, where m=[c;b;a]
% So to solve for m, we just do p1/G = m
% But it just looks a bit more complicated because we have to use valid
% matrix operations:

m=inv(G'*G)*G'*p1;

disp(['c = ' num2str(m(1)) ' dbar']) % offset between our pressures
disp(['b = ' num2str(m(2)) ' dbar/day']) % difference in slope
disp(['a = ' num2str(m(3))]) % scaling factor between pressure amplitudes

% we can then create what's called a "model" for the formation pressure
% (p1) from a, b, c, and p2:

p1_model = G*m; % you'll recognize that this is our equation from above
% another way to write the same thing is:
p1_model_alt = m(1)*G(:,1) + m(2)*G(:,2) + m(3)*G(:,3);
% which is the same as:
p1_model_alt_2 = m(1)*ones(size(t)) + m(2)*tinv + m(3)*p2;
% This is useful if you only want to correct your data with some part of
% the model rather than the entire thing. For example, you might want to
% leave in the constant offset and linear part of the original signal and
% just correct for the wellhead pressure:
p1_corrected = p1-m(3)*p2;

figure(12); clf; hold on
plot(t,p1,'o-')
plot(t,p1_model,'g^-')
plot(t,p1-p1_model+mean(p1),'kx-')
legend('formation pressure','modeled formation pressure','residuals')
ylabel('Pressure (dbar)')
datetick('x')
set(gca,'fontsize',14)
box on; grid on

% The "residuals" are the difference between the observed formation
% pressure (p1) and our model (p1_model). On the plot, you can see that
% these are very small compared to the variation in p1, so we've come up
% with a good model. There are various ways to quantify the residuals into
% a single goodness-of-fit number. A good one is the root mean square
% (RMS). The closer to zero the better

disp(['Misfit RMS = ' num2str(rms(p1-p1_model)) ' dbar'])

%-----------------------%

% now we'll just do the above for the entire dataset

tinv=datenum(tE-tE(1)); tinv=tinv/max(tinv);
% Here, I'm choosing a subset of the data that excludes the majority of the
% exponential seen at the beginning. Since we don't have an exponential
% term in our model, we have no ability to fit/correct it.
Gl=[ones(size(tinv(1000000:2500000))),tinv(1000000:2500000),E(1000000:2500000,7)];
ml=inv(Gl'*Gl)*Gl'*E(1000000:2500000,1);

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

%% final thoughts

% You can see in the plots of the corrected pressure that there are still a
% number of features that could be cleaned up, like the exponential and the
% blips where the different data files were stitched together. I'd be
% happy to chat about these things down the road, but I think this is
% plenty for you to go off of for now.

% Use the above code as a guide to do the same import and processing for
% the 1519 data and reach out if you hit any snags.