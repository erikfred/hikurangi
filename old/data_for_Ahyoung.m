% data_for_Ahyoung.m
%
% Load up relevant data, decimate as needed, write to text file for easy
% sharing
%

clear; close all

dest='../processed_data/for_Ahyoung/';

%% 2023-2024 GNS data

load('../processed_data/GONDOR_II.mat');

% destination to write
sdir=[dest 'GNS'];
if ~exist(sdir,'dir')
    mkdir(sdir)
end

% write 1/hr data to text file
for i=1:length(staname)
    sname=[sdir '/' staname{i}];

    tstr=datestr(tf{i}',30);
    X=table(tstr,pf{i}'*100,Tf{i}');
    X.Properties.VariableNames(1:3)={'t (yyyymmddTHHMMSS)','P (Pa)','T (C)'};
    writetable(X,sname)
end

%% 2023-2024 GNS data

POBS_dir=dir('../../A-0-A/stitched_data_Y2/drift_corrected');
flist={POBS_dir.name}';
fcheck=cellfun(@(v)v(1),flist);
i_list=find(eq(fcheck,'P'));
POBS_list=flist(i_list);
POBS_list(strcmp(POBS_list,'POBS09.mat') | strcmp(POBS_list,'POBS10.mat'))=[];

staname={'LDE23-AB';'LDE23-AN';'LDE23-AE';'LDE23-AD';'LDE23-AK';'LDE23-AM';...
    'LDE23-AP';'LDE23-AC';'LDE23-AL';'LDE23-AJ';'LDE23-AI';'LDE23-AG'};

% destination to write
sdir=[dest 'LDEO'];
if ~exist(sdir,'dir')
    mkdir(sdir)
end

% write 1/hr data to text file
for i=1:length(POBS_list)
    sname=[sdir '/' staname{i}];

    load([POBS_dir(1).folder '/' POBS_list{i}])

    tstr=datestr(dataf.tf,30);
    sname1=[sname '_gauge1'];
    X=table(tstr,dataf.p1_dcor*100,dataf.T1f);
    X.Properties.VariableNames(1:3)={'t (yyyymmddTHHMMSS)','P (Pa)','T (C)'};
    writetable(X,sname1)
    sname2=[sname '_gauge2'];
    Y=table(tstr,dataf.p2_dcor*100,dataf.T2f);
    Y.Properties.VariableNames(1:3)={'t (yyyymmddTHHMMSS)','P (Pa)','T (C)'};
    writetable(Y,sname2)
end