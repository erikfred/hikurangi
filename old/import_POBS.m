% import_POBS.m
%
%

clear; close all

POBS_dir=dir('../apg_data/2018_2019_Gisborne_HawkesBay/POBSdataforLaura/APGmat3/');

A=[];
for i=3:13 %length(POBS_dir)
    load([POBS_dir(i).folder '/' POBS_dir(i).name])

    segA=table2array(pObsData);
    segA(:,1)=segA(:,1)/2000; % converts to seconds

    A=cat(1,A,segA);
end

figure(32); clf; hold on
plot(A(:,1),A(:,2))
plot(A(:,1),A(:,4))