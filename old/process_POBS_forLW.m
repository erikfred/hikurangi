% process_POBS.m
%
% Convert POBS data from pseudo-binary
%

%%-----SETUP-----%%

importAPG=true;
importTAX=true;

% my directories are structured as /POBS-data/POBS1/, /POBS-data/POBS2, etc.
rawdir='/Volumes/Spare/POBS-data';
instruments={'POBS1','POBS5'};
debugEnable = false;

% Apg data format
% Typical: ---,-@5---,*0001,14.06051,22.859413419,14.13115,22.899507
% With Time Stamp: ---,-@5---,*0001,V,01/01/70 00:00:03.000,14.06051,22.859413419,14.13115,22.899507
% All data frames end with <CR><LF> i.e. 0x0D0A

apgFile.formatSpec = '%s%s%s%f%f%f%f';
apgFile.numFields = 7;
apgFile.nonDataFields = 3;
apgFile.type = 'APG*';
apgFile.fieldNames = {'APG1_Pres','APG1_Temp','APG2_Pres','APG2_Temp'};

% Barometer data format
% Typical: ---,[V5---,*0001,14.44610,24.2647
% With Time Stamp: ---,[V5---,*0001,V,01/01/70 00:00:03.000,14.44610,24.2647
% All data frames end with <CR><LF> i.e. 0x0D0A

baroFile.formatSpec = '%s%s%s%f%f';
baroFile.numFields = 5;
baroFile.nonDataFields = 3;
baroFile.type = 'BAR*';
baroFile.fieldNames = {'Baro_Pres','Baro_Temp'};

% Triaxial accelerometer data format
% Typical: ---,gF5---,*0001,-9.6824911,-1.4282716,-.5340717,24.648580880
% WIth Time Stamp: ---,gF5---,*0001,V,01/01/70 00:00:03.000,-9.6824911,-1.4282716,-.5340717,24.648580880
% All data frames end with <CR><LF> i.e. 0x0D0A

taxFile.formatSpec = '%s%s%s%f%f%f%f';
taxFile.numFields = 7;
taxFile.nonDataFields = 3;
taxFile.type = 'TAX*';
taxFile.fieldNames = {'TAX_XAcc','TAX_YAcc','TAX_ZAcc','TAX_Temp'};
taxFile.startFileNum = startfilenum;

%%-----END SETUP-----%%

for i=1:length(instruments)
    
    %---import APG
    % pare down file list to simplify loops
    flist_all=dir([rawdir '/' instruments{i}]);
    fname_all={flist_all.name}';
    ia=strncmp(fname_all,'A',1);
    flist=flist_all(ia);

    fbytes=cat(1,flist.bytes);
    ib=fbytes~=0;
    flist=flist(ib); % this is now the list of all APG files with data

    if importAPG
        j1=1; % first non-empty file
        j2=length(flist); % last non-empty file

        fwork=[flist(j1).folder '/' flist(j1).name];
        
        apgFile.inst=instruments{i};
        apgFile.loc=fwork;

        % external function handles the reading
        apg_temp=readPobsText_EKF(apgFile,debugEnable);
        
        % convert tables to matrices for manipulation
        apg_temp=table2array(apg_temp);
        apg_temp(:,1)=apg_temp(:,1)/2000/60/60/24; % converts 2k count to days

        figure(54); clf;
        subplot(211); hold on
        plot(apg_temp(:,1),apg_temp(:,2),'.')
        ylabel('P1')
        yyaxis right
        plot(apg_temp(:,1),apg_temp(:,3),'.')
        ylabel('T1')
        title([instruments{i} ', ' flist(j1).name ' [APG]'])
        subplot(212); hold on
        plot(apg_temp(:,1),apg_temp(:,4),'.')
        ylabel('P2')
        yyaxis right
        plot(apg_temp(:,1),apg_temp(:,5),'.')
        ylabel('T2')
        
        %---same thing for the final datafile

        fwork=[flist(j2).folder '/' flist(j2).name];
        
        apgFile.inst=instruments{i};
        apgFile.loc=fwork;

        % external function handles the reading
        apg_temp2=readPobsText_EKF(apgFile,debugEnable);
        
        % convert tables to matrices for manipulation
        apg_temp2=table2array(apg_temp2);
        apg_temp2(:,1)=apg_temp2(:,1)/2000/60/60/24; % converts 2k count to days

        figure(54); clf;
        subplot(211); hold on
        plot(apg_temp2(:,1),apg_temp2(:,2),'.')
        ylabel('P1')
        yyaxis right
        plot(apg_temp2(:,1),apg_temp2(:,3),'.')
        ylabel('T1')
        title([instruments{i} ', ' flist(j2).name ' [APG]'])
        subplot(212); hold on
        plot(apg_temp2(:,1),apg_temp2(:,4),'.')
        ylabel('P2')
        yyaxis right
        plot(apg_temp2(:,1),apg_temp2(:,5),'.')
        ylabel('T2')
    end

    %---import TAX
    % pare down file list to simplify loops
    ia=strncmp(fname_all,'T',1);
    flist=flist_all(ia);

    fbytes=cat(1,flist.bytes);
    ib=fbytes~=0;
    flist=flist(ib);

    if importTAX
        j1=1; % first non-empty file
        j2=length(flist); % last non-empty file

        fwork=[flist(j1).folder '/' flist(j1).name];
        
        taxFile.inst=instruments{i};
        taxFile.loc=fwork;

        % external function handles the reading
        tax_temp=readPobsText_EKF(taxFile,debugEnable);
        
        % convert tables to matrices for manipulation
        tax_temp=table2array(tax_temp);
        tax_temp(:,1)=tax_temp(:,1)/2000/60/60/24; % converts 2k count to days

        keyboard % I don't recall how the TAX data are structured, so you'll have to look

        %---same thing for the final datafile

        fwork=[flist(j2).folder '/' flist(j2).name];
        
        taxFile.inst=instruments{i};
        taxFile.loc=fwork;

        % external function handles the reading
        tax_temp2=readPobsText_EKF(taxFile,debugEnable);
        
        % convert tables to matrices for manipulation
        tax_temp2=table2array(tax_temp2);
        tax_temp2(:,1)=tax_temp2(:,1)/2000/60/60/24; % converts 2k count to days

        keyboard % I don't recall how the TAX data are structured, so you'll have to look
    end
end