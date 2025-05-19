% process_POBS_2024_additions.m
%
% Convert POBS data from pseudo-binary to text/csv at 1 Hz
%
% This is just a dupe of the original with directory tweaks to point only
% to new data that was recovered in June 2024
%

clear; close all

%%-----SETUP-----%%

% easier to list these manually so I don't have to re-run all of them every
% single time
instruments={'POBS09','POBS10'};

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

%---FIGURE OUT LATER FOR DAILY FILES
%---DOESN'T MAKE SENSE TO DECIMATE SEISMIC DATA
% % Triaxial accelerometer data format
% % Typical: ---,gF5---,*0001,-9.6824911,-1.4282716,-.5340717,24.648580880
% % WIth Time Stamp: ---,gF5---,*0001,V,01/01/70 00:00:03.000,-9.6824911,-1.4282716,-.5340717,24.648580880
% % All data frames end with <CR><LF> i.e. 0x0D0A
% taxFile.formatSpec = '%s%s%s%f%f%f%f';
% taxFile.numFields = 7;
% taxFile.nonDataFields = 3;
% taxFile.type = 'TAX*';
% taxFile.fieldNames = {'TAX_XAcc','TAX_YAcc','TAX_ZAcc','TAX_Temp'};
% taxFile.startFileNum = startfilenum;

%%-----END SETUP-----%%

for i=1:length(instruments)
    % make a directory for things to go in
    if ~exist(['../processed_POBS/' instruments{i}],'dir')
        mkdir(['../processed_POBS/' instruments{i}])
    end
    
    %---import GPS
    try
        Logs=importGPS_EKF('/Volumes/Spare/POBS-data_2024_update',instruments{i});
    catch
        warning('backtrace','off')
        warning('ImportPobs:GPSReadFail','ImportPobs: Failed to parse GPS files!')
        warning('backtrace','on')
    end
    
    %--- Determine clock/count conversion (overkill, but eases troubleshooting)

    % extract relevant columns
    clocktype=table2array(Logs(:,1));
    lc=floor(length(clocktype)/2);
    ct1=clocktype(1:lc);
    ct2=clocktype(lc+1:end);

    iG1=find(ct1=='G',1,'last');
    iG2=find(ct2=='G')+length(ct1);

    if length(iG2)>1
        nosync=false;
        iG2a=iG2(1); iG2b=iG2(2);
        % sometimes first GPS sync upon recovery writes cached (?) time from deployment
        Logs(iG2a,:)=[]; % just delete that row
        iG2b=iG2b-1; % account for deletion
    elseif length(iG2)==1
        nosync=false;
        iG2b=iG2;
    elseif isempty(iG2)
        nosync=true;
        iG2b=height(Logs);
    end

    tg=datenum(table2array(Logs(iG1:iG2b,2))); % GPS clock
    c=table2array(Logs(iG1:iG2b,4))/2000/60/60/24; % 2K count converted to d
    tc=tg(1)+c-c(1); % pseduo-clock from 2K count

    % temp=table2array(Logs(iG1:iG2b,3)); % the unknown column

    % remove NaNs, ensure parity
    inan_g=find(isnan(tg));
    inan_c=find(isnan(tc));
    if any(inan_c~=inan_g)
        warning('NaNs did not match; forced equal')
        inan_temp=unique([inan_g;inan_c]);
        inan_g=inan_temp;
        inan_c=inan_temp;
    end
    tg(inan_g)=[];
    tc(inan_c)=[];
    c(inan_c)=[];
    % temp(inan_c)=[];

    datestart=tg(1);
    dateend=tg(end);
    countstart=c(1);
    countend=c(end);
    
    %---import APG
    % pare down file list to simplify loops
    flist_all=dir(['/Volumes/Spare/POBS-data_2024_update/' instruments{i}]);
    fname_all={flist_all.name}';
    ia=strncmp(fname_all,'A',1);
    flist=flist_all(ia);

    fbytes=cat(1,flist.bytes);
    ib=fbytes~=0;
    flist=flist(ib);

    % declare decimation variables
    ta=[];
    a1=[]; a2=[];
    Ta1=[]; Ta2=[];
    apg_hold=[];

    % loop through APG files
    for j=1:length(flist)
        if strcmp(flist(j).name,'APG.zip')
            continue
        end

        fwork=[flist(j).folder '/' flist(j).name];
        
        apgFile.inst=instruments{i};
        apgFile.loc=fwork;

        % external function handles the reading
        apg_temp=readPobsText_EKF(apgFile,debugEnable);
        
        % convert tables to matrices for manipulation
        apg_temp=table2array(apg_temp);
        apg_temp(:,1)=datestart+(apg_temp(:,1)-countstart)/2000/60/60/24;

        % remove anomalous entries and note indices
        if i==1 && j<4 % I don't think these data are real, and they cause problems
            continue
        end
        timecheck=flip(find(diff(apg_temp(:,1))<0),1);
        if ~isempty(timecheck)
            for k=1:length(timecheck)
                warning('i = %d, j = %d. apg_temp(%d,1) = %f manually removed. \n',i,j,timecheck(k),apg_temp(timecheck(k),1))
                apg_temp(timecheck(k),:)=[];
            end
        end
        
        if ~isempty(apg_hold)
            apg_temp=cat(1,apg_hold,apg_temp);
            apg_hold=[];
        end

        % Decimate to 1Hz
        i1 = 1;
        nextday = floor(apg_temp(1,1))+1;
        while i1<length(apg_temp)
            i2 = find(apg_temp(:,1)>=nextday,1);
            if isempty(i2)
                if j~=length(flist)
                    apg_hold=apg_temp(i1:end,:);
                    break
                else
                    i2=length(apg_temp);
                end
            end
            [segt,segA,~,~] = downsample_uneven(apg_temp(i1:i2-1,1),apg_temp(i1:i2-1,2:5),1/24/60/60);
            if length(segt)>24*60*60
                keyboard
            end
            ta=cat(1,ta,segt);
            a1=cat(1,a1,segA(:,1));
            Ta1=cat(1,Ta1,segA(:,2));
            a2=cat(1,a2,segA(:,3));
            Ta2=cat(1,Ta2,segA(:,4));
            i1=i2;
            nextday=floor(apg_temp(i2,1))+1;
        end
    end

    % correct clock drift, if any
    if ~nosync
        cdrift=(tc(end)-tc(1))-(dateend-datestart); % in days
        ta=ta-linspace(0,cdrift,length(ta))';
        
        % interpolate onto 1 Hz timestamp
        ttemp=(ta(1):datenum(0,0,0,0,0,1):ta(end))';
        a1=interp1(ta,a1,ttemp);
        a2=interp1(ta,a2,ttemp);
        Ta1=interp1(ta,Ta1,ttemp);
        Ta2=interp1(ta,Ta2,ttemp);
        ta=ttemp;
    end

    %---import BARO
    % pare down file list to simplify loops
    ia=strncmp(fname_all,'B',1);
    flist=flist_all(ia);

    fbytes=cat(1,flist.bytes);
    ib=fbytes~=0;
    flist=flist(ib);

    % declare decimation variables
    tb=[]; b=[]; Tb=[];

    % loop through BARO files
    for j=1:length(flist)
        if strcmp(flist(j).name,'BAR.zip')
            continue
        end

        fwork=[flist(j).folder '/' flist(j).name];
        
        baroFile.inst=instruments{i};
        baroFile.loc=[fwork(1:end-12) 'BAR' fwork(end-8:end)];
        % taxFile.inst=instruments{i};
        % taxFile.loc=[fwork(1:end-12) 'TAX' fwork(end-8:end)];

        % external function handles the reading
        bar_temp=readPobsText_EKF(baroFile,debugEnable);
        if isempty(bar_temp)
            continue
        end
        % tax_temp=readPobsText(taxFile,debugEnable);
        
        % convert tables to matrices for manipulation
        bar_temp=table2array(bar_temp);
        bar_temp(:,1)=datestart+(bar_temp(:,1)-countstart)/2000/60/60/24;
        % tax_temp=table2array(tax_temp);
        % tax_temp(:,1)=datestart+(tax_temp(:,1)-countstart)/2000/60/60/24;

        % Barometer already at ~1Hz
        tb=cat(1,tb,bar_temp(:,1));
        b=cat(1,b,bar_temp(:,2));
        Tb=cat(1,Tb,bar_temp(:,3));
    end

    %--- barometer more complicated because it isn't continuous
    if ~nosync
        cdrift=(tc(end)-tc(1))-(dateend-datestart); % in days
        btemp=b; Tbtemp=Tb;

        % find "irregular" time steps, make NaNs
        ireg=find(abs(diff(tb))>0.0001);
        btemp(ireg)=NaN;
        Tbtemp(ireg)=NaN;

        % now interpolate to 1 Hz (NaNs will propagate)
        ttemp=(tb(1):datenum(0,0,0,0,0,1):tb(end))';
        b=interp1(tb,btemp,ttemp);
        Tb=interp1(tb,Tbtemp,ttemp);

        % correct drift at 1 Hz
        ttemp=ttemp-linspace(0,cdrift,length(ttemp))';

        % interpolate back to 1 Hz
        tb=(ttemp(1):datenum(0,0,0,0,0,1):ttemp(end))';
        b=interp1(ttemp,b,tb);
        Tb=interp1(ttemp,Tb,tb);
        
        % remove NaNs to get back to "normal" sampling
        tb(isnan(b))=[];
        Tb(isnan(b))=[];
        b(isnan(b))=[];
    end
    
    % save decimated data locally
    save(['../processed_POBS/' instruments{i} '/' instruments{i} '.mat'],'ta','a1','a2','Ta1','Ta2','tb','b','Tb')
end