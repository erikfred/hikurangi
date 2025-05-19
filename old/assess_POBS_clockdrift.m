% assess_POBS_clockdrift.m
%
% Look at 2K counting clock and compare against nominal timestamp from GPS
% clock. Depending on results, may have to build something into
% 'process_POBS.m' to account for this effect
%

clear; close all

%%-----SETUP-----%%

% easier to list these manually so I don't have to re-run all of them every
% single time
instruments={'POBS1','POBS2','POBS3','POBS4','POBS5','POBS7_co-located-w-QA15',...
    'POBS8','POBS11','POBS12','POBS13','POBS14_QA15','POBS15','POBS16'};

instruments2={'POBS09','POBS10'};

%%-----END SETUP-----%%

for i=1:length(instruments)
    % make a directory for things to go in
    if ~exist(['../processed_POBS/' instruments{i}],'dir')
        mkdir(['../processed_POBS/' instruments{i}])
    end
    
    %---import GPS
    try
        Logs=importGPS_EKF('/Volumes/Spare/POBS-data',instruments{i});
    catch
        warning('backtrace','off')
        warning('ImportPobs:GPSReadFail','ImportPobs: Failed to parse GPS files!')
        warning('backtrace','on')
    end

    %--- Determine clock/count conversion

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

    if nosync
        disp([newline instruments{i} ':' newline 'No final GPS sync; assumed no clock drift'...
            newline 'ti = ' datestr(tc(1)) newline 'tf = ' datestr(tc(end))])
        cdrift=0;
        keyboard
    else
        disp([newline instruments{i} ':' newline 'ti = ' datestr(tc(1))...
            newline 'tf = ' datestr(tc(end)) newline 'Clock delta = ' num2str(dateend-datestart)...
            newline 'Count delta = ' num2str(tc(end)-tc(1)) newline 'delta delta = '...
            num2str(((dateend-datestart)-(tc(end)-tc(1)))*86400) ' seconds'])
        cdrift=(tc(end)-tc(1))-(dateend-datestart); % in days
        keyboard
    end

    % this is how I currently determine time
    % apg_temp(:,1)=datestart+(apg_temp(:,1)-countstart)/2000/60/60/24;
    % but I should instead do this to correct clockdrift (if any)
    % apg_temp(:,1)=datestart-countstart+apg_temp(:,1)/2000/60/60/24;
    % apg_temp(:,1)=apg_temp(:,1)-linspace(0,cdrift,length(apg_temp))';
end

for i=1:length(instruments2)
    % make a directory for things to go in
    if ~exist(['../processed_POBS/' instruments2{i}],'dir')
        mkdir(['../processed_POBS/' instruments2{i}])
    end
    
    %---import GPS
    try
        Logs=importGPS_EKF('/Volumes/Spare/POBS-data_2024_update',instruments2{i});
    catch
        warning('backtrace','off')
        warning('ImportPobs:GPSReadFail','ImportPobs: Failed to parse GPS files!')
        warning('backtrace','on')
    end

    %--- Determine clock/count conversion

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

    if nosync
        disp([newline instruments2{i} ':' newline 'No final GPS sync; assumed no clock drift'...
            newline 'ti = ' datestr(tc(1)) newline 'tf = ' datestr(tc(end))])
        cdrift=0;
        keyboard
    else
        disp([newline instruments2{i} ':' newline 'ti = ' datestr(tc(1))...
            newline 'tf = ' datestr(tc(end)) newline 'Clock delta = ' num2str(dateend-datestart)...
            newline 'Count delta = ' num2str(tc(end)-tc(1)) newline 'delta delta = '...
            num2str(((dateend-datestart)-(tc(end)-tc(1)))*86400) ' seconds'])
        cdrift=(tc(end)-tc(1))-(dateend-datestart); % in days
        keyboard
    end

    % this is how I currently determine time -- do I need to change?
    % apg_temp(:,1)=datestart+(apg_temp(:,1)-countstart)/2000/60/60/24;
    % but I should instead do this to correct clockdrift (if any)
    % apg_temp(:,1)=datestart-countstart+apg_temp(:,1)/2000/60/60/24;
    % apg_temp(:,1)=apg_temp(:,1)-linspace(0,cdrift,length(apg_temp))';
end