function Logs=importGPS_EKF(inpath,inname)
%% This funciton concatenates all the gpslog files into single file
% Any existing files are overwritten
%% Persistent lastFilePath variable sotes recently used path for use as starting path
persistent lastFilePath;
if isempty(lastFilePath)
    lastFilePath = pwd;
elseif lastFilePath == 0
    lastFilePath = pwd;
end

%% Housekeeping - Check inputs and prompt user if needed

startFileNum = 0;

fList = dir([inpath,filesep,inname,filesep,'GPS*']);
% Remove any directories from the list
fList = fList(~[fList.isdir]);
% Remove any empty files from the list
fList = fList([fList.bytes]~=0);

if(numel(fList) == 0)
    warning('gpslogcat:noGPSLOGFiles','No GPSLOG files found in %s',[inpath,filesep,inname]);
    return;
end

% Extract file number from the Tax file name 
% e.g. for file TAX00004.001, file number is 00004
fnum = zeros(size(fList));
for ff = 1:numel(fList)
    [~,fileNameOnly,~] = fileparts(fList(ff).name);
    fnum(ff) = str2double(fileNameOnly(4:end));
end

% Remove any files with unrecognized file name format from the list
if(any(isnan(fnum)))
    warning('readPobsText:BadFileName','File name format not Recognized %s',fList(isnan(fnum)).name);
    fList(isnan(fnum)) = [];
    fnum(isnan(fnum)) = [];
end

% Sort the files in numerical order
[fnum,sx] = sort(fnum);
fList = fList(sx);

% Remove any empty files from the list
if(any([fList.bytes] == 0))
    emptyfiles = find([fList.bytes] == 0);
    for eef = emptyfiles
        warning('backtrace','off');
        warning('readPobsText:ZeroSizeFile','Found empty file: %s',fList(eef).name);
        warning('backtrace','off');
    end
    fList(emptyfiles) = [];
    fnum(emptyfiles) = [];
end

% Remove files before startfile from the list
fList(fnum<startFileNum) = [];
fnum(fnum<startFileNum) = [];

numfiles = numel(fnum);
fprintf('%d GPS files to read.\n',numfiles);

%% Start Processing
[~,firstFileNameOnly,ext] = fileparts(fList(1).name);
batchfilename = ['../processed_POBS_Y2/',inname,filesep,'templogs',ext];
fileIDwr = fopen(batchfilename,'w');
writefilecleaner = onCleanup(@() fclose(fileIDwr));

% dumps parsed data into a single file to simplify reading/converting later
for fileNum = 1:length(fList)
    filename = [inpath,filesep,inname,filesep,fList(fileNum).name];
    fprintf('%s\t file %d of %d.\n',fList(fileNum).name,fileNum,numfiles)

    fileIDrd = fopen(filename,'r');
    % Create cleanup object in case of crash
    readfilecleaner = onCleanup(@() fclose(fileIDrd));
    fwrite(fileIDwr,fread(fileIDrd));
    clear readfilecleaner;
end
clear writefilecleaner;

% Initialize variables.
maxYears = 5;  % filter out erroneous values above this
delimiter = ',';
if nargin<=2
    startRow = 1;
    endRow = inf;
end

% Read columns of data as text:
formatSpec = '%s%s%s%s%s%s%s%[^\n\r]';
% R/G,RTCTIME/GPSTIME,PPS_RE_2kCountH,L,LastSampleTime2kCountH,L,StatusString

%% Open the text file.
fileID = fopen(batchfilename,'r');

% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

% Close the text file.
fclose(fileID);

% Convert the contents of columns containing numeric text to numbers.
% Replace non-numeric text with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

% Convert the contents of columns with dates to MATLAB datetimes using the
% specified date format.
try
    dates{2} = datetime(dataArray{2}, 'Format', 'yyMMddHHmmss', 'InputFormat', 'yyMMddHHmmss');
catch
    try
        % Handle dates surrounded by quotes
        dataArray{2} = cellfun(@(x) x(2:end-1), dataArray{2}, 'UniformOutput', false);
        dates{2} = datetime(dataArray{2}, 'Format', 'yyMMddHHmmss', 'InputFormat', 'yyMMddHHmmss');
    catch
        dates{2} = repmat(datetime([NaN NaN NaN]), size(dataArray{2}));
    end
end

dates = dates(:,2);

% Split data into numeric and string columns.
rawNumericColumns = {};
rawStringColumns = string(raw(:, [1,3,4,5,6,7]));


% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
rawNumericColumns(R) = {NaN}; % Replace non-numeric cells

% Make sure any text containing <undefined> is properly converted to an <undefined> categorical
idx = (rawStringColumns(:, 1) == "<undefined>");
rawStringColumns(idx, 1) = "";

% Create output variable
firstLineDataEmpty = cellfun(@isempty,rawStringColumns(1,:));
lastLineDataEmpty = cellfun(@isempty,rawStringColumns(end,:));

if all(lastLineDataEmpty(4:6))
    rawStringColumns(end,:) = [];
    dates{1}(end) = [];
end
if all(firstLineDataEmpty(4:6))
    rawStringColumns(1,:) = [];
    dates{1}(1) = [];
end

Logs = table;
Logs.ClockType = categorical(rawStringColumns(:, 1));
Logs.Time = dates{:, 1};
Logs.PPS_RE2kCnt = anybaseString2num(82,rawStringColumns(:,2),45)*2^32 + anybaseString2num(82,rawStringColumns(:,3),45);
Logs.Last_Sample2kCnt = anybaseString2num(82,rawStringColumns(:,4),45)*2^32 + anybaseString2num(82,rawStringColumns(:,5),45);
Logs.Status = rawStringColumns(:, 6);

% Some QC
Logs.ClockType(isundefined(Logs.ClockType)) = 'Status';
Logs.PPS_RE2kCnt(Logs.PPS_RE2kCnt <= 0 | Logs.PPS_RE2kCnt > 2000*60*60*24*365*maxYears) = nan;
Logs.Last_Sample2kCnt(Logs.Last_Sample2kCnt <= 0 | Logs.PPS_RE2kCnt > 2000*60*60*24*365*maxYears) = nan;

save(['../processed_POBS_Y2/',inname,filesep,firstFileNameOnly,'.mat'],'Logs');