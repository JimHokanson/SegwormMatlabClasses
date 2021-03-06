function worm2stats(stats_output_file_path, wormFiles, varargin)
%worm2stats  Convert worms to a set of statistics.
%
%   seg_worm.w.stats.worm2stats(stats_output_file_path, wormFiles, *CONTROLFILES, *ISOLDCONTROL, *VERBOSE)
%
%   The stats are basically a subset of the histogram. For
%   non-histogram data (summary statstics for events), the data are
%   maintained, otherwise 'histogram' structures are replaced with
%   'statistics' structures.
%
%   The function that computes statistics (worm2StatsInfo) can use
%   either the histogram or statistics file. As such I think this file
%   MIGHT BE KIND OF USELESS.
%
%   Inputs:
%   =======================================================================
%       stats_output_file_path     - the file name for the statistics;
%                      the file includes:
%
%                      wormInfo    = the worm information
%                      worm        = the worm statistics
%                      controlInfo = the control information (if it exists)
%                      worm        = the control statistics (if they exist)
%
%       wormFiles    - the histogram files to use for the worm(s)
%       controlFiles - the histogram files to use for the control(s);
%                      if empty, the worm has no new control
%       isOldControl - are we adding the old controls?
%                      the default is yes (true)
%       isVerbose    - verbose mode display the progress;
%                      the default is no (false)
%
%   See also:
%   seg_worm.w.stats.worm2histogram
%   seg_worm.util.histogram,
%   WORMDISPLAYINFO, 
%   WORMDATAINFO
%
%
% � Medical Research Council 2012
% You will not remove any copyright or other notices from the Software;
% you must reproduce all copyright notices and other proprietary
% notices on any copies of the Software.

%Input Handling
%--------------------------------------------------------------------------
% Do we have a control?
controlFiles = [];
if ~isempty(varargin)
    controlFiles = varargin{1};
end

% Are we adding the old controls?
isOldControl = true;
if length(varargin) > 1
    isOldControl = varargin{2};
end

% Are we in verbose mode?
isVerbose = false;
if length(varargin) > 2
    isVerbose = varargin{3};
end
%--------------------------------------------------------------------------


% Organize the worm files.
if ~iscell(wormFiles)
    wormFiles = {wormFiles};
end
if ~isempty(controlFiles) && ~iscell(controlFiles)
    controlFiles = {controlFiles};
end

% Delete the file if it already exists.
if exist(stats_output_file_path, 'file')
    delete(stats_output_file_path);
end

% Save the worm information.
%---------------------------------------------------------------
if isVerbose
    disp('Combining "wormInfo" ...');
end
newWormInfo = [];
for i = 1:length(wormFiles)
    wormInfo = [];
    load(wormFiles{i}, 'wormInfo');
    newWormInfo = cat(1, newWormInfo, wormInfo);
end
wormInfo = newWormInfo;
save(stats_output_file_path, 'wormInfo');
clear wormInfo;




% Collect the new control information.
%---------------------------------------------------------------
if isVerbose
    disp('Combining "controlInfo" ...');
end
newControlInfo = [];
if ~isempty(controlFiles)
    for i = 1:length(controlFiles)
        wormInfo = [];
        load(controlFiles{i}, 'wormInfo');
        newControlInfo =  cat(1, newControlInfo, wormInfo);
    end
end

% Collect the old control information.
%---------------------------------------------------------------
if isOldControl
    for i = 1:length(wormFiles)
        controlInfo = who('-FILE', wormFiles{i}, 'controlInfo');
        if ~isempty(controlInfo)
            load(wormFiles{i}, 'controlInfo');
            newControlInfo =  cat(1, newControlInfo, controlInfo);
        end
    end
end

% Save the control information.
%---------------------------------------------------------------
if ~isempty(newControlInfo)
    controlInfo = newControlInfo;
    save(stats_output_file_path, 'controlInfo', '-append');
end
clear controlInfo;



%--------------------------------------------------------------------------

% Initialize the data information.
dataInfo = seg_worm.feature.roots();

% Save the worm statistics.
%--------------------------------------------------------------------------
h__saveStatistics(stats_output_file_path, wormFiles, dataInfo, 'worm', 'worm', isVerbose);

% Are we adding the old controls?
if isOldControl
    
    % Initialize the new control names.
    controlNames = cell(length(controlFiles), 1);
    for i = 1:length(controlFiles)
        controlNames{i} = 'worm';
    end
    
    % Initialize the old control names.
    for i = 1:length(wormFiles)
        control = who('-FILE', wormFiles{i}, 'control');
        if ~isempty(control)
            controlFiles{end + 1} = wormFiles{i};
            controlNames{end + 1} = 'control';
        end
    end
    
    % Initialize the new control names.
else
    controlNames = 'worm';
end

% Save the control statistics.
if ~isempty(controlFiles)
    h__saveStatistics(stats_output_file_path, controlFiles, dataInfo, controlNames, 'control', isVerbose);
end
end



%% Load worm data from files.
function data = h__loadWormFiles(filenames, wormName, field)
%
%   INPUTS
%   -----------------------------------------------------------------------
%
%   OUTPUTS
%   -----------------------------------------------------------------------
%   data : (cell)
%

% Fix the data.
if ~iscell(wormName)
    wormName = {wormName};
end

% Load each worm by name.
if length(wormName) > 1
    data = cell(length(wormName), 1);
    for i = 1:length(wormName)
        data{i} = seg_worm.util.loadStructField(filenames{i}, wormName{i}, field);
    end
    
    % Load all worms using the same name.
else
    data = cellfun(@(x) seg_worm.util.loadStructField(x, wormName{1}, field), filenames, 'un', 0);
end
end



%% Save the worm statistics.
function h__saveStatistics(filename, wormFiles, dataInfo, loadName, saveName, isVerbose)
%
%
%
%   THE TWO CALLS TO THIS ARE:
%   -----------------------------------------------------------------------
%   h__saveStatistics(filename, controlFiles, dataInfo, controlNames, 'control', isVerbose);
%   h__saveStatistics(filename, wormFiles,    dataInfo, 'worm',       'worm',    isVerbose);
%   h__saveStatistics(filename, wormFiles,    dataInfo, loadName,      saveName, isVerbose)
%
%   NOTE: The controlNames indicates that in a new file, the field that
%   will be the control is the 'worm' entry, in old files, it is identified
%   as 'control
%
%   i.e. features A - histogram A - 'worm'
%        features B - histogram B - 'worm'
%   
%   now, B might be our control for A
%
%   Alternatively, when doing histograms, B might be identified as a
%   control at that time, and subsequently needs to be accessed by
%   'control'
%
%   INPUTS
%   =======================================================================
%   dataInfo : 
%
%   See Also
%   seg_worm.feature.roots 


%NOTE: .field



% Combine the statistics.
for i = 1:length(dataInfo)
    field = dataInfo(i).field;
    if isVerbose
        disp(['Combining "' field '" ...']);
    end
    switch dataInfo(i).type
        
        % Combine simple statistics.
        case 's'
            data = h__addStatistics(wormFiles, loadName, field);
            
            eval([saveName '.' field '=data;']);
            
            % Combine motion statistics.
        case 'm'
            data = h__addMotionStatistics(wormFiles, loadName, field);
            
            eval([saveName '.' field '=data;']);
            
            % Combine event statistics.
        case 'e'
            
            % Combine the event data statistics.
            subFields = dataInfo(i).subFields.summary;
            for j = 1:length(subFields);
                subField = [field '.' subFields{j}];
                
                data = h__addEventStatistics(wormFiles, loadName, subField);
                
                eval([saveName '.' subField '=data;']);
            end
            
            % Combine the event statistics.
            subFields = dataInfo(i).subFields.data;
            for j = 1:length(subFields);
                subField = [field '.' subFields{j}];
                
                data = h__addStatistics(wormFiles, loadName, subField);
                
                eval([saveName '.' subField '=data;']);
            end
    end
end

% Save the statistics.
save(filename, saveName, '-append');
end



%% Combine statistics.
function data = h__addStatistics(wormFiles, wormName, field)
%
%
%   INPUTS
%   --------------------------------------------------------------
%   wormFiles :
%   wormName  : 
%   

%For simple data types ...

%???? what does addData contain?????
addData = h__loadWormFiles(wormFiles, wormName, field);

data(length(addData{1})).statistics = [];

for i = 1:length(addData{1})
    
    subData = cellfun(@(x) x(i).histogram, addData, 'un', 0);
    
    data(i).statistics = h__addStatisticsData(subData);
end

end



%% Combine motion statistics.
function data = h__addMotionStatistics(wormFiles, wormName, field)
%
%
%   Motion Statistics:
%   forward, paused, backward


% Initialize the locomotion modes.
motionNames = { ...
    'forward', ...
    'paused', ...
    'backward'};

% Initialize the data.
data.statistics = [];
for i = 1:length(motionNames)
    data.(motionNames{i}).statistics = [];
end

% Get the data.
addData = h__loadWormFiles(wormFiles, wormName, field);
if isempty(addData)
    return;
end

% Combine the statistics.
data(length(addData{1})).statistics = [];
for i = 1:length(addData{1})
    
    % Combine the data statistics.
    subData = cellfun(@(x) x(i).histogram, addData, 'un', 0);
    data(i).statistics = h__addStatisticsData(subData);
    
    % Combine the motion statistics.
    for j = 1:length(motionNames)
        subData = cellfun(@(x) x(i).(motionNames{j}).histogram, addData, 'un', 0);
        data(i).(motionNames{j}).statistics = h__addStatisticsData(subData);
    end
end
end



%% Combine statistics.
function data = h__addStatisticsData(addData)
%
%   INPUTS
%   ================================================
%   addData : 
%

% Is the data signed?
data     = [];
isSigned = [];
numSets  = 0;
for i = 1:length(addData)
    
    
    if isempty(addData{i})
        % Add the set.
        numSets = numSets + 1; %???? Why would we add 1
        %Just because nothing is an observation????
        
        
    else  % Sign the data.
        
        % Add the sets.
        numSets = numSets + length(addData{i});
        
        % Sign the data.
        %
        %
        %   This little bit of code is fairly confusing
        %
        %   true   - if not defined yet and the observation is signed
        %   false  - if any observations are false
        %
        %   Why would these not be internally consistent????
        if ~isempty(addData{i}.isSigned) && ~isnan(addData{i}.isSigned)
            if isempty(isSigned) %i.e. hasn't been defined yet ..
                isSigned = addData{i}.isSigned;
            elseif ~addData{i}.isSigned
                isSigned = false;
            end
        end
    end
end

% Is there any data?
% -> could be empty if any signed entries are empty or NaN
%
%????? - When would we have null data like this??
if isempty(isSigned)
    data = seg_worm.util.nanHistogram(numSets);
    data = rmfield(data, {'sets', 'allData', 'PDF', 'bins', 'resolution', 'isZeroBin'});
    
    %i.e. keep
    %    .isSigned
    %and .data (This contains stats)
    
    return;
end

% Initialize the statistics.
% Note: this must match the field order in worm2histogram.
data.data = [];
data.isSigned = isSigned;

% Combine the data.
data.data.samples    = [];
data.data.mean.all   = [];
data.data.stdDev.all = [];
for i = 1:length(addData)
    if isempty(addData{i})
        data.data.samples    = cat(1, data.data.samples, 0);
        data.data.mean.all   = cat(1, data.data.mean.all, NaN);
        data.data.stdDev.all = cat(1, data.data.stdDev.all, NaN);
    else
        data.data.samples    = cat(1, data.data.samples,    addData{i}.data.samples);
        data.data.mean.all   = cat(1, data.data.mean.all,   addData{i}.data.mean.all);
        data.data.stdDev.all = cat(1, data.data.stdDev.all, addData{i}.data.stdDev.all);
    end
end

if isSigned
    data.data.mean.abs   = [];
    data.data.stdDev.abs = [];
    data.data.mean.pos   = [];
    data.data.stdDev.pos = [];
    data.data.mean.neg   = [];
    data.data.stdDev.neg = [];
    for i = 1:length(addData)
        if isempty(addData{i})
            data.data.mean.abs   = cat(1, data.data.mean.abs,   NaN);
            data.data.stdDev.abs = cat(1, data.data.stdDev.abs, NaN);
            data.data.mean.pos   = cat(1, data.data.mean.pos,   NaN);
            data.data.stdDev.pos = cat(1, data.data.stdDev.pos, NaN);
            data.data.mean.neg   = cat(1, data.data.mean.neg,   NaN);
            data.data.stdDev.neg = cat(1, data.data.stdDev.neg, NaN);
        else
            data.data.mean.abs   = cat(1, data.data.mean.abs,   addData{i}.data.mean.abs);
            data.data.stdDev.abs = cat(1, data.data.stdDev.abs, addData{i}.data.stdDev.abs);
            data.data.mean.pos   = cat(1, data.data.mean.pos,   addData{i}.data.mean.pos);
            data.data.stdDev.pos = cat(1, data.data.stdDev.pos, addData{i}.data.stdDev.pos);
            data.data.mean.neg   = cat(1, data.data.mean.neg,   addData{i}.data.mean.neg);
            data.data.stdDev.neg = cat(1, data.data.stdDev.neg, addData{i}.data.stdDev.neg);
        end
    end
end

end



%% Combine event data statistics.
function data = h__addEventStatistics(wormFiles, wormName, field)

% Initialize the combined statistics.
addData = h__loadWormFiles(wormFiles, wormName, field);
data    = [];
if isempty(addData)
    data.data = NaN;
    return;
end

% Combine the data.
data.data = [];
for i = 1:length(addData)
    if isempty(addData{i})
        data.data = cat(1, data.data, NaN);
    else
        data.data = cat(1, data.data, addData{i}.data);
    end
end
end
