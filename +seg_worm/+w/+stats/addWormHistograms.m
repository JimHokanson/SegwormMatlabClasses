function addWormHistograms(hist_output_filepath, wormFiles, varargin)
%addWormHistograms  Combine worm histograms.
%
%   addWormHistograms(hist_output_filepath, wormFiles, *controlFiles, *add_old_controls, *VERBOSE)
%
%   seg_worm.w.stats.addWormHistograms
%
%
%   This function is meant to take histogram files originally created
%   using  'seg_worm.w.stats.worm2histogram' and merges them. Each of the
%   old files contains a histogram (consisting of one or more sets/videos).
%   The old files may also contain control histograms (which also may have
%   one or more sets of videos).
%
%
%   Inputs:
%       filename     - the file name for the histograms;
%                      the file includes:
%
%                      wormInfo    = the worm information
%                      worm        = the worm histograms
%                      controlInfo = the control information (if it exists)
%                      worm        = the control histograms (if they exist)
%
%       wormFiles    - the histogram files to use for the worm(s)
%       controlFiles - the histogram files to use for the control(s);
%                      if empty, the worm has no new control
%       add_old_controls - are we adding the old controls?
%                      the default is yes (true)
%       isVerbose    - verbose mode display the progress;
%                      the default is no (false)
%
%   See also:
%   seg_worm.w.stats.worm2histogram     %Original function that created wormFiles
%      
%
%   HISTOGRAM
%   WORM2CSV
%   WORMDISPLAYINFO,
%   WORMDATAINFO
%
%
% © Medical Research Council 2012
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
add_old_controls = true;
if length(varargin) > 1
    add_old_controls = varargin{2};
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
if exist(hist_output_filepath, 'file')
    delete(hist_output_filepath);
end

% Save the worm information.
if isVerbose
    disp('Combining "wormInfo" ...');
end
newWormInfo = [];
for i = 1:length(wormFiles)
    h = load(wormFiles{i}, 'wormInfo');
    wormInfo    = h.wormInfo;
    newWormInfo = cat(1, newWormInfo, wormInfo);
end
wormInfo = newWormInfo; %#ok<NASGU>
save(hist_output_filepath, 'wormInfo', '-v7.3');
clear wormInfo;

% Collect the new control information.
if isVerbose
    disp('Combining "controlInfo" ...');
end
newControlInfo = [];
if ~isempty(controlFiles)
    for i = 1:length(controlFiles)
        h = load(controlFiles{i}, 'wormInfo');
        wormInfo = h.wormInfo;
        newControlInfo =  cat(1, newControlInfo, wormInfo);
    end
end

% Collect the old control information.
if add_old_controls
    for i = 1:length(wormFiles)
        controlInfo = who('-FILE', wormFiles{i}, 'controlInfo');
        if ~isempty(controlInfo)
            h              = load(wormFiles{i}, 'controlInfo');
            controlInfo    = h.controlInfo;
            newControlInfo =  cat(1, newControlInfo, controlInfo);
        end
    end
end

% Save the control information.
if ~isempty(newControlInfo)
    controlInfo = newControlInfo; %#ok<NASGU>
    save(hist_output_filepath, 'controlInfo', '-append', '-v7.3');
end
clear controlInfo;

% Initialize the data information.
dataInfo = seg_worm.feature.roots();

% Save the worm histograms.
h__saveHistogram(hist_output_filepath, wormFiles, dataInfo, 'worm', 'worm', isVerbose);

% Are we adding the old controls?
%
%   ??? What does this mean???
%
if add_old_controls
    
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

% Save the control histograms.
if ~isempty(controlFiles)
    h__saveHistogram(hist_output_filepath, controlFiles, dataInfo, controlNames,'control', isVerbose);
end
end



%% Load worm data from files.
function data = h__loadWormFiles(filenames, wormName, field)

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
    data = cellfun(@(x) seg_worm.util.loadStructField(x, wormName{1}, field), ...
        filenames, 'UniformOutput', false);
end
end



%% Save the worm histograms.
function h__saveHistogram(filename, wormFiles, dataInfo, loadName, ...
    saveName, isVerbose)

% Initialize the locomotion modes.
motionNames = { ...
    'forward', ...
    'paused', ...
    'backward'};

% Combine the histograms.
for i = 1:length(dataInfo)
    field = dataInfo(i).field;
    if isVerbose
        disp(['Combining "' field '" ...']);
    end
    switch dataInfo(i).type
        
        % Combine simple histograms.
        case 's'
            data = h__addHistograms(wormFiles, loadName, field);
            eval([saveName '.' field '=data;']);
            
        % Combine motion histograms.
        case 'm'
            data = h__addMotionHistograms(wormFiles, loadName, field, motionNames);
            eval([saveName '.' field '=data;']);
            
        % Combine event histograms.
        case 'e'
            
            % Combine the event data.
            subFields = dataInfo(i).subFields.summary;
            for j = 1:length(subFields);
                subField = [field '.' subFields{j}];
                data = h__addEventData(wormFiles, loadName, subField);
                eval([saveName '.' subField '=data;']);
            end
            
            % Combine the event histograms.
            subFields = dataInfo(i).subFields.data;
            for j = 1:length(subFields);
                subField = [field '.' subFields{j}];
                data = h__addHistograms(wormFiles, loadName, subField);
                eval([saveName '.' subField '=data;']);
            end
    end
end

% Save the histograms.
save(filename, saveName, '-append', '-v7.3');
end



%% Combine histograms.
function data = h__addHistograms(wormFiles, wormName, field)
addData = h__loadWormFiles(wormFiles, wormName, field);
data(length(addData{1})).histogram = [];
for i = 1:length(addData{1})
    subData = cellfun(@(x) x(i).histogram, addData, 'UniformOutput', false);
    data(i).histogram = h__addHistogramData(subData);
end
end



%% Combine motion histograms.
function data = h__addMotionHistograms(wormFiles, wormName, field, motionNames)

% Get the data.
addData = h__loadWormFiles(wormFiles, wormName, field);

% Combine the histograms.
data(length(addData{1})).histogram = [];
for i = 1:length(addData{1})
    
    % Combine the data histograms.
    subData = cellfun(@(x) x(i).histogram, addData, 'UniformOutput', false);
    data(i).histogram = h__addHistogramData(subData);
    
    % Combine the motion histograms.
    for j = 1:length(motionNames)
        subData = cellfun(@(x) x(i).(motionNames{j}).histogram, addData,'UniformOutput', false);
        data(i).(motionNames{j}).histogram = h__addHistogramData(subData);
    end
end
end



%% Combine histograms.
function data = h__addHistogramData(addData)

% Initialize the histogram information and check for consistency.
data       = [];
resolution = [];
isZeroBin  = [];
isSigned   = [];
numSets    = 0;
for i = 1:length(addData)
    
    % Add the set(s).
    if isempty(addData{i})
        numSets = numSets + 1;
    elseif isnan(addData{i}.resolution)
        numSets = numSets + addData{i}.sets.samples;
        isSigned = addData{i}.isSigned;
    else
        numSets = numSets + addData{i}.sets.samples;
        
        % Initialize the histogram information.
        if isempty(resolution)
            resolution = addData{i}.resolution;
            isZeroBin  = addData{i}.isZeroBin;
            isSigned   = addData{i}.isSigned;
            
        % Check the histograms for consistency.
        elseif resolution ~= addData{i}.resolution || isZeroBin ~= addData{i}.isZeroBin
            warning('addWormHistograms:UnequalBins', ...
                'Two or more histograms have inconsistent bins');
            data = seg_worm.util.nanHistogram(isSigned, numSets);
            return;

        % Sign the data.
        elseif addData{i}.isSigned
            isSigned = true;
        end
    end
end

% Is there any data?
if isempty(resolution)
    data = seg_worm.util.nanHistogram(isSigned, numSets);
    return;
end

% Combine the data counts.
numCounts = 0;
counts    = cell(length(addData), 1);
bins      = cell(length(addData), 1);
for i = 1:length(addData)
    
    % No data.
    if isempty(addData{i})
        numCounts = numCounts + 1;
        counts{numCounts} = 0;
        bins{numCounts} = NaN;
        
    % Combine the data.
    else
        for j = 1:size(addData{i}.data.counts, 1)
            numCounts = numCounts + 1;
            counts{numCounts} = addData{i}.data.counts(j,:);
            bins{numCounts} = addData{i}.bins;
        end
    end
end

% Normalize the data and bins.
[counts,bins] = normBinData(counts, bins, resolution);
if isempty(bins) || any(isnan(bins))
    data = [];
    return;
end

% Initialize the histogram.
% Note: this must match the field order in worm2histogram.
data.sets    = [];
data.data    = [];
data.allData = [];
data.PDF     = [];
data.bins    = [];
data.resolution = [];
data.isZeroBin  = [];
data.isSigned   = [];

% Combine all the data.
data.allData.counts   = nansum(counts, 1);
data.allData.samples  = 0;
data.allData.mean.all = 0;
for i = 1:length(addData)
    if ~isempty(addData{i}) && addData{i}.allData.samples > 0
        addSamples = addData{i}.allData.samples;
        data.allData.samples = data.allData.samples + addSamples;
        data.allData.mean.all = data.allData.mean.all + ...
            addData{i}.allData.mean.all * addSamples;
    end
end
if data.allData.samples > 0
    data.allData.mean.all = data.allData.mean.all / data.allData.samples;
else
    data.allData.mean.all = NaN;
end
data.allData.stdDev.all = NaN; % unrecoverable

% Combine all the signed data.
if isSigned
    if data.allData.samples > 0
        data.allData.mean.abs = 0;
        for i = 1:length(addData)
            if ~isempty(addData{i}) && addData{i}.allData.samples > 0
                data.allData.mean.abs = data.allData.mean.abs + ...
                    addData{i}.allData.mean.abs * ...
                    addData{i}.allData.samples;
            end
        end
        data.allData.mean.abs = data.allData.mean.abs / ...
            data.allData.samples;
    else
        data.allData.mean.abs = NaN;
    end
    data.allData.stdDev.abs = NaN; % unrecoverable
    data.allData.mean.pos   = NaN; % unrecoverable
    data.allData.stdDev.pos = NaN; % unrecoverable
    data.allData.mean.neg   = NaN; % unrecoverable
    data.allData.stdDev.neg = NaN; % unrecoverable
end

%This is the part that really matters ...

% Combine the data.
%==========================================================================
data.data.counts     = counts;
data.data.samples    = [];
data.data.mean.all   = [];
data.data.stdDev.all = [];
for i = 1:length(addData)
    if isempty(addData{i})
        data.data.samples    = cat(1, data.data.samples, 0);
        data.data.mean.all   = cat(1, data.data.mean.all, NaN);
        data.data.stdDev.all = cat(1, data.data.stdDev.all, NaN);
    else
        data.data.samples = ...
            cat(1, data.data.samples, addData{i}.data.samples);
        data.data.mean.all = ...
            cat(1, data.data.mean.all, addData{i}.data.mean.all);
        data.data.stdDev.all = ...
            cat(1, data.data.stdDev.all, addData{i}.data.stdDev.all);
    end
end

% Combine the signed data.
if isSigned
    data.data.mean.abs   = [];
    data.data.stdDev.abs = [];
    data.data.mean.pos   = [];
    data.data.stdDev.pos = [];
    data.data.mean.neg   = [];
    data.data.stdDev.neg = [];
    for i = 1:length(addData)
        if isempty(addData{i})
            data.data.mean.abs   = cat(1, data.data.mean.abs, NaN);
            data.data.stdDev.abs = cat(1, data.data.stdDev.abs, NaN);
            data.data.mean.pos   = cat(1, data.data.mean.pos, NaN);
            data.data.stdDev.pos = cat(1, data.data.stdDev.pos, NaN);
            data.data.mean.neg   = cat(1, data.data.mean.neg, NaN);
            data.data.stdDev.neg = cat(1, data.data.stdDev.neg, NaN);
        else
            data.data.mean.abs = ...
                cat(1, data.data.mean.abs, addData{i}.data.mean.abs);
            data.data.stdDev.abs = ...
                cat(1, data.data.stdDev.abs, addData{i}.data.stdDev.abs);
            data.data.mean.pos = ...
                cat(1, data.data.mean.pos, addData{i}.data.mean.pos);
            data.data.stdDev.pos = ...
                cat(1, data.data.stdDev.pos, addData{i}.data.stdDev.pos);
            data.data.mean.neg = ...
                cat(1, data.data.mean.neg, addData{i}.data.mean.neg);
            data.data.stdDev.neg = ...
                cat(1, data.data.stdDev.neg, addData{i}.data.stdDev.neg);
        end
    end
end

%==========================================================================


% Combine the sets.
data.sets.samples = sum(~isnan(data.data.mean.all));
data.sets.mean.all = nanmean(data.data.mean.all);
data.sets.stdDev.all = nanstd(data.data.mean.all);
if isSigned
    data.sets.mean.abs = nanmean(data.data.mean.abs);
    data.sets.stdDev.abs = nanstd(data.data.mean.abs);
    data.sets.mean.pos = nanmean(data.data.mean.pos);
    data.sets.stdDev.pos = nanstd(data.data.mean.pos);
    data.sets.mean.neg = nanmean(data.data.mean.neg);
    data.sets.stdDev.neg = nanstd(data.data.mean.neg);
end

% Combine the PDFs.
data.bins = bins;
data.PDF = zeros(1, length(bins));
numSets = 0;
for i = 1:length(data.data.samples)
    if data.data.samples(i) > 0
        data.PDF = data.PDF + data.data.counts(i,:) ./ data.data.samples(i);
        numSets = numSets + 1;
    end
end
if numSets > 0
    data.PDF = data.PDF ./ numSets;
end

% Set the histogram information.
data.resolution = resolution;
data.isZeroBin  = isZeroBin;
data.isSigned   = isSigned;
end



%% Normalize binned data.
function [normData,normBins] = normBinData(data, bins, resolution)

% Initialize the normalized data and bins.
normData = [];
normBins = [];
if isempty(data)
    return;
end

% Compute the normalized bins.
minBin = min(cellfun(@(x) x(1), bins));
if isnan(minBin)
    return;
end
maxBin   = max(cellfun(@(x) x(end), bins));
numBins  = round((maxBin - minBin) / resolution) + 1;
normBins = linspace(minBin, maxBin, numBins);

% Normalize the data.
normData = zeros(length(data), numBins);
for i = 1:length(data)
    
    % No data.
    if isnan(bins{i})
        continue;
    end
    
    % Copy the data.
    startI = round((bins{i}(1) - minBin) / resolution) + 1;
    endI = round((bins{i}(end) - minBin) / resolution) + 1;
    normData(i, startI:endI) = data{i};
end
end



%% Combine event data.
function data = h__addEventData(wormFiles, wormName, field)

% Initialize the combined histogram.
addData = h__loadWormFiles(wormFiles, wormName, field);
data = [];
if isempty(addData)
    data.samples = 0;
    data.data = 0;
    data.mean = 0;
    data.stdDev = 0;
    return;
end

% Combine the data.
data.data = [];
for i = 1:length(addData)
    if isempty(addData{i})
        data.data = cat(1, data.data, 0);
    else
        addData{i}.data(isnan(addData{i}.data)) = 0;
        data.data = cat(1, data.data, addData{i}.data);
    end
end

% Combine the sets.
data.samples = length(data.data);
data.mean = nanmean(data.data);
data.stdDev = nanstd(data.data);
end
