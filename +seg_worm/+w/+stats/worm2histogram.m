function worm2histogram(hist_output_file_path, wormFiles, varargin)
%worm2histogram   Convert worm features to their histogram.
%
%   seg_worm.w.stats.worm2histogram(hist_output_file_path, wormFiles, *CONTROLFILES, *VERBOSE, *PROGFUNC, *PROGSTATE)
%
%
%   Called by:
%   worm2Histogram_GUI
%   This also seems like it could be a top level file
%
%
%   hist_output_file_path:
%       saved in this file will be the following:
%           .worm     : output of seg_worm.util.histogram
%           .wormInfo : structure array of the 'info' entry from each feature
%           MIGHT ALSO CONTAIN - these are the same as above, just for 
%           feature files that have been designated "controls"
%           .control     :
%           .controlInfo :
%
%   Inputs:
%
%           ????? - is this going to be an output????
%       filename     - the file name for the histograms;
%                      the file includes:
%
%                      wormInfo    = the worm information
%                      worm        = the worm histograms
%                      controlInfo = the control information (if it exists)
%                      worm        = the control histograms (if they exist)
%
%
%
%       wormFiles    - the feature files to use for the worm
%                       
%                   These contain the 'info' and 'worm' structures.
%
%   
%       controlFiles - the feature files to use for the control;
%                      if empty, the worm has no control
%       isVerbose    - verbose mode display the progress;
%                      the default is no (false)
%       progFunc     - a function to update on the progress
%       progState    - a state for the progress function
%
%       Note: the progress function signature should be
%
%       FUNCSTATE = PROGFUNC(PERCENT, MSG, FUNCSTATE)
%
%       Arguments:
%          funcState = a progress function state
%          percent   = the progress percent (0 to 100%)
%          msg       = a message on our progress (to display)
%
% See also:
%   ADDWORMHISTOGRAMS, 
%   seg_worm.util.histogram  %This is a key function
%   WORM2CSV, 
%   WORMDISPLAYINFO,
%   WORMDATAINFO
%
%   See Also:
%   seg_worm.features.wormDataInfo
%   
%
%
% © Medical Research Council 2012
% You will not remove any copyright or other notices from the Software; 
% you must reproduce all copyright notices and other proprietary 
% notices on any copies of the Software.

% Do we have a control?
controlFiles = [];
if ~isempty(varargin)
    controlFiles = varargin{1};
end

% Are we in verbose mode?
isVerbose = false;
if length(varargin) > 1
    isVerbose = varargin{2};
end

% Initialize the progress function.
progFunc  = [];
progState = [];
if length(varargin) > 2
    progFunc = varargin{3};
    if length(varargin) > 3 && ~isempty(progFunc)
        progState = varargin{4};
    end
end

% Organize the worm files.
if ~iscell(wormFiles)
    wormFiles = {wormFiles};
end
if ~isempty(controlFiles) && ~iscell(controlFiles)
    controlFiles = {controlFiles};
end

% Delete the file if it already exists.
if exist(hist_output_file_path, 'file')
    delete(hist_output_file_path);
end

% Initialize the worm data information.
%--------------------------------------------------------------------------
histInfo = seg_worm.feature.displayInfo();
dataInfo = seg_worm.feature.roots();

% Initilize the progress.
progCount = 0;
progSize  = length(dataInfo) + 2;

% Save the worm information.
if isVerbose
    disp('Saving "wormInfo" ...');
end
[progState, progCount] = h__progress('Saving "wormInfo" ...',progFunc, progState, progCount, progSize);

wormInfo = cellfun(@(x) load(x, 'info'), wormFiles); %Yikes!!!
%This loads the 'info' entry for each worm file
wormInfo = [wormInfo.info];

n_frames_per_video = arrayfun(@(x) x.video.length.frames, wormInfo);
if isempty(controlFiles)
    save(hist_output_file_path, 'wormInfo', '-v7.3');
else
    if isVerbose
        disp('Saving "controlInfo" ...');
    end
    controlInfo = cellfun(@(x) load(x, 'info'), controlFiles);
    controlInfo = [controlInfo.info];
    if size(controlInfo, 1) < size(controlInfo, 2)
        controlInfo = controlInfo';
    end
    controlFrames = arrayfun(@(x) x.video.length.frames, controlInfo);
    save(hist_output_file_path, 'wormInfo', 'controlInfo', '-v7.3');
end

% Free memory.
clear('wormInfo', 'controlInfo');

% Save the worm histograms.
h__saveHistogram(hist_output_file_path, wormFiles, n_frames_per_video, histInfo, dataInfo, ...
    'worm', isVerbose, progFunc, progState, progCount, progSize);

% Save the control histograms.
if ~isempty(controlFiles)
    h__saveHistogram(hist_output_file_path, controlFiles, controlFrames, histInfo, ...
        dataInfo, 'control', isVerbose, ...
        progFunc, progState, progCount, progSize);
end
end



%% Load worm data from files.
function data = h__loadWormFiles(filenames, field)
    %From the worm structure in the file, load the desired field
    %
    %   NOTE: The files contain two structures:
    %   worm
    %   info
    data = cellfun(@(x) seg_worm.util.loadStructField(x, 'worm', field), filenames, 'un', 0);
end



%% Save the worm histograms.
function h__saveHistogram(filename, wormFiles, frames, histInfo, dataInfo, ...
    wormName, isVerbose, progFunc, progState, progCount, progSize)

% Determine the locomotion modes.
motionModes = h__loadWormFiles(wormFiles, 'locomotion.motion.mode');
motionNames = {'forward'; 'paused'; 'backward'};
motionEvents = { ...
    cellfun(@(x) x == 1, motionModes, 'un', 0), ...
    cellfun(@(x) x == 0, motionModes, 'un', 0), ...
    cellfun(@(x) x == -1, motionModes, 'un', 0)};

% Check the locomotion modes.
for i = 1:length(motionEvents)
    for j = 1:length(motionEvents{i})
        if all(motionEvents{i}{j} ~= 1)
            warning('worm2histogram:NoMotionData', ...
                ['There are no ' motionNames{i} ' motion frames in "' wormFiles{j} '"']);
        end
    end
end

% Combine the histograms.
for i = 1:length(dataInfo)
    field = dataInfo(i).field;
    if isVerbose
        disp(['Saving "' field '" ...']);
    end
    [progState, progCount] = h__progress(['Saving "' field '" ...'], ...
        progFunc, progState, progCount, progSize);
    
    
    subFields = dataInfo(i).subFields;
    
    
    %wormName - either 'worm' or 'control'
    
    %This is where the magic happens
    %----------------------------------------------------------------------
    switch dataInfo(i).type
        
        % Compute the simple histogram.
        case 's'
            data = h__data2histogram(wormFiles, field, subFields, histInfo);
            eval([wormName '.' field '=data;']);
            
            %field -> path.duration.worm
            %
            %data =>
            %    histogram: [1x1 struct]
            
            
        % Compute the motion histogram.
        case 'm'
            %????? - where would this be the case ???
            if ~isempty(dataInfo(i).subFields)
                field = [field '.' dataInfo(i).subFields{1}];
            end
            data = h__motion2histograms(wormFiles, field, subFields, histInfo, motionNames, motionEvents);
            eval([wormName '.' field '=data;']);
            
%             if strcmp(field,'posture.wavelength.secondary')
%                keyboard 
%             end
            
            %field -> morphology.length
            %
            %data =>
            %     histogram: [1x1 struct]
            %       forward: [1x1 struct]
            %        paused: [1x1 struct]
            %      backward: [1x1 struct]
            %
            %   data.forward =>
            %     histogram: [1x1 struct]
            
            
        % Compute the event histogram.
        case 'e'
            data = h__event2histograms(wormFiles, frames, field, subFields.summary, subFields.data, subFields.sign, histInfo);
            
            %field -> posture.coils
            %
            %data =>
            %             frequency: [1x1 struct] <- summary field
            %             timeRatio: [1x1 struct] <- summary field
            %                  time: [1x1 struct] <- data field
            %             interTime: [1x1 struct] <- data field
            %         interDistance: [1x1 struct] <- data field
            %
            %   data.frequency => NOTE: summary fields have this form
            %            data: [4x1 double]
            %         samples: 4
            %            mean: 0.0047
            %          stdDev: 0.0033
            %
            %   data.time      =>
            %           histogram: [1x1 struct]
            %   data.time.histogram =>
            %           sets: [1x1 struct]
            %           data: [1x1 struct]
            %        allData: [1x1 struct]
            %            PDF: [1x32 double]
            %           bins: [1x32 double]
            %     resolution: 0.1000
            %      isZeroBin: 0
            %       isSigned: 0 
            
            eval([wormName '.' field '=data;']);
    end
end

% Save the histograms.
save(filename, wormName, '-append', '-v7.3');
end



%% Convert data to a histogram.

function histData = h__data2histogram(wormFiles, field, subFields, histInfo)
%
%   This function is for 'simple' data types.
%
%   See: seg_worm.feature.roots
%
%   example: path.duration.worm
%
%   
%
%
%   Called by: h__saveHistogram
%

% Get the histogram information.
resolution = [];
isZeroBin  = [];
isSigned   = [];
info       = seg_worm.util.getStructField(histInfo, field);
if isempty(info)
    warning('worm2histogram:NoInfo', ...
        ['There is no information for "' field '"']);
else
    resolution = info.resolution;
    isZeroBin  = info.isZeroBin;
    isSigned   = info.isSigned;
end

% Get the data.
dataField = field;
if ~isempty(subFields)
    dataField = [dataField '.' subFields{1}];
end
data = h__loadWormFiles(wormFiles, dataField);

% Check the data.
for i = 1:length(data)
    if isempty(data{i})
        warning('worm2histogram:NoData', ['"' field '" in "' ...
            wormFiles{i} '" contains no data']);
    end
end

% Compute the histogram.
histData.histogram = seg_worm.util.histogram(data, resolution, isZeroBin, isSigned);
end



%% Convert motion data to a set of histograms.
function histData = h__motion2histograms(wormFiles, field, subFields, histInfo, motionNames, motionEvents)
%
%
%   Called by: h__saveHistogram
%   
%
%     wormFiles : 
%     field     : ex. morphology.length
%     subFields : ????
%     histInfo  : See seg_worm.feature.displayInfo
%   
%           Important fields that in histInfo (at 'field'):
%
%                   i.e. histInfo.morphology.length =>
%             resolution: 1
%              isZeroBin: 0
%               isSigned: 0
%                   name: 'Length'
%              shortName: 'Length'
%                   unit: 'Microns'

%This is passsed in ...
% motionNames = {'forward'; 'paused'; 'backward'};
% motionEvents = { ...
%     cellfun(@(x) x == 1, motionModes, 'un', 0), ...
%     cellfun(@(x) x == 0, motionModes, 'un', 0), ...
%     cellfun(@(x) x == -1, motionModes, 'un', 0)};



% Get the histogram information.
resolution = [];
isZeroBin  = [];
isSigned   = [];
info       = seg_worm.util.getStructField(histInfo, field);
if isempty(info)
    warning('worm2histogram:NoInfo', ['There is no information for "' field '"']);
else
    resolution = info.resolution;
    isZeroBin  = info.isZeroBin;
    isSigned   = info.isSigned;
end

% Get the data.
dataField = field;
if ~isempty(subFields)
    dataField = [dataField '.' subFields{1}];
end
data = h__loadWormFiles(wormFiles, dataField);

% Check the data.
for i = 1:length(data)
    if isempty(data{i})
        warning('worm2histogram:NoData', ['"' field '" in "' wormFiles{i} '" contains no data']);
    end
end

% Initialize the histogram data.
histData(size(data{1}, 1)).histogram = [];
for i = 1:size(data{1}, 1)
    for j = 1:length(motionNames)
        histData(i).(motionNames{j}).histogram = [];
    end
end

% Go through the data.
for i = 1:size(data{1}, 1) %Why would this have more than one entry ?????
    %??? - it looks like size(data{1},1) might always be 1 ...
    
    % Compute the data histogram.
    subData = cellfun(@(x) x(i,:), data, 'un', 0);
    histData(i).histogram = seg_worm.util.histogram(subData, resolution, isZeroBin, isSigned);
    
    % Compute the motion histograms.
    for j = 1:length(motionEvents)
        
        
        %motionEvents
        %
        %   {1x4 cell}    {1x4 cell}    {1x4 cell}
        %
        %   4 experiments, indices (3 total), represent
        %   forward, paused, backward
        
        
        % Organize the motion data.
        isData     = false;
        motionData = cell(length(subData),1);
        for k = 1:length(subData) %For each experiment ...
            %motionEvents{j}{k}
            motionData{k} = subData{k}(motionEvents{j}{k}); 
            if ~isempty(motionData{k})
                isData = true;
            end
        end
        
%    motionData => {4 x 1}
%     [1x1750 double]
%     [1x1955 double]
%     [1x3102 double]
%     [1x2635 double]
        
        
        % Compute the motion histogram.
        if isData
            histData(i) = seg_worm.util.setStructField(...
                histData(i), ...
                [motionNames{j} '.histogram'], ...
                seg_worm.util.histogram(motionData, resolution, isZeroBin, isSigned));
        end
    end
end
end



%% Convert event data to a set of histograms.
function histData = h__event2histograms(wormFiles, frames, field, ...
    statFields, histFields, signField, histInfo)
%
%   See Also: seg_worm.feature.roots
%   field       : ex. 'posture.coils'
%   statFields  : ex. {'frequency' 'timeRatio'} aka 'summary'
%   histFields  : ex. {'time' 'interTime' 'interDistance'} aka 'data'
%


% Get the data.
data = h__loadWormFiles(wormFiles, field);

% Remove partial events.
for i = 1:length(data)
    data{i}.frames = seg_worm.feature.removePartialEvents(data{i}.frames, frames(i));
end

% Check the data.
for i = 1:length(data)
    if isempty(data{i})
        warning('worm2histogram:NoEventData', ['"' field ...
            '" in "' wormFiles{i} '" contains no data (excluding ' ...
            'partial events at the start and end of the video)']);
    end
end

% Initialize the histogram data.
histData = [];
for i = 1:length(statFields)
    histData = seg_worm.util.setStructField(histData, statFields{i}, []);
end

% Compute the event statistics. 'summary'
%--------------------------------------------------------------------------
for i = 1:length(statFields)

    % Organize the event data.
    eventData = zeros(length(data),1);
    for j = 1:length(data)
        subData = seg_worm.util.getStructField(data{j}, statFields{i});
        if ~isempty(subData)
            eventData(j) = subData;
        end
    end
    
    % Compute the event statistics.
    %
    %   ??? - length - what is an example where this is more than one per
    %   data set???
    %   
    %   ??? - nanmean, nanstd - doesn't length become innacurate??
    %
    %
    histData = seg_worm.util.setStructField(histData, [statFields{i} '.data'],      eventData);
    histData = seg_worm.util.setStructField(histData, [statFields{i} '.samples'],   length(eventData));
    histData = seg_worm.util.setStructField(histData, [statFields{i} '.mean'],      nanmean(eventData));
    histData = seg_worm.util.setStructField(histData, [statFields{i} '.stdDev'],    nanstd(eventData));
end

% Create the histogram structures.
for i = 1:length(histFields)
    histData = seg_worm.util.setStructField(histData, [histFields{i} '.histogram'], []);
end

% Compute the event histograms.
%--------------------------------------------------------------------------
for i = 1:length(histFields)
    
    %Example: - 'time' 'interTime' 'interDistance'
    
    % Get the histogram information.
    resolution = [];
    isZeroBin  = [];
    isSigned   = [];
    subField   = [field '.' histFields{i}];
    info = seg_worm.util.getStructField(histInfo, subField);
    if isempty(info)
        %??? - when would this be called ...
        warning('worm2histogram:NoInfo', ['There is no information for "' subField '"']);
    else
        resolution = info.resolution;
        isZeroBin  = info.isZeroBin;
        isSigned   = info.isSigned;
    end

    % Organize the event data.
    eventData = cell(length(data),1);
    for j = 1:length(data)
        subData = data{j}.frames; %All the event data is in .frames ...
        if ~isempty(field)
            eventData{j} = seg_worm.util.getStructField(subData, histFields{i});
            
            % Sign the event.
            if ~isempty(signField)
                signs = seg_worm.util.getStructField(subData, signField);
                %??? certain indices are flipped???
                %signs - logical array
                eventData{j}(signs) = -eventData{j}(signs);
            end
        end
    end
    %At this point, we have a cell array where each element
    %is an experiment, and each element contains
    %the values from the specified field
    %
    %   At this point it is similar to other data and is ready
    %   to be put into the histogram function ...
    
    % Compute the event histogram.
    histData = seg_worm.util.setStructField(histData, [histFields{i} '.histogram'], ...
        seg_worm.util.histogram(eventData, resolution, isZeroBin, isSigned));
end
end



%% Progress.
function [progState, progCount] = ...
    h__progress(msg, progFunc, progState, progCount, progSize)

% Are we monitoring progression?
if ~isempty(progFunc)
    
    % Progress.
    progCount = progCount + 1;
    progState = progFunc(progCount / progSize, msg, progState);
end
end
