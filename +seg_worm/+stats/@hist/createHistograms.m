function createHistograms(feature_file_path)
%worm2histogram   Convert worm features to their histogram.
%
%   seg_worm.stats.hist.createHistograms
%
%   OLD_NAME: seg_worm.w.stats.worm2histogram 
%
%
%   JAH: Currently working on this file ...
%
%   Code snippet for running in:
%   seg_worm.stats
%   
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

WORM_NAME = 'h.worm.';

% Initialize the worm data information.
%--------------------------------------------------------------------------
h = load(feature_file_path);


%This code should eventually be translated into a different function
%
%   seg_worm.stats.hist.getSpecs() or something like that ...
%
%We will ultimately be working with the specs ...
histInfo  = seg_worm.feature.displayInfo();
data_info = seg_worm.feature.roots();

n_fields = length(data_info);

temp_specs = cell(1,n_fields);
for iField = 1:n_fields
   cur_data_info = data_info(iField);
   if cur_data_info.type ~= 'e'
      cur_hist_info = eval(['histInfo.' cur_data_info.field]);
      temp_specs{iField} = seg_worm.stats.hist_specs(cur_data_info,cur_hist_info);
   end
end

specs = [temp_specs{:}];
%--------------------------------------------------------------------------



%Motion Expanded Features
%--------------------------------------------------------------------------
%Always contains at least 4 values, for during forward, backward, paused,
%and all motions ...
%
%Might contain 


n_specs = length(specs);

%feature_file_path

motion_modes = h.worm.locomotion.motion.mode;

m.is_forward  = motion_modes == 1;
m.is_paused   = motion_modes == 0;
m.is_backward = motion_modes == -1;

% Combine the histograms.
for iField = 1:n_specs
    
    cur_spec       = specs(iField);
    cur_field      = specs(iField).field;    
    cur_sub_fields = specs(iField).sub_fields;
        
    %wormName - either 'worm' or 'control'
    
    %This is where the magic happens
    %----------------------------------------------------------------------
    switch specs(iField).type
        
        % Compute the simple histogram.
        case 's'
            
            %data = eval('WORM_NAME
            
            %
            %
            %   data - field and sub fields
            
            %Examples
            %-----------------------------------------
            continue
            keyboard
            
            temp = seg_worm.stats.hist(data,resolution,is_zero_bin,is_signed);
            
            
            data = h__data2histogram(wormFiles, cur_field, cur_sub_fields, histInfo);
            eval([wormName '.' cur_field '=data;']);
            
            %field -> path.duration.worm
            %
            %data =>
            %    histogram: [1x1 struct]
            
            
        % Compute the motion histogram.
        case 'm'

            data = eval([WORM_NAME cur_spec.field]);
            
            %TODO: Add index check ...
            
            obj = seg_worm.stats.hist(data,cur_spec.resolution,cur_spec.is_zero_bin,cur_spec.is_signed);

            
            %{
            
            temp = seg_worm.stats.hist(data,resolution,is_zero_bin,is_signed);

            
            
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
            
            keyboard
            
            
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
            
            
            %}
    end
end

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


%Where is this done ??????

%{

Im = find([r.type] == 'm');

I = find([r.type] == 'm' & arrayfun(@(x) ~isempty(x.subFields),r));

%We never have subFields for movement events

%}




data = h__loadWormFiles(wormFiles, dataField);


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

function h__actualCreateHistogram()

function histData = histogram(data, varargin)
%histogram  Compute the data's histogram.
%   
%   histData = seg_worm.util.histogram(data, *resolution, *isZeroBin, *isSigned, *verbose)
%
%   Inputs:
%       data       - the data for the histogram(s) (a cell array of observed sets)
%       resolution - the data resolution;
%                    if empty, the resolution is estimated from the square
%                    root of the data's length (a very popular method)
%       isZeroBin  - is there a bin center at 0?
%                    if empty, data with values greater and less than 0 is
%                    centered at 0
%       isSigned   - is the data signed (+/-)?
%                    if empty, the data is tested to see whether its signed
%       verbose    - verbose mode shows the results in a figure
%
%   Outputs:
%
%           Example:
%       worm.morphology.length.histogram
%                       sets: [1x1 struct]
%                       data: [1x1 struct]
%                    allData: [1x1 struct]
%                        PDF: [1x418 double]
%                       bins: [1x418 double]
%                 resolution: 1
%                  isZeroBin: 0
%                   isSigned: 0
%   -----------------------------------------------------------------------
%       histData - the histogram(s) data. A struct with fields:
%                  Note: the absolute, positive, and negative value
%                  statistics are only available when the data is signed
%                  (see wormDisplayInfo).
%
%                  PDF        = the PDF (probability density of each bin)
%                  bins       = the center of each bin
%                  resolution = the data resolution
%                  isZeroBin  = is there a bin center at 0
%                  isSigned   = is the data signed (+/-)?
%                  --------------------------------------------------------
%                  sets       = a struct with fields: (Statistics computed
%                           on the stats in 'data', i.e mean(data.mean)
%                     samples    = the number of set samples
%                     mean.all   = the mean of the sets
%                     stdDev.all = the standard deviation of the sets
%                     mean.abs   = the mean of the set magnitudes
%                     stdDev.abs = the deviation of the set magnitudes
%                     mean.pos   = the mean of the positive set values
%                     stdDev.pos = the deviation of the positive set values
%                     mean.neg   = the mean of the negative set values
%                     stdDev.neg = the deviation of the negative set values
%
%                  --------------------------------------------------------
%                  NOTE: This field ('data') is considered to be the 
%                  "Stats" that is referred to in other functions.
%                  --------------------------------------------------------
%                  data       = a struct with fields: (Statistics computed
%                    on a video by video basis, i.e. [mean(set1) mean(set2) mean(set3)])
%
%                     counts     = the count per set per bin (sets x bins)
%                     samples    = the number of data samples per set
%                     mean.all   = the mean of the data per set
%                     stdDev.all = the standard deviation the data per set
%                     mean.abs   = the mean of the data magnitudes per set
%                     stdDev.abs = the deviation of the data magnitudes per set
%                     mean.pos   = the mean of the positive data per set
%                     stdDev.pos = the deviation of the positive data per set
%                     mean.neg   = the mean of the negative data per set
%                     stdDev.neg = the deviation of the negative data per set
%
%                  --------------------------------------------------------
%                  allData    = a struct with fields: (Statistics computed
%                       on all data i.e. mean(all_data))
%
%                     counts     = the count per bin
%                     samples    = the total number of all data samples
%                     mean.all   = the mean of all the data
%                     stdDev.all = the standard deviation of all the data
%                     mean.abs   = the mean of all the data magnitudes
%                     stdDev.abs = the deviation of all the data magnitudes
%                     mean.pos   = the mean of all the positive data
%                     stdDev.pos = the deviation of all the positive data
%                     mean.neg   = the mean of all the negative data
%                     stdDev.neg = the deviation of all the negative data
%
%
%   IMPROVEMENTS:
%   1) It would be nice to know what histogram data we were binning 
%
%
%
% © Medical Research Council 2012
% You will not remove any copyright or other notices from the Software; 
% you must reproduce all copyright notices and other proprietary 
% notices on any copies of the Software.

% Organize the data.
if ~iscell(data)
    data = {data};
end

% Determine the data resolution.
resolution = [];
if ~isempty(varargin)
    resolution = varargin{1};
end

% Is there a bin center at 0?
isZeroBin = [];
if length(varargin) > 1
    isZeroBin = varargin{2};
end

% Is the data signed (+/-)?
isSigned = [];
if length(varargin) > 2
    isSigned = varargin{3};
end

% Are we displaying the results?
verbose = false;
if length(varargin) > 3
    verbose = varargin{4};
end

% Combine the data into a single array of observations.
if size(data, 1) ~= size(data{1}, 1) && size(data, 2) ~= size(data{1}, 2)
    data = data';
end
allData = cell2mat(data);
allData = allData(:);

% Remove empty data.
samples = cellfun(@(x) sum(~isnan(x)), data);

%JAH QUESTION: Is commenting this out an appropriate fix or not?
%This can cause a problem ...
%------------------------------------------------
%This means that although we asked for n sets,
%we might not get n sets of data out ...
%
% % % % keepI = find(samples > 0);
% % % % if length(data) ~= length(keepI)
% % % %     keyboard
% % % %     data = data(keepI);
% % % % end
%
%   Problem seen with a backwards event ...
%   worm2StatsInfo/worm2data
%   

% Remove infinite data measurements.
for i = 1:length(data)
    data{i}(data{i} == -inf) = NaN;
    data{i}(data{i} == inf)  = NaN;
end

% Is there any data?
if isempty(data)
    warning('histogram:NoData', 'There is no data to bin');
    histData = seg_worm.util.nanHistogram(isSigned);
    return;
end
empty = cellfun(@(x) isempty(x) || all(isnan(x)), data);
if all(empty)
    warning('histogram:NoData', 'There is no data to bin');
    histData = seg_worm.util.nanHistogram(isSigned, length(empty));
    return;
end

% Compute the data range.
minDataLength = min(cellfun('length', data));
minDataLength = max(minDataLength,1);
minData = min(allData);
maxData = max(allData);

% Compute the data resolution.
if isempty(resolution)
    resolution = (maxData - minData) / sqrt(minDataLength);
end

% Compute the bin center.
if isempty(isZeroBin)
    isZeroBin = sign(minData) ~= sign(maxData);
end

% Compute the data sign.
if isempty(isSigned)
    isSigned = any(allData < 0);
end

% Compute the padding.
if minData < 0
    minPad = resolution - abs(rem(minData, resolution));
else
    minPad = abs(rem(minData, resolution));
end
if maxData < 0
    maxPad = abs(rem(maxData, resolution));
else
    maxPad = resolution - abs(rem(maxData, resolution));
end

% Translate the bins by half the resolution to create a zero bin.
% Note: we compute just enough padding to capture the data.
halfResolution = resolution / 2;
if isZeroBin
    if minPad > halfResolution
        minPad = minPad - halfResolution;
    else
        minPad = minPad + halfResolution;
    end
    if maxPad > halfResolution
        maxPad = maxPad - halfResolution;
    else
        maxPad = maxPad + halfResolution;
    end
end

% Compute the edge range.
minEdge = minData - minPad;
maxEdge = maxData + maxPad;

% Compute the bins and their edges.
% Note: histc fills all bins with edges(k) <= data < edges(k + 1).
% The last bin is filled with data == edges(end).
% Therefore, we keep the last bin empty and throw it away to preserve
% equal bin sizes. For this reason the bin centers are spaced for
% their final composition (what they will look like after tossing away
% the empty last bin).
numBins = round((maxEdge - minEdge) / resolution);
bins    = linspace(minEdge + halfResolution, maxEdge - halfResolution, numBins);
edges   = bins - halfResolution;
edges(end + 1) = edges(end) + resolution;

% Fix the zero bin.
% Note: IEEE floating point issues may shift us just off zero.
if isZeroBin
    [zeroBin, zeroI] = min(abs(bins));
    if zeroBin < halfResolution / 2
        bins(zeroI) = 0;
    end
end
    
% Compute the histogram counts for all the data.
counts(1,:) = histc(allData, edges);
if length(edges) > 1
    
    % Add the very last bin.
    if counts(1,end) > 0
        bins(end + 1) = bins(end) + resolution;
        
    % Strip off the empty last bin.
    else
        counts(end) = [];
    end
    
    % Strip off the empty first bin.
    if counts(1,1) == 0
        bins(1) = [];
        counts(1) = [];
    end
end

% Compute the counts for the data set.
allCounts(1,:) = counts;
if length(data) == 1
    pdfs(1,:) = counts ./ sum(counts);
    
% Compute the normalized histogram for the data sets.
else
    counts = zeros(length(data), length(edges));
    pdfs   = zeros(length(data), length(edges));
    for i = 1:length(data)
        if ~empty(i)
            counts(i,:) = histc(data{i}, edges);
            pdfs(i,:)   = counts(i,:) ./ sum(counts(i,:));
        end
    end
    pdfs = mean(pdfs, 1);
    
    % Strip off the empty bin.
    if length(edges) > 1
        
        % Add a bin (this should never happen).
        if any(counts(:,end) > 0)
            bins(end + 1) = bins(end) + resolution;
            warning('histogram:LastBinNotEmpty', ...
                ['The last bin in the histogram is not empty like it ' ...
                'should be. Please contact the programmer)']);
            
        % Strip off the empty bin.
        else
            counts(:,end) = [];
            pdfs(:,end) = [];
        end
    end
end

% Compute the means of the data sets.
means   = cellfun(@nanmean, data);
stdDevs = cellfun(@nanstd, data);
if isSigned
    
    % Compute the absolute value statisitics.
    absData    = cellfun(@(x) abs(x), data, 'UniformOutput', false);
    absMeans   = cellfun(@(x) nanmean(x), absData);
    absStdDevs = cellfun(@(x) nanstd(x), absData);
    
    % Compute the positive value statisitics.
    posData    = cellfun(@(x) x(x > 0), data, 'UniformOutput', false);
    posMeans   = cellfun(@(x) nanmean(x), posData);
    posStdDevs = cellfun(@(x) nanstd(x), posData);
    
    % Compute the negative value statisitics.
    negData    = cellfun(@(x) x(x < 0), data, 'UniformOutput', false);
    negMeans   = cellfun(@(x) nanmean(x), negData);
    negStdDevs = cellfun(@(x) nanstd(x), negData);
end

% Organize the histogram set data.
histData.sets.samples    = length(data);
histData.sets.mean.all   = nanmean(means);
histData.sets.stdDev.all = nanstd(means);
if isSigned
    
    % Compute the absolute value statisitics.
    histData.sets.mean.abs   = nanmean(absMeans);
    histData.sets.stdDev.abs = nanstd(absMeans);
    
    % Compute the positive value statisitics.
    histData.sets.mean.pos   = nanmean(posMeans);
    histData.sets.stdDev.pos = nanstd(posMeans);
    
    % Compute the negative value statisitics.
    histData.sets.mean.neg   = nanmean(negMeans);
    histData.sets.stdDev.neg = nanstd(negMeans);
end

% Organize the histogram data sets.
histData.data.counts          = counts;
histData.data.samples(:,1)    = samples; %Yikes, I think is is to force a colum vector ...
histData.data.mean.all(:,1)   = means;
histData.data.stdDev.all(:,1) = stdDevs;
if isSigned
    
    % Compute the absolute value statisitics.
    histData.data.mean.abs(:,1)   = absMeans;
    histData.data.stdDev.abs(:,1) = absStdDevs;
    
    % Compute the positive value statisitics.
    histData.data.mean.pos(:,1)   = posMeans;
    histData.data.stdDev.pos(:,1) = posStdDevs;
    
    % Compute the negative value statisitics.
    histData.data.mean.neg(:,1)   = negMeans;
    histData.data.stdDev.neg(:,1) = negStdDevs;
end

% Organize the histogram total data.
histData.allData.counts     = allCounts;
histData.allData.samples    = sum(samples);
histData.allData.mean.all   = nanmean(allData);
histData.allData.stdDev.all = nanstd(allData);
if isSigned
    
    % Compute the absolute value statisitics.
    absAllData = abs(allData);
    histData.allData.mean.abs   = nanmean(absAllData);
    histData.allData.stdDev.abs = nanstd(absAllData);
    
    % Compute the positive value statisitics.
    posAllData = allData(allData > 0);
    histData.allData.mean.pos   = nanmean(posAllData);
    histData.allData.stdDev.pos = nanstd(posAllData);
    
    % Compute the negative value statisitics.
    negAllData = allData(allData < 0);
    histData.allData.mean.neg   = nanmean(negAllData);
    histData.allData.stdDev.neg = nanstd(negAllData);
end

% Organize the histogram.
histData.PDF  = pdfs;
histData.bins = bins;
histData.resolution = resolution;
histData.isZeroBin  = isZeroBin;
histData.isSigned   = isSigned;

% Show the results.
if verbose
    
    % Create a figure.
    h = figure;
    set(h, 'units', 'normalized', 'position', [0 0 1 1]);
    hold on;
    
    % Plot the histogram.
    if length(data) > 1
        subplot(1, 2, 1);
    end
    histColor = [0 0 0];
    statColor = [1 0 0];
    plotHistogram(histData, 'Histogram', 'Value', 'Data', '', ...
        histColor, statColor, 2);
    
    % Plot the contributing histograms.
    if length(data) > 1
        
        % Construct the contributing histograms.
        histDatas(1:length(data)) = histData;
        for i = 1:length(data)
            
            % Set the statistics.
            histDatas(i).sets.samples    = 1;
            histDatas(i).sets.mean       = histData.data.mean(i);
            histDatas(i).sets.stdDev     = 0;
            histDatas(i).data.counts     = histData.data.counts(i,:);
            histDatas(i).data.samples    = histData.data.samples(i);
            histDatas(i).data.mean       = histData.data.mean(i);
            histDatas(i).data.stdDev     = histData.data.stdDev(i);
            histDatas(i).allData.samples = histData.data.samples(i);
            histDatas(i).allData.mean    = histData.data.mean(i);
            histDatas(i).allData.stdDev  = histData.data.stdDev(i);
            
            % Set the PDF.
            counts = histDatas(i).data.counts;
            histDatas(i).PDF = counts ./ sum(counts);
        end
        
        % Plot the contributing histograms.
        subplot(1, 2, 2);
        dataNames(1:length(data)) = {'Data'};
        colors = lines(length(data));
        seg_worm.util.plotHistogram(histDatas, 'Histogram', 'Value', dataNames, colors, colors);
    end
end
end



end
