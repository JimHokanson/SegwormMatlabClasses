function all_hists = createHistograms(feature_file_paths)
%worm2histogram   Convert worm features to their histogram.
%
%   seg_worm.stats.hist.createHistograms
%
%   OLD_NAME: 
%   - seg_worm.w.stats.worm2histogram 
%
%
%   JAH: Currently working on this file ...
%
%   Code snippet for running in:
%   seg_worm.stats.testing_code
%   
%
%   See also:
%   ADDWORMHISTOGRAMS, 
%   seg_worm.util.histogram  %This is a key function
%   WORM2CSV, 
%   WORMDISPLAYINFO,
%   WORMDATAINFO
%
%   See Also:
%   seg_worm.features.wormDataInfo

if ischar(feature_file_paths)
    feature_file_paths = {feature_file_paths};
end

n_files = length(feature_file_paths);
hist_cell_array = cell(1,n_files);

for iFile  = 1:n_files
    hist_cell_array{iFile} = h__getHistsFromFeaturePath(feature_file_paths{iFile});
end

all_hists = seg_worm.stats.hist.mergeObjects(hist_cell_array);

end

function hist_objs = h__getHistsFromFeaturePath(feature_file_path)

% Initialize the worm data information.
%--------------------------------------------------------------------------
h = load(feature_file_path);

%Movement histograms - DONE
m_hists = h_computeMHists(h,seg_worm.stats.movement_specs.getSpecs);

%Simple histograms - DONE
s_hists = h_computeSHists(h,seg_worm.stats.simple_specs.getSpecs);

%Event histograms - DONE
e_hists = h_computeEHists(h,seg_worm.stats.event_specs.getSpecs);

hist_objs = [m_hists s_hists e_hists];

end


%% Convert data to a histogram.

function e_hists = h_computeEHists(h,specs)

n_specs    = length(specs);
temp_hists = cell(1,n_specs);

for iSpec = 1:n_specs
   
   cur_specs = specs(iSpec);

   %NOTE: Because we are doing structure array indexing, we need to capture
   %multiple outputs using [], otherwise we will only get the first value
   %...
   cur_data  = eval(['h.' cur_specs.feature_field]);
    
   if ~isempty(cur_data) && ~isempty(cur_specs.sub_field)
      %This will go from:
      %   frames (structure array) 
      %to:
      %   frames.time
      %for example.
      %
      %It is also used for event.ratio.time and event.ratio.distance
      %going from:
      %ratio (structure or [])
      %to:
      %   ratio.time
      %   ratio.distance
      %
      %
      cur_data = [cur_data.(cur_specs.sub_field)];
   end
   
   temp_hists{iSpec} = seg_worm.stats.hist(cur_data,cur_specs,'event','all','all');
   
end

e_hists = [temp_hists{:}];


end

function m_hists = h_computeMHists(h,specs)
%
%
%   We still need to populate the meta data
%   
%   m_spec needs to inherit from spec, which will
%   have its valued copied by hist ...
%

MAX_NUM_HIST_OBJECTS = 1000;

motion_modes = h.worm.locomotion.motion.mode;

m.is_forward  = motion_modes == 1;
m.is_paused   = motion_modes == 0;
m.is_backward = motion_modes == -1;

n_frames = length(motion_modes);

indices_by_motion = {1:n_frames ...
    find(motion_modes == 1) ...
    find(motion_modes == 0) ...
    find(motion_modes == -1)};

motion_types = {'all' 'forward' 'paused' 'backward'};
data_types   = {'all' 'absolute' 'positive' 'negative'};

all_hist_objects = cell(1,MAX_NUM_HIST_OBJECTS);
hist_count = 0;

n_specs = length(specs);

for iSpec = 1:n_specs
   
   cur_specs = specs(iSpec);
   cur_data  = eval(['h.' cur_specs.feature_field]);
    
   if ~isnan(cur_specs.index)
      %This is basically for eigenprojections
      %JAH: I really don't like the orientation: [Dim x n_frames]
      cur_data = cur_data(cur_specs.index,:);
   end
   
   for iMotion = 1:4 
      cur_motion_type = motion_types{iMotion};
       
      hist_count = hist_count + 1;
      temp_data = cur_data(indices_by_motion{iMotion});
      all_hist_objects{hist_count} = seg_worm.stats.hist(temp_data,cur_specs,'motion',cur_motion_type,data_types{1});
      if cur_specs.is_signed
         all_hist_objects{hist_count+1} = seg_worm.stats.hist(abs(temp_data),          cur_specs,'motion',cur_motion_type,data_types{2});
         all_hist_objects{hist_count+2} = seg_worm.stats.hist(temp_data(temp_data > 0),cur_specs,'motion',cur_motion_type,data_types{3});
         all_hist_objects{hist_count+3} = seg_worm.stats.hist(temp_data(temp_data < 0),cur_specs,'motion',cur_motion_type,data_types{4});
         hist_count = hist_count + 3;     
      end
   end
end

m_hists = [all_hist_objects{1:hist_count}];

end

function s_hists = h_computeSHists(h,specs)

n_specs = length(specs);

temp_hists = cell(1,n_specs);

for iSpec = 1:n_specs
   
   cur_specs = specs(iSpec);
   cur_data  = eval(['h.' cur_specs.feature_field]);
    
   temp_hists{iSpec} = seg_worm.stats.hist(cur_data,cur_specs,'simple','all','all');
   
end

s_hists = [temp_hists{:}];


end

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

