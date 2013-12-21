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

%Movement histograms - DONE
m_hists = h_computeMHists(h,seg_worm.stats.movement_specs.getSpecs);

%Simple histograms
s_hists = h_computeSHists(h,asdf);







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

function m_hists = h_computeMHists(h,m_specs)
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

n_specs = length(m_specs);

for iSpec = 1:n_specs
   
   cur_specs = m_specs(iSpec);
    
   cur_data = eval(['h.' cur_specs.feature_field]);
    
   if ~isnan(cur_specs.index)
      %This is basically for eigenprojections
      %JAH: I really don't like the orientation: [Dim x n_frames]
      cur_data = cur_data(cur_specs.index,:);
   end
   
   for iMotion = 1:4 
      cur_motion_type = motion_types{iMotion};
       
      hist_count = hist_count + 1;
      temp_data = cur_data(indices_by_motion{iMotion});
      all_hist_objects{hist_count} = seg_worm.stats.hist(temp_data,cur_specs,cur_motion_type,data_types{1});
      if cur_specs.is_signed
         all_hist_objects{hist_count+1} = seg_worm.stats.hist(abs(temp_data),cur_specs,cur_motion_type,data_types{2});
         all_hist_objects{hist_count+2} = seg_worm.stats.hist(temp_data(temp_data > 0),cur_specs,cur_motion_type,data_types{3});
         all_hist_objects{hist_count+3} = seg_worm.stats.hist(temp_data(temp_data < 0),cur_specs,cur_motion_type,data_types{4});
         hist_count = hist_count + 3;     
      end
   end
end

m_hists = [all_hist_objects{1:hist_count}];

end

function s_hists = h_computeSHists(h,asdf)

%Not yet implemented


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

