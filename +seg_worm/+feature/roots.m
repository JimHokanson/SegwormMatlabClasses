function info = roots()
%roots  Get root definitions for calculating all features.
%
%   info = seg_worm.feature.roots()
%
%       This provides a base set from which all features can be calculated.
%
%   Old Name: wormDataInfo
%
%   seg_worm.w.stats.wormStatsInfo <- This seems like a wrapper around
%   this function ...
%
%   seg_worm.feature.displayInfo
%
%   Called by:
%   filterWormHist
%   plotWormStats
%   worm2TOCOntologyTIF
%   worm2detailsPDF2
%   worm2detailsTIF
%   worm2summaryTIF
%   addWormHistograms
%   addWormStats
%   worm2NormTimes
%   seg_worm.w.stats.worm2histogram
%   seg_worm.w.stats.worm2stats
%   wormDataInfoIndices
%   wormStats2Matrix2
%   wormStatsInfo
%   
%
%   Output 
%   -----------------------------------------------------------------------
%   info : (structure array)
% 
%      field        = the feature's path
%                     This can be nested:
%                     'posture.bends.head.mean'
%
%      subFields    = the feature's subfields (for complex features);
%                     for an event type:
% 
%                     summary = the summary subfields
%                     data    = the data subfields
%                     sign    = a field describing the sign
%                               (true/false = +/-);
%                               if empty, the event is not signed
%
%               EXAMPLE:
%               turnEventFields.summary = {'frequency', 'timeRatio'};
%               turnEventFields.data    = {'time', 'interTime', 'interDistance'};
%               turnEventFields.sign    = 'isVentral';
%
% 
%      type         = the feature's type where:
% 
%                     s = a simple data array
%                     m = a data array which can be further divided
%                         into forward/paused/backward motion
%                     e = an event data array
% 
%      category     = the feature's category where:
% 
%                     m = morphology
%                     s = posture (shape)
%                     l = locomotion
%                     p = path
%
%      isTimeSeries = is the data a time series?
%
%   See also:
%   seg_worm.feature.displayInfo
%   seg_worm.w.stats.worm2histogram, 
%   ADDWORMHISTOGRAMS, 
%   seg_worm.w.stats.worm2stats,
%   WORM2NORMTIMES, 
%   WORM2CSV
%
%
% © Medical Research Council 2012
% You will not remove any copyright or other notices from the Software; 
% you must reproduce all copyright notices and other proprietary 
% notices on any copies of the Software.

%Q&A
%--------------------------------------------------------------------------
%{

1. Why do we care of something is a time-series
   - only the following functions care:
        - plotWormStats.m
        - worm2NormTimes.m


%}


% Initialize the data types.
simpleDataType = 's';
motionDataType = 'm';
eventDataType  = 'e';

% Initialize the data categories.
morphologyType = 'm';
postureType    = 's';
locomotionType = 'l';
pathType       = 'p';

% Initialize the duration data field.

%Simple Data Type
%----------------------------------------------------------------
% 'path.duration.worm'
% 'path.duration.head'
% 'path.duration.midbody'
% 'path.duration.tail'};
DURATION_FIELD = {'times'};  %?? This doesn't follow the same format ...

%EVENT NOTES
%----------------------------------------------------------------
%'summary' - stat fields
%
%   summary - provides some overall summary of the event data
%
% - data
% - length(data)
% - nanmean(data)
% - nanstd(data)
%
%'data' - 
%
%   Very similar to the other fields, just they are present
%   in the 'frames' struct
%
%   example:
%
%   posture.coils.frames.time
%   posture.coils.frames.interTime
%   posture.coils.frames.interDistance
%

%interTime - time to next event, last event value NaN (I think this is right)
%interDistance - ????

% Initialize the coil event data fields.
%-----------------------------------------------------------------
% posture.coils
coilEventFields.summary = {'frequency', 'timeRatio'};
coilEventFields.data    = {'time', 'interTime', 'interDistance'};
coilEventFields.sign    = [];

% Initialize the turn event data fields.
%-----------------------------------------------------------------
% locomotion.turns.omegas
% locomotion.turns.upsilons
turnEventFields.summary = {'frequency', 'timeRatio'};
turnEventFields.data    = {'time', 'interTime', 'interDistance'};
turnEventFields.sign    = 'isVentral';

% Initialize the motion event data fields.
%-----------------------------------------------------------------
% locomotion.motion.forward
% locomotion.motion.paused
% locomotion.motion.backward
motionEventFields.summary = {'frequency', 'ratio.time', 'ratio.distance'};
motionEventFields.data    = {'time', 'distance', 'interTime', 'interDistance'};
motionEventFields.sign    = [];

% Initialize the features.
info = struct('field',{},'subFields',{},'type',{},'category',{},'isTimeSeries',{});

morphology_fields = {
    'morphology.length'
    'morphology.width.head'
    'morphology.width.midbody'
    'morphology.width.tail'
    'morphology.area'
    'morphology.areaPerLength'
    'morphology.widthPerLength'};

for iField = 1:length(morphology_fields)
   info(end+1).field      = morphology_fields{iField}; %#ok<AGROW>
   info(end).type         = motionDataType;
   info(end).category     = morphologyType;
   info(end).isTimeSeries = true; 
end

posture_fields = {
    'posture.bends.head.mean'
    'posture.bends.neck.mean'
    'posture.bends.midbody.mean'
    'posture.bends.hips.mean'
    'posture.bends.tail.mean'
    'posture.bends.head.stdDev'
    'posture.bends.neck.stdDev'
    'posture.bends.midbody.stdDev'
    'posture.bends.hips.stdDev'
    'posture.bends.tail.stdDev'
    'posture.amplitude.max'
    'posture.amplitude.ratio'
    'posture.wavelength.primary'
    'posture.wavelength.secondary'
    'posture.tracklength'
    'posture.eccentricity'
    'posture.kinks'
    'posture.directions.tail2head'
    'posture.directions.head'
    'posture.directions.tail'
    'posture.eigenProjection'};

for iField = 1:length(posture_fields)
   info(end+1).field      = posture_fields{iField}; %#ok<AGROW>
   info(end).type         = motionDataType;
   info(end).category     = postureType;
   info(end).isTimeSeries = true; 
end

info(end+1).field      = 'posture.coils';
info(end).subFields    = coilEventFields;
info(end).type         = eventDataType;
info(end).category     = postureType;
info(end).isTimeSeries = true;


locomotion_motion_fields = {
    'locomotion.velocity.headTip.speed'
    'locomotion.velocity.head.speed'
    'locomotion.velocity.midbody.speed'
    'locomotion.velocity.tail.speed'
    'locomotion.velocity.tailTip.speed'
    'locomotion.velocity.headTip.direction'
    'locomotion.velocity.head.direction'
    'locomotion.velocity.midbody.direction'
    'locomotion.velocity.tail.direction'
    'locomotion.velocity.tailTip.direction'
    'locomotion.bends.foraging.amplitude'
    'locomotion.bends.foraging.angleSpeed'
    'locomotion.bends.head.amplitude'
    'locomotion.bends.midbody.amplitude'
    'locomotion.bends.tail.amplitude'
    'locomotion.bends.head.frequency'
    'locomotion.bends.midbody.frequency'
    'locomotion.bends.tail.frequency'};


for iField = 1:length(locomotion_motion_fields)
   info(end+1).field      = locomotion_motion_fields{iField}; %#ok<AGROW>
   info(end).type         = motionDataType;
   info(end).category     = locomotionType;
   info(end).isTimeSeries = true; 
end

%LOCOMOTION MOTION --------------------------------------------------------
locomotion_event_fields = {
    'locomotion.motion.forward'
    'locomotion.motion.paused'
    'locomotion.motion.backward'};

for iField = 1:length(locomotion_event_fields)
   info(end+1).field      = locomotion_event_fields{iField}; %#ok<AGROW>
   info(end).subFields    = motionEventFields;
   info(end).type         = eventDataType;
   info(end).category     = locomotionType;
   info(end).isTimeSeries = true; 
end

%LOCOMOTION EVENTS --------------------------------------------------------
locomotion_turn_fields = {
    'locomotion.turns.omegas'
    'locomotion.turns.upsilons'};

for iField = 1:length(locomotion_turn_fields)
   info(end+1).field      = locomotion_turn_fields{iField}; %#ok<AGROW>
   info(end).subFields    = turnEventFields;
   info(end).type         = eventDataType;
   info(end).category     = locomotionType;
   info(end).isTimeSeries = true; 
end

%PATH MOTION --------------------------------------------------------------
path_motion_fields = {
    'path.range'
    'path.curvature'};

for iField = 1:length(path_motion_fields)
   info(end+1).field      = path_motion_fields{iField}; %#ok<AGROW>
   info(end).type         = motionDataType;
   info(end).category     = pathType;
   info(end).isTimeSeries = true; 
end

%PATH SIMPLE -------------------------------------------------------------
path_simple_fields = {
    'path.duration.worm'
    'path.duration.head'
    'path.duration.midbody'
    'path.duration.tail'};

for iField = 1:length(path_simple_fields)
   info(end+1).field      = path_simple_fields{iField}; %#ok<AGROW>
   info(end).subFields    = DURATION_FIELD;
   info(end).type         = simpleDataType;
   info(end).category     = pathType;
   info(end).isTimeSeries = false; 
end