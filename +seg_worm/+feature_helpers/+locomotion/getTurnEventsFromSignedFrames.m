function turn_events = getTurnEventsFromSignedFrames(signed_frames,midbody_distance,FPS)
%
%
%   turn_event = seg_worm.feature_helpers.locomotion.getTurnEventFromSignedFrames(signed_frames,midbody_distance,FPS)


INTER_DATA_SUM_NAME = 'interDistance';
DATA_SUM_NAME       = '';

ef = seg_worm.feature.event_finder;
ef.include_at_thr = true;

%seg_worm.feature.event_finder.getEvents
frames_dorsal  = ef.getEvents(signed_frames,1,[]);
frames_ventral = ef.getEvents(signed_frames,[],-1);

% Unify the ventral and dorsal turns.
%--------------------------------------------------------------------------
[frames_merged,is_ventral] = frames_ventral.merge(frames_dorsal);

temp = seg_worm.feature.event(frames_merged,FPS,midbody_distance,DATA_SUM_NAME,INTER_DATA_SUM_NAME); 
turn_events = temp.getFeatureStruct;

%Add extra field, isVentral ...
%---------------------------------------------------------------
n_events = length(turn_events.frames);
for iEvent = 1:n_events
   turn_events.frames(iEvent).isVentral = is_ventral(iEvent); 
end


end