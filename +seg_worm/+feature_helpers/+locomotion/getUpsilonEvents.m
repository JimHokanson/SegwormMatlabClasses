function upsilon_events = getUpsilonEvents(upsilon_frames,midbody_distance,FPS)
%
%
%   upsilon_events = seg_worm.feature_helpers.locomotion.getUpsilonEvents(upsilon_frames)
%

upsilon_events = seg_worm.feature_helpers.locomotion.getTurnEventsFromSignedFrames(upsilon_frames,midbody_distance,FPS);

end