function getUpsilonEvents(obj,upsilon_frames,midbody_distance,FPS)
%
%   seg_worm.features.locomotion.getUpsilonEvents
%

obj.turns.upsilons = obj.getTurnEventsFromSignedFrames(upsilon_frames,midbody_distance,FPS);

end