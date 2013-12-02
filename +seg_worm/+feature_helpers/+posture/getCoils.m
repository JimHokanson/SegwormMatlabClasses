function coils = getCoils(n_frames,frame_codes)
%
%   coils = seg_worm.feature_helpers.posture.getCoils(n_frames,frame_codes)
%
%                   UNFINISHED   UNFINISHED    UNFINISHED
%
%
%   Nature Methods Description
%   =======================================================================
%   Coils. Worm coiling (touching) events are found by scanning the video
%   frame annotations. During segmentation, every frame that cannot be
%   segmented is annotated with a cause for failure. Two of these
%   annotations reflect coiling events. First, if we find fewer than two
%   sharp ends on the contour (reflecting the head and tail) then the head
%   and/or tail are obscured in a coiling event. Second, if the length
%   between the head and tail on one side of the contour is more than
%   double that of the other side, the worm has either assumed an omega
%   bend or is crossed like a wreath. Empirically, less than 1/5 of a
%   second is a very fast touch and not usually reflective of coiling.
%   Therefore, when a period of unsegmented video frames exceeds 1/5 of a
%   second, and either of the coiling annotations are found, we label the
%   event coiling.


FPS = 20;


%TODO: Midbody speed needs to be passed in ..., or just the distance ...

%{

speed = velocity.midbody.speed;

% Compute the distance. This avoids both segmentation noise and having to
% re-interpolate the distance when frames are missing.

distance = abs(speed / fps);

%}
distance = rand(1,n_frames); %This is temporary ...



%NOTE: Output is zero based currently :/
% TODO: This should be changed, use 1-based for Matlab ...
coiled_frames = seg_worm.feature_helpers.posture.wormTouchFrames(frame_codes, FPS);

%TODO: I might want to create a static method of the class so I can merge these two lines ...
coiled_events = seg_worm.feature.event(coiled_frames,FPS,distance,[],'interDistance');
coils = coiled_events.getFeatureStruct;