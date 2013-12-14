function coils = getCoils(frame_codes,midbody_distance,FPS)
%
%   seg_worm.feature_helpers.posture.getCoils
%
%
%   Old Name:
%   - wormTouchFrames
%
%
%   Nature Methods Description
%   =======================================================================
%   Coils. 
%   ----------------------------------------------
%   Worm coiling (touching) events are found by scanning the video frame
%   annotations. During segmentation, every frame that cannot be segmented
%   is annotated with a cause for failure. Two of these annotations reflect
%   coiling events. 
%
%   First, if we find fewer than two sharp ends on the contour (reflecting
%   the head and tail) then the head and/or tail are obscured in a coiling
%   event.
%
%   Second, if the length between the head and tail on one side of the
%   contour is more than double that of the other side, the worm has either
%   assumed an omega bend or is crossed like a wreath.
%
%   Empirically, less than 1/5 of a second is a very fast touch and not
%   usually reflective of coiling. Therefore, when a period of unsegmented
%   video frames exceeds 1/5 of a second, and either of the coiling
%   annotations are found, we label the event coiling.

INTER_DATA_NAME = 'interDistance';
DATA_NAME = [];

coiled_frames = seg_worm.feature_helpers.posture.wormTouchFrames(frame_codes, FPS);

coiled_events = seg_worm.feature.event(coiled_frames,FPS,midbody_distance,DATA_NAME,INTER_DATA_NAME);
coils = coiled_events.getFeatureStruct;

end

function touchFrames = h__getWormTouchFrames(frameCodes, fps)

COIL_FRAME_THRESHOLD = round(1/5 * fps);

%TODO: Push this up ...
COIL_START_CODES = [105 106];
FRAME_SEGMENTED  = 1; %Go back 1 frame, this is the end of the coil ...

%Algorithm: Whenever a new start is found, find the first segmented frame, 
%that's the end.

%Add on a frame to allow closing a coil at the end ...
coil_start_mask = [frameCodes == COIL_START_CODES(1) | frameCodes == COIL_START_CODES(2) false];

%NOTE: These are not guaranteed ends, just possible ends ...
end_coil_mask   = [frameCodes == FRAME_SEGMENTED true];

in_coil = false;
coil_frame_start = 0;

n_coils = 0;

n_frames = length(frameCodes) + 1;

for iFrame = 1:n_frames
    if in_coil
        if end_coil_mask(iFrame)
            
            n_coil_frames = iFrame - coil_frame_start;
            if n_coil_frames >= COIL_FRAME_THRESHOLD
                n_coils = n_coils + 1;
                
                %NOTE: The output is 0 based (for now ...) :/
                touchFrames(n_coils).start = coil_frame_start - 1; %#ok<AGROW>
                touchFrames(n_coils).end   = iFrame - 2; %#ok<AGROW> %:/
            end
            in_coil = false;
        end
    elseif coil_start_mask(iFrame)
        in_coil = true;
        coil_frame_start = iFrame;
    end
end


end