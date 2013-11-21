function touchFrames = wormTouchFrames(frameCodes, fps)
%
%   touchFrames = seg_worm.feature_helpers.posture.wormTouchFrames(frameCodes,fps)
%
%   frameCodes: [1 n]

% Coils. Worm coiling (touching) events are found by scanning the video frame
% annotations. During segmentation, every frame that cannot be segmented is
% annotated with a cause for failure. Two of these annotations reflect coiling events.
% First, if we find fewer than two sharp ends on the contour (reflecting the head and
% tail) then the head and/or tail are obscured in a coiling event. Second, if the length
% between the head and tail on one side of the contour is more than double that of
% the other side, the worm has either assumed an omega bend or is crossed like a
% wreath. Empirically, less than 1/5 of a second is a very fast touch and not usually
% reflective of coiling. Therefore, when a period of unsegmented video frames
% exceeds 1/5 of a second, and either of the coiling annotations are found, we label
% the event coiling.


COIL_FRAME_THRESHOLD = round(1/5 * fps);

%TODO: Push this up ...
COIL_START_CODES = [105 106];
FRAME_SEGMENTED  = 1; %Go back 1 frame, this is the end of the coil ...


% tic
% touchFrames2 = oldwormTouchFrames(frameCodes, fps);
% toc

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

%{

function touchFrames = oldwormTouchFrames(frameCodes, fps)


%WORMTOUCHFRAMES Find the frames where the worm touches itself (i.e., coils).
%
%   Used for detecting coil events
%
%   JAH: I'm currently rewriting this code ...
%
%   TODO: This should be renamed to getWormCoilFrames
%
%
%   touchFrames = seg_worm.feature_helpers.posture.wormTouchFrames(frameCodes, fps)
%
%   EXAMPLES
%   =======================================================================
%   fps = 5; %Make the threshold really low for testing cases ...
%   frameCodes = [2 104 105 3 1]
%   start 2, end 2
%
%   frameCodes = [105 105 105 105 1]
%   start 0, end 3
%
%   frameCodes = [2 105 105 1]
%   start 1, end 1
%
%   frameCodes = [105 2 106 1]
%   start 0, end 2
%
%   SURPRISES: output is 0 based ...
%
%
%   TOUCHFRAMES = WORMTOUCHFRAMES(FRAMECODES, FPS)
%
%   Input:
%       frameCodes - the frame codes annotating the worm segmentation per
%                    video frame
%       fps        - the video's frames/second
%
%   Output:
%       touchFrames - a struct containing the frames where the worm
%                     touches itself; the fields are:
%
%                     start = the starting frame wherein the worm
%                             initiates the touch
%                     end   = the ending frame wherein the worm
%                             terminates the touch
%
%
% © Medical Research Council 2012
% You will not remove any copyright or other notices from the Software;
% you must reproduce all copyright notices and other proprietary
% notices on any copies of the Software.

% Initialize the frame codes.
codes          = helper__wormFrameAnnotation();
successful_segmentation  = codes(1).id; %1 %Successful segmentation

stage_motion   = codes(2).id;           %2
dropped_video_frame = codes(3).id;      %3
%noWormCode    = codes(4).id;
%boundaryCode  = codes(5).id;
%smallWormCode = codes(6).id;

%These initiate a touch
%--------------------------------------------------------------------------
tooFewEndsCode       = codes(8).id;           %105
doubleLengthSideCode = codes(9).id;     %106


%Summary:
%--------------------------------------------
%1) A good segmentation wipes the slate clean
%
%2) 105 and 106 initiate a coiling
%
%Any other issues are counted if following a coiling



% Compute the touch duration threshold.
touchThr = round(1/5 * fps);

% Find the touch frames.
touchFrames = struct( ...
    'start',    [], ...
    'end',      []);
touchFramesI = 0;
numTouches   = 0;
numVideoErrs = 0;
numOtherErrs = 0;
for i = 1:length(frameCodes)
    
    switch frameCodes(i)
        
        % Do we have a potential touch frame?
        case {tooFewEndsCode, doubleLengthSideCode}
            
            % Absorb any intervening frames.
            numTouches   = numTouches + numVideoErrs + numOtherErrs;
            numVideoErrs = 0;
            numOtherErrs = 0;
            numTouches   = numTouches + 1;
            
            % Do we have an intervening video issue?
        case {stage_motion, dropped_video_frame}
            if numTouches > 0
                numVideoErrs = numVideoErrs + 1;
            end
            
            % Do we have a potential non-touch frame?
        case successful_segmentation %, noWormCode, boundaryCode, smallWormCode}
            
            % Do we have enough potential touch frames?
            if numTouches + numVideoErrs + numOtherErrs >= touchThr
                numFrames    = numTouches + numVideoErrs + numOtherErrs;
                numTouches   = numTouches + numOtherErrs;
                touchFramesI = touchFramesI + 1;
                touchFrames(touchFramesI).start = i - numFrames - 1;
                touchFrames(touchFramesI).end   = touchFrames(touchFramesI).start + numTouches - 1;
            end
            
            % Intialize the frame counts.
            numTouches   = 0;
            numVideoErrs = 0;
            numOtherErrs = 0;
            
            % Do we have an intervening segmentation error?
        otherwise
            
            % Absorb any video issues.
            if numTouches > 0
                numOtherErrs = numOtherErrs + numVideoErrs;
                numVideoErrs = 0;
                numOtherErrs = numOtherErrs + 1;
            end
    end
end

% At the end of the video, do we have enough potential touch frames?
%NOTE: This is the same as the code in the loop ...
if numTouches + numVideoErrs + numOtherErrs >= touchThr
    numFrames    = numTouches + numVideoErrs + numOtherErrs;
    numTouches   = numTouches + numOtherErrs;
    touchFramesI = touchFramesI + 1;
    touchFrames(touchFramesI).start = i - numFrames - 1;
    touchFrames(touchFramesI).end   = touchFrames(touchFramesI).start + numTouches - 1;
end

% Did we find any touch frames?
if isempty(touchFrames(1).start)
    touchFrames = [];
end

end

function annotation = helper__wormFrameAnnotation()

annotation = struct( ...
    ... % Annotation IDs.
    'id', ...
    { 1, ...
    2, ...
    3, ...
    101, ...
    102, ...
    103, ...
    104, ...
    105, ...
    106, ...
    107, ...
    108, ...
    109, ...
    110, ...
    111, ...
    1001}, ...
    ... % Annotation function signatures.
    'function', ...
    { 'segWorm:Success', ...
    'findStageMovement:StageMovement', ...
    'segWorm:DroppedFrame', ...
    'segWorm:NoWormFound', ...
    'segWorm:ContourTouchesBoundary', ...
    'segWorm:ContourTooSmall', ...
    'segWorm:TooManyEnds', ...
    'segWorm:TooFewEnds', ...
    'segWorm:DoubleLengthSide', ...
    'segWorm:DoubleHeadWidth', ...
    'segWorm:DoubleTailWidth', ...
    'segWorm:SmallTail', ...
    'segWorm:SmallHead', ...
    'segWorm:SmallHeadTail', ...
    'normWorms:TooShort'}, ...
    ... % Annotation function signatures.
    'message', ...
    { 'The worm was successfully segmented.', ...
    'The video frame contains stage motion.', ...
    'The video frame was dropped.', ...
    'No worm was found in the video frame.', ...
    'The worm contour touches the image boundary.', ...
    'The worm contour is too small.', ...
    ['The worm has 3 or more low-frequency sampled convexities ' ...
    'sharper than 90 degrees (possible head/tail points).'], ...
    ['The worm contour has less than 2 high-frequency sampled ' ...
    'convexities sharper than 60 degrees (the head and tail). ' ...
    'Therefore, the worm is coiled or obscured and cannot be ' ...
    'segmented.'], ...
    ['The worm length, from head to tail, is more than twice as ' ...
    'large on one side than it is on the other. Therefore, the worm ' ...
    'is coiled or obscured and cannot be segmented.'], ...
    ['The worm more than doubles its width from the end of its head. ' ...
    'Therefore, the worm is coiled, laid an egg, and/or is ' ...
    'significantly obscured and cannot be segmented'], ...
    ['The worm more than doubles its width from the end of its tail. ' ...
    'Therefore, the worm is coiled, laid an egg, and/or is ' ...
    'significantly obscured and cannot be segmented.'], ...
    ['The worm tail is less than half the size of its head. ' ...
    'Therefore, the worm is significantly obscured and cannot be ' ...
    'segmented.'], ...
    ['The worm head is less than half the size of its tail. ' ...
    'Therefore, the worm is significantly obscured and cannot be ' ...
    'segmented.'], ...
    ['The worm head and tail are less than 1/4 the size of its ' ...
    'remaining body. Therefore, the worm is significantly obscured ' ...
    'and cannot be segmented.'], ...
    'The worm is shorter than the sampling points requested.'});
end

%}

