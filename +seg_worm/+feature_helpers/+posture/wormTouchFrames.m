function touchFrames = wormTouchFrames(frameCodes, fps)
%WORMTOUCHFRAMES Find the frames where the worm touches itself (i.e., coils).
%
%   Used for detecting coil events
%
%
%   seg_worm.feature_helpers.posture.wormTouchFrames    
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
% See also SEGWORM, WORMFRAMEANNOTATION
%
%
% � Medical Research Council 2012
% You will not remove any copyright or other notices from the Software; 
% you must reproduce all copyright notices and other proprietary 
% notices on any copies of the Software.

% Initialize the frame codes.
codes          = helper__wormFrameAnnotation();
successful_segmentation  = codes(1).id; %Successful segmentation

stage_motion   = codes(2).id;
dropped_video_frame = codes(3).id;
%noWormCode    = codes(4).id;
%boundaryCode  = codes(5).id;
%smallWormCode = codes(6).id;
tooFewEndsCode = codes(8).id;
doubleLengthSideCode = codes(9).id;


%Descriptions
%   too_few_ends 



%Rules
%-----------------------------------------------
%1) 



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


