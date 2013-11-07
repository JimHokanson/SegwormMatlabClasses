classdef frame_errors < sl.obj.handle_light
    %
    %   Class:
    %   seg_worm.parsing.frame_errors
    %
    %   There is a file that ends with the extension:
    %
    %   _failedFrames
    %
    %   It contains the variable 'failedFrames'
    %
    %   Column 1: frame number
    %   Column 2: warning id
    %
    %
    %   NOTE: This class 
    
    properties
       frame_numbers
       id
    end
    
    methods
    end
    
end

    { 1 1       'segWorm:Success'                   'The worm was successfully segmented.' ...
      2 2       'findStageMovement:StageMovement'   'The video frame contains stage motion.' ...
      3 3       'segWorm:DroppedFrame'              'The video frame was dropped.'
      4 101     'segWorm:NoWormFound'               'No worm was found in the video frame.'
      5 102     'segWorm:ContourTouchesBoundary'    'The worm contour touches the image boundary.'
      6 103     'segWorm:ContourTooSmall'           'The worm contour is too small.'
      7 104     'segWorm:TooManyEnds'               'The worm has 3 or more low-frequency sampled convexities sharper than 90 degrees (possible head/tail points).'
      8 105, 	'segWorm:TooFewEnds'                'The worm contour has less than 2 high-frequency sampled convexities sharper than 60 degrees (the head and tail). Therefore, the worm is coiled or obscured and cannot be segmented.', ...
      9 106,    'segWorm:DoubleLengthSide'          'The worm length, from head to tail, is more than twice as large on one side than it is on the other. Therefore, the worm is coiled or obscured and cannot be segmented.', ...
      10 107,   'segWorm:DoubleHeadWidth'           'The worm more than doubles its width from the end of its head. Therefore, the worm is coiled, laid an egg, and/or is significantly obscured and cannot be segmented'], ...
      11 108,   'segWorm:DoubleTailWidth'           'The worm more than doubles its width from the end of its tail. Therefore, the worm is coiled, laid an egg, and/or is significantly obscured and cannot be segmented.', ...
      12 109,   'segWorm:SmallTail'                 'The worm tail is less than half the size of its head. Therefore, the worm is significantly obscured and cannot be segmented.'
      13 110, 	'segWorm:SmallHead'                 'The worm head is less than half the size of its tail. Therefore, the worm is significantly obscured and cannot be segmented.'
      14 111,   'segWorm:SmallHeadTail'             'The worm head and tail are less than 1/4 the size of its remaining body. Therefore, the worm is significantly obscured and cannot be segmented.'
      15 1001   'normWorms:TooShort'}






      'The worm is shorter than the sampling points requested.'});