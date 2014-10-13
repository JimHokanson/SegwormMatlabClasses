classdef worm < sl.obj.handle_light
    %
    %   Class:
    %   seg_worm.worm
    
    properties (Dependent)
       parse_error 
    end
    methods
        function value = get.parse_error(obj)
           value =  obj.error_handler.error_found;
        end
    end
    
    properties
        error_handler  %seg_worm.parsing.parse_error
        orientation_handler
        frame_number   %(frame)
        original_image 
        
        contour  %seg_worm.worm.contour
        skeleton %seg_worm.worm.skeleton
        
        head %seg_worm.worm.head
        tail %seg_worm.worm.tail
        left_side  %seg_worm.worm.left
        right_side %seg_worm.worm.right

        manipulations = struct('code',{},'desc',{},'details',{})
    end
    
    methods
        function obj = worm(img, frame_number, isNormalized, verbose, varargin)
            
           %TODO: Create error object here, pass into intialization object 
           obj.error_handler       = seg_worm.parsing.parse_error(img,frame_number,verbose);
           
           obj.original_image      = img;
           obj.frame_number        = frame_number;
           
           %Initialization
           %---------------------------------------------------------------
           %seg_worm.worm.initialize
           obj.intialize(img, isNormalized, varargin{:}) 

           if obj.error_handler.error_found, return; end
           
           %TODO: orientation_handler sounds weird, why not just hold onto
           %an orientation object?????
           obj.orientation_handler = seg_worm.worm.orientation(obj);
        end
        function addManipulation(obj,code,description,details)
            %
            %   This can be used to document code manipulations ...
            %
           obj.manipulations = [obj.manipulations ...
               struct('code',code,'desc',description,'details',details)]; 
        end
    end
    
end

%{

        vWorm = worm2struct(frame_number, contour, [], [], [], lfCAngles, [], ...
            [], cCCLengths, [], [], [], [], [], [], [], [], [], [], [], ...
            [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], ...
            [], [], [], [], [], [], [], [], [], 0, [], [], 0, [], []);

worm = worm2struct(frame, ...
    cPixels, cTouchI, cInI, cOutI, cAngles, cHeadI, cTailI, cCCLengths, ...
    sPixels, sTouchI, sInI, sOutI, sInOutI, sAngles, ...
    sLength, sCCLengths, sWidths, ...
    hlcBounds, hrcBounds, hsBounds, hPixels, hArea, hCDF, hStdev, ...
    tlcBounds, trcBounds, tsBounds, tPixels, tArea, tCDF, tStdev, ...
    lcBounds, lsBounds, lPixels, lArea, lCDF, lStdev, ...
    rcBounds, rsBounds, rPixels, rArea, rCDF, rStdev, ...
    isHeadTailFlipped, hConfidence, tConfidence,...
    isVulvaClockwiseFromHead, vConfidence, nvConfidence)




%}


%
%       * Head *
%
%       hlcBounds - the worm head's, left-side (counter clockwise from the head),
%                   contour bounds (the start and end indices of the segment)
%       hrcBounds - the worm head's, right-side (clockwise from the head),
%                   contour bounds (the start and end indices of the segment)
%       hsBounds  - the worm head's, skeleton bounds (the start and end
%                   indices of the segment)
%                   Note: due to the clockwise ordering of the worm contour
%                   and the head-to-tail ordering of the worm skeleton,
%                   the bounds enclose the head as
%                   [hsBounds(1), hrcBounds(1:2), hsBounds(2), hlcBounds(1:2)]
%       hPixels   - the worm head's circularly continuous contour pixels
%       hArea     - the worm head's pixel area
%       hCDF      - the worm head's pixel-intensity, cumulative distribution
%                   function at 2.5%, 25%, 50%, 75%, and 97.5%
%       hStdev    - the worm head's pixel-intensity standard deviation
%
%       * Tail *
%
%       tlcBounds - the worm tail's, left-side (counter clockwise from the head),
%                   contour bounds (the start and end indices of the segment)
%       trcBounds - the worm tail's, right-side (clockwise from the head),
%                   contour bounds (the start and end indices of the segment)
%       tsBounds  - the worm tail's, skeleton bounds (the start and end
%                   indices of the segment)
%                   Note: due to the clockwise ordering of the worm contour
%                   and the head-to-tail ordering of the worm skeleton,
%                   the bounds enclose the tail as
%                   [tsBounds(1), trcBounds(1:2), tsBounds(2), tlcBounds(1:2)]
%       tPixels   - the worm tail's circularly continuous contour pixels
%       tArea     - the worm tail's pixel area
%       tCDF      - the worm tail's pixel-intensity, cumulative distribution
%                   function at 2.5%, 25%, 50%, 75%, and 97.5%
%       tStdev    - the worm tail's pixel-intensity standard deviation
%
%       * Left Side (Counter Clockwise from the Head) *
%
%       lcBounds - the worm's, left-side (counter clockwise from the head),
%                  contour bounds (the start and end indices of the segment)
%       lsBounds - the worm's, left-side (counter clockwise from the head),
%                  skeleton bounds (the start and end indices of the segment)
%                  Note: due to the clockwise ordering of the worm contour
%                  and the head-to-tail ordering of the worm skeleton,
%                  the bounds enclose the left side as
%                  [lcBounds(1:2), lsBounds(2:1)]
%       lPixels  - the worm's left-side (counter clockwise from the head)
%                 circularly continuous contour pixels
%       lArea    - the worm's left-side pixel area
%       lCDF     - the worm's left-side (counter clockwise from the head)
%                  pixel-intensity, cumulative distribution function at
%                  2.5%, 25%, 50%, 75%, and 97.5%
%       lStdev   - the worm's left-side (counter clockwise from the head)
%                  pixel-intensity standard deviation
%
%       * Right Side (Clockwise from the Head) *
%
%       rcBounds - the worm's, right-side (clockwise from the head),
%                  contour bounds (the start and end indices of the segment)
%       rsBounds - the worm's, right-side (clockwise from the head),
%                  skeleton bounds (the start and end indices of the segment)
%                  Note: due to the clockwise ordering of the worm contour
%                  and the head-to-tail ordering of the worm skeleton,
%                  the bounds enclose the left side as
%                  [rcBounds(1:2), rsBounds(2:1)]
%       rPixels  - the worm's right-side (clockwise from the head)
%                  circularly continuous contour pixels
%       rArea    - the worm's right-side pixel area
%       rCDF     - the worm's right-side (clockwise from the head)
%                  pixel-intensity, cumulative distribution function at
%                  2.5%, 25%, 50%, 75%, and 97.5%
%       rStdev   - the worm's right-side (clockwise from the head)
%                  pixel-intensity standard deviation
%
%       * Orientation *
%
%       isHeadTailFlipped        - are the head and tail flipped?
%                                  Note 1: the head and tail may be
%                                  incorrectly assigned. This flag, allows
%                                  the assignment to be easily flipped.
%                                  Note 2: this flag also flips the
%                                  skeleton orientation.
%       hConfidence              - how much confidence do we have in our
%                                  head choice as the worm's head?
%       tConfidence              - how much confidence do we have in our
%                                  tail choice as the worm's head?
%       isVulvaClockwiseFromHead - is the vulva on the side clockwise from
%                                  the head?
%       vConfidence              - how much confidence do we have in our
%                                  vulval-side choice as the worm's
%                                  vulval side?
%       nvConfidence             - how much confidence do we have in our
%                                  non-vulval-side choice as the worm's
%                                  vulval side?
