classdef tail < seg_worm.worm.head_tail
    %
    %   Class:
    %   seg_worm.worm.tail
    
    properties
       type = 'tail' 
    end
    
    methods
        function obj = tail(parent)
           obj@seg_worm.worm.head_tail(parent); 
        end
    end
    
end

%{

[tail,tlcBounds,trcBounds,tsBounds] = ...
    ...
    worm2poly(size(sPixels, 1), ...
    chainCodeLength2Index(sCCLengths(end) - htSSegLength, sCCLengths), ...
    sPixels, headI, tailI, cPixels, false, sCCLengths, cCCLengths);

% Determine the tail's MER (minimum enclosing rectangle).
tMinY = min(tail(:,1));
tMaxY = max(tail(:,1));
tMinX = min(tail(:,2));
tMaxX = max(tail(:,2));

% Measure the tail statistics.
merTImg = oImg(tMinY:tMaxY, tMinX:tMaxX);
merTail = [tail(:,1) - tMinY + 1, tail(:,2) - tMinX + 1];
[merTMask, merTailI] = inPolyMask(merTail, [], size(merTImg));
tColors = single(merTImg(merTMask));
tArea   = length(tColors);
tCDF    = prctile(tColors,[2.5 25 50 75 97.5]);
tStdev  = std(tColors);





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
%
%
%       tArea     - the worm tail's pixel area
%       tCDF      - the worm tail's pixel-intensity, cumulative distribution
%                   function at 2.5%, 25%, 50%, 75%, and 97.5%
%       tStdev    - the worm tail's pixel-intensity standard deviation


%}
