classdef head < seg_worm.worm.head_tail
    %
    %   Class:
    %   seg_worm.worm.head
    
    properties
       type = 'head'
    end
    

    
    methods
        function obj = head(parent)
           %
           %    head = seg_worm.worm.head(contour,skeleton)
           obj@seg_worm.worm.head_tail(parent); 
        end
    end
    
end

%{

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



%}



%{

%[head,hlcBounds,hrcBounds,hsBounds] = ...
%     ...
%     worm2poly(1, chainCodeLength2Index(htSSegLength, sCCLengths), ...
%     sPixels, headI, tailI, cPixels, false, sCCLengths, cCCLengths);

% Determine the head's MER (minimum enclosing rectangle).
hMinY = min(head(:,1));
hMaxY = max(head(:,1));
hMinX = min(head(:,2));
hMaxX = max(head(:,2));

% Measure the head statistics.
oImg    = img;
merHImg = oImg(hMinY:hMaxY, hMinX:hMaxX);
merHead = [head(:,1) - hMinY + 1, head(:,2) - hMinX + 1];

[merHMask, merHeadI] = inPolyMask(merHead, [], size(merHImg));

hColors = single(merHImg(merHMask));
hArea   = length(hColors);
hCDF    = prctile(hColors,[2.5 25 50 75 97.5]);
hStdev  = std(hColors);


%}





