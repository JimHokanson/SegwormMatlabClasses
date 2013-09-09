function [worm, errNum, errMsg] = ...
    intialize(obj,img, frame_number, isNormalized, verbose, varargin)
%SEGWORM Segment the worm in an image and organize the information in a
%   structure.
%
%   WORM = SEGWORM(IMG, FRAME, ISNORMALIZED, VERBOSE, *NUMERODE, *NUMDILATE, *SAMPLES, *ISINTERP)
%
%   Inputs:
%       img          - the image to segment
%       frame        - the frame number (if the image comes from video)
%       isNormalized - is the image already normalized (i.e., all pixel
%                      values are between 0 to 1, inclusive)?
%       verbose      - verbose mode shows the results in a figure
%       numErode     - the number of time to erode the binary worm image;
%                      if empty or 0, the binary image is left unchanged
%       numDilate    - the number of time to dilate the binary worm image;
%                      if empty or 0, the binary image is left unchanged
%       samples      - the number of samples to use in verbose mode;
%                      if empty, all the worm is used.
%       isInterp     - when downsampling, should we interpolate the missing
%                      data or copy it from the original worm;
%                      if empty, we interpolate the missing data.
%
%   Output:
%       worm - the worm information organized in a structure
%              This structure contains 8 sub-structures,
%              6 sub-sub-structures, and 4 sub-sub-sub-structures:
%
%              * Video *
%              video = {frame}
%
%              * Contour *
%              contour = {pixels, touchI, inI, outI, angles, headI, tailI}
%
%              * Skeleton *
%              skeleton = {pixels, touchI, inI, outI, inOutI, angles,
%                          length, chainCodeLengths, widths}
%
%              Note: positive skeleton angles bulge towards the side
%              clockwise from the worm's head (unless the worm is flipped).
%
%              * Head *
%              head = {bounds, pixels, area,
%                      cdf (at [2.5% 25% 50% 75% 97.5%]), stdev}
%              head.bounds{contour.left (indices for [start end]),
%                          contour.right (indices for [start end]),
%                          skeleton indices for [start end]}
%
%              * Tail *
%              tail = {bounds, pixels, area,
%                      cdf (at [2.5% 25% 50% 75% 97.5%]), stdev}
%              tail.bounds{contour.left (indices for [start end]),
%                          contour.right (indices for [start end]),
%                          skeleton indices for [start end]}
%
%              * Left Side (Counter Clockwise from the Head) *
%              left = {bounds, pixels, area,
%                      cdf (at [2.5% 25% 50% 75% 97.5%]), stdev}
%              left.bounds{contour (indices for [start end]),
%                          skeleton (indices for [start end])}
%
%              * Right Side (Clockwise from the Head) *
%              right = {bounds, pixels, area,
%                       cdf (at [2.5% 25% 50% 75% 97.5%]), stdev}
%              right.bounds{contour (indices for [start end]),
%                           skeleton (indices for [start end])}
%
%              * Orientation *
%              orientation = {head, vulva}
%              orientation.head = {isFlipped,
%                                  confidence.head, confidence.tail}
%              orientation.vulva = {isClockwiseFromHead,
%                                  confidence.vulva, confidence.nonVulva}
%
%       errNum - the error number if segmentation failed
%                (see also WORMFRAMEANNOTATION)
%       errMsg - the error message if segmentation failed
%
%   See also WORM2STRUCT, NORMWORMS
%
%   See Also:
%   segWormFrames
%
%
%
% © Medical Research Council 2012
% You will not remove any copyright or other notices from the Software;
% you must reproduce all copyright notices and other proprietary
% notices on any copies of the Software.

import seg_worm.cv.*

%INPUT HANDLING
%--------------------------------------------------------------------------
% Are we eroding and/or dilating the worm?
num_erode  = [];
num_dilate = [];
if ~isempty(varargin)
    num_erode = varargin{1};
end
if length(varargin) > 1
    num_dilate = varargin{2};
end

if ~exist('isNormalized','var') || isempty(isNormalized)
    isNormalized = true;
end
%--------------------------------------------------------------------------
error_handler = obj.error_handler;

contour_obj = seg_worm.worm.contour(obj,img,num_erode,num_dilate,isNormalized);
if error_handler.error_found, return; end

obj.contour = contour_obj;

skeleton_obj = seg_worm.worm.skeleton(contour_obj);
obj.skeleton = skeleton_obj;

obj.head = seg_worm.worm.head(obj);
obj.tail = seg_worm.worm.tail(obj);

error_handler.checkHeadTailArea(obj);
if error_handler.error_found, return; end

[obj.left_side,obj.right_side] = seg_worm.worm.left_right.createSides(obj);

error_handler.headTailSmallOrBodyLarge(obj);
if error_handler.error_found, return; end


%Some temporary plotting I am messing around with ...
%{
cplot = seg_worm.video.pixel_list_image;
cplot.addList('skeleton',[1 1 1],obj.skeleton.pixels);
cplot.addList('main_contour',[0 0 1],obj.contour.pixels);
plot(cplot)

cplot = seg_worm.video.pixel_list_image;
cplot.addList('main_contour',[1 1 1],contour_obj.pixels);
cplot.addList('head',[1 0 0],obj.head.contour_pixels);
cplot.addList('tail',[0 1 0],obj.tail.contour_pixels);
plot(cplot)
%}

obj.checkIfWormIsCoiled();
if error_handler.error_found, return; end


% Organize the available worm information.
%{
worm = worm2struct(frame_number, cPixels, [], [], [], lfCAngles, ...
headI, tailI, cCCLengths, sPixels, [], [], [], [], ...
sAngles, sLength, sCCLengths, cWidths, ...
hlcBounds, hrcBounds, hsBounds, head, hArea, hCDF, hStdev, ...
tlcBounds, trcBounds, tsBounds, tail, tArea, tCDF, tStdev, ...
lcBounds, sBounds, lSide, lArea, lCDF, lStdev, ...
rcBounds, sBounds, rSide, rArea, rCDF, rStdev, ...
isHeadTailFlipped, hConfidence, tConfidence, ...
isVulvaClockwiseFromHead, vConfidence, nvConfidence);
%}



%=================================================================================


end

function helper__unknownCode(obj) %#ok<INUSD,DEFNU>
%
%
%   This looks like temporary code that was trying to resolve an inner
%   contour but that was never finished.

% Get the inner contour, if it exists.
if ~verbose && isempty(worm)
    
    warning('segWorm:CannotSegment', ...
        'Frame %d: The worm cannot be segmented', frame_number);
    return;
    
    % Create a small image of the worm's complement to locate its inner loop.
    minY = min(contour(:,1));
    maxY = max(contour(:,1));
    minX = min(contour(:,2));
    maxX = max(contour(:,2));
    cImg = ~img(minY:maxY, minX:maxX);
    
    % Order the connected components by size.
    cc = bwconncomp(cImg);
    ccSizes = zeros(length(cc.PixelIdxList), 1);
    for i = 1:length(cc.PixelIdxList)
        ccSizes(i) = length(cc.PixelIdxList{i});
    end
    [~, o] = sort(ccSizes, 1, 'descend');
    
    % Choose the largest connected component that doesn't touch any image edges.
    iWormPixels = [];
    wMinX = 0;
    wMaxX = 0;
    wMinY = 0;
    wMaxY = 0;
    for i = o'
        [y x] = ind2sub(size(cImg), cc.PixelIdxList{i});
        wMinX = min(x);
        wMaxX = max(x);
        wMinY = min(y);
        wMaxY = max(y);
        if (wMinY > 1 && wMaxY < size(cImg,1) && ...
                wMinX > 1 && wMaxX < size(cImg,2))
            iWormPixels = [y x];
            break;
        end
    end
    
    % Use the inner loop's maxima/minima to find a point on the worm's
    % inner contour. Try to avoid 1-pixel width sections as they may end up
    % tracing the outer contour.
    if ~isempty(iWormPixels)
        if wMinX > 2
            i = find(iWormPixels(:,2) == wMinX);
            x = wMinX + minX - 2;
            y = iWormPixels(i(1),1) + minY - 1;
        elseif wMinY > 2
            i = find(iWormPixels(:,2) == wMinY);
            x = iWormPixels(i(1),2) + minX - 1;
            y = wMinY + minY - 2;
        elseif wMaxX < size(cImg, 2) - 1
            i = find(iWormPixels(:,2) == wMaxX);
            x = wMaxX + minX;
            y = iWormPixels(i(1),1) + minY - 1;
        else % default to wMaxY < size(cImg, 1) - 1
            i = find(iWormPixels(:,2) == wMaxY);
            x = iWormPixels(i(1),2) + minX - 1;
            y = wMaxY + minY;
        end
        
        % Trace the contour counter clockwise.
        iContour = bwClockTrace(img, [x y], false);
        
        % Correct the worm segment size.
        % Note: a looped worm can take the shape of an upper-case omega and
        % omicron or, lower-case alpha and delta. An omega and omicron
        % hide, approximately, at least 2 segments of the worm's contour. A
        % delta hides, approximately, at least 3 segments. And, and alpha
        % hides, approximately, at least 4 segments. Therefore, we
        % approximate the new contour size as the inner and outer contour
        % sizes plus 3 additional, hidden segments.
        wormSegLength = round((size(iContour, 1) + size(contour, 1)) / ...
            (cWormSegs - 3));
        
        % Clean up the worm's inner contour.
        if verbose
            roughIContour = iContour;
        end
        [iContour] = cleanWorm(iContour, wormSegLength);
        
        % Are the inner and outer contour identical or switched?
        % Note: if we begin tracing either contour at a 1-pixel wide worm
        % section we may end up tracing the wrong contour when we fork.
        if size(iContour, 1) >= size(contour, 1)
            iMinX = min(iContour(:,2));
            iMaxX = max(iContour(:,2));
            iMinY = min(iContour(:,1));
            iMaxY = max(iContour(:,1));
            oMinX = min(contour(:,2));
            oMaxX = max(contour(:,2));
            oMinY = min(contour(:,1));
            oMaxY = max(contour(:,1));
            
            % The contours are switched.
            if iMinX < oMinX || iMaxX > oMaxX || ...
                    iMinY < oMinY || iMaxY > oMaxY
                
                % Switch the contours.
                tmp = flipud(contour);
                contour = flipud(iContour);
                iContour = tmp;
                
                % The contours are identical.
            elseif iMinX == oMinX && iMaxX == oMaxX && ...
                    iMinY == oMinY && iMaxY == oMaxY
                warning('segWorm:IdenticalContours', ...
                    ['Frame ' num2str(frame) ...
                    ': The inner and outer contour cannot be '...
                    'distinguished from each other']);
                return;
            end
        end
        
        % Compute the worm's contour and skeleton.
        worm = coiledSkeleton(headI, tailI, contour, iContour, wormSegLength);
        
        % Orient the contour and angles at the maximum curvature (the head or tail).
        % FIXME!!!
        if 0
            contour = [contour(headI:end,:); contour(1:(headI - 1),:)];
            hfCAngles = [hfCAngles(headI:end), hfCAngles(1:(headI - 1))];
            mhfCAngles = [mhfCAngles(headI:end), mhfCAngles(1:(headI - 1))];
            
            keyboard
            %NOTE: note, if headI and tailI change, I need
            %to update them in the skeleton (or contour?)
            
            if headI <= tailI
                tailI = tailI - headI + 1;
            else
                tailI = tailI + size(contour, 1) - headI + 1;
            end
            headI = 1;
        end
    end
end


end
