function [worm, errNum, errMsg] = ...
    intialize(obj,img, frame_number, isNormalized, verbose, varargin)
%SEGWORM Segment the worm in an image and organize the information in a
%   structure.
%
%   WORM = SEGWORM(IMG, FRAME, ISNORMALIZED, VERBOSE)
%
%   WORM = SEGWORM(IMG, FRAME, ISNORMALIZED, VERBOSE, NUMERODE, NUMDILATE)
%
%   WORM = SEGWORM(IMG, FRAME, ISNORMALIZED, VERBOSE, NUMERODE, NUMDILATE,
%                  SAMPLES, ISINTERP)
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


% The worm is roughly divided into 24 segments of musculature (i.e., hinges
% that represent degrees of freedom) on each side. Therefore, 48 segments
% around a 2-D contour.
% Note: "In C. elegans the 95 rhomboid-shaped body wall muscle cells are
% arranged as staggered pairs in four longitudinal bundles located in four
% quadrants. Three of these bundles (DL, DR, VR) contain 24 cells each,
% whereas VL bundle contains 23 cells." - www.wormatlas.org
S_WORM_SEGS = 24;

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

contour_obj = seg_worm.worm.contour(img,frame_number,verbose,num_erode,num_dilate,isNormalized);
if contour_obj.parse_error
   keyboard 
end

obj.contour = contour_obj;

skeleton_obj = seg_worm.worm.skeleton(contour_obj);
obj.skeleton = skeleton_obj;


%Temporary, move to using the variables eventually
sCCLengths = skeleton_obj.cc_lengths;
cCCLengths = contour_obj.cc_lengths;
cPixels    = contour_obj.pixels;
sPixels    = skeleton_obj.pixels;
headI      = contour_obj.head_I;
tailI      = contour_obj.tail_I;

obj.head = seg_worm.worm.head(obj);
obj.tail = seg_worm.worm.tail(obj);

[obj.left_side,obj.right_side] = seg_worm.worm.left_right.createSides(obj);

%JAH NOTE: I'm at this point although
%I need to change the error handler slightly and 
%I might also want to implement the head and tail confidence measure
%- see head_tail class

keyboard


% Compute the contour's local low-frequency curvature minima.
%[lfCMinP, lfCMinI] = minPeaksCircDist(lfCAngles, lfAngleEdgeLength, cCCLengths);

%{
cplot = seg_worm.video.pixel_list_image;
cplot.addList('main_contour',[1 1 1],contour_obj.pixels);
cplot.addList('head',[1 0 0],head);
cplot.addList('tail',[0 1 0],tail);
%}

%JAH NOTE: I just threw this function together
%to try and encapsulate things. It is likely that it is missing
%a lot of properties
helper__checkWormCoiled(contour_obj);





% Measure the skeleton angles (curvature).
lfAngleEdgeLength = sCCLengths(end) * 2 / S_WORM_SEGS;
sAngles           = curvature(sPixels, lfAngleEdgeLength, sCCLengths);






% Are the head and tail too small (or the body too large)?
% Note: earlier, the head and tail were each chosen to be 4/24 = 1/6
% the body length of the worm. The head and tail are roughly shaped
% like rounded triangles with a convex taper. And, the width at their
% ends is nearly the width at the center of the worm. Imagine they were
% 2 triangles that, when combined, formed a rectangle similar to the
% midsection of the worm. The area of this rectangle would be greater
% than a 1/6 length portion from the midsection of the worm (the
% maximum area per length in a worm is located at its midsection). The
% combined area of the right and left sides is 4/6 of the worm.
% Therefore, the combined area of the head and tail must be greater
% than (1/6) / (4/6) = 1/4 the combined area of the left and right
% sides.
if 4 * (hArea + tArea) < lArea + rArea
    errNum = 111;
    errMsg = ['The worm head and tail are less than 1/4 the size ' ...
        'of its remaining body. Therefore, the worm is ' ...
        'significantly obscured and cannot be segmented.'];
    
    % Defer organizing the available worm information.
    if verbose
        warning('segWorm:SmallHeadTail', ['Frame %d: ' errMsg], frame_number);
        vWorm = 0;
    else
        return;
    end
end

% How much confidence do we have in our vulva orientation?
% Note: generally, the vulval side contains less white pixels (a lower
% 50% and 75% CDF for color) and more gray pixels (a lower variance and
% 25% to 75% interquartile range) than the opposing side. We give each
% probability equal weight, then compare. Also, in the absence of
% information, we assume the vulva is on the left side (and use a trick
% to avoid reciprocals in our equations).
isVulvaClockwiseFromHead = 0; % default orientation
vConfidenceScale = 1048576; % 2^20
vConfidence = (rCDF(3) * rCDF(4) * rStdev * (rCDF(4) - rCDF(2))) ...
    / vConfidenceScale;
nvConfidence = (lCDF(3) * lCDF(4) * lStdev * (lCDF(4) - lCDF(2))) ...
    / vConfidenceScale;



keyboard

vWorm = [];
% Organize the available worm information.
if isempty(vWorm)
    worm = worm2struct(frame_number, cPixels, [], [], [], lfCAngles, ...
        headI, tailI, cCCLengths, sPixels, [], [], [], [], ...
        sAngles, sLength, sCCLengths, cWidths, ...
        hlcBounds, hrcBounds, hsBounds, head, hArea, hCDF, hStdev, ...
        tlcBounds, trcBounds, tsBounds, tail, tArea, tCDF, tStdev, ...
        lcBounds, sBounds, lSide, lArea, lCDF, lStdev, ...
        rcBounds, sBounds, rSide, rArea, rCDF, rStdev, ...
        isHeadTailFlipped, hConfidence, tConfidence, ...
        isVulvaClockwiseFromHead, vConfidence, nvConfidence);
else
    vWorm = worm2struct(frame_number, cPixels, [], [], [], lfCAngles, ...
        headI, tailI, cCCLengths, sPixels, [], [], [], [], ...
        sAngles, sLength, sCCLengths, cWidths, ...
        hlcBounds, hrcBounds, hsBounds, head, hArea, hCDF, hStdev, ...
        tlcBounds, trcBounds, tsBounds, tail, tArea, tCDF, tStdev, ...
        lcBounds, sBounds, lSide, lArea, lCDF, lStdev, ...
        rcBounds, sBounds, rSide, rArea, rCDF, rStdev, ...
        isHeadTailFlipped, hConfidence, tConfidence, ...
        isVulvaClockwiseFromHead, vConfidence, nvConfidence);
end



%=================================================================================





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

function helper__checkWormCoiled(obj)


keyboard

% Is the worm coiled?
% If there are no large concavities, the worm is not coiled.
lfCBendI = lfCMinI(lfCMinP < -30);
if ~isempty(lfCBendI)
    
    % Find concavities near the head. If there are any concavities
    % near the tail, the head may be portruding from a coil; in
    % which case, the width at the end of the head may be
    % inaccurate.
    if hlcBounds(1) < hrcBounds(2)
        hBendI = lfCBendI(lfCBendI > hlcBounds(1) & lfCBendI < hrcBounds(2));
    else
        hBendI = lfCBendI(lfCBendI > hlcBounds(1) | lfCBendI < hrcBounds(2));
    end
    
    cWidths = skeleton_obj.c_widths;
    
    % Does the worm more than double its width from the head?
    % Note: if the worm coils, its width will grow to more than
    % double that at the end of the head.
    maxWidth = max(cWidths);
    if isempty(hBendI)
        if maxWidth / cWidths(hsBounds(2)) > 2
            errNum = 107;
            errMsg = ['The worm more than doubles its width ' ...
                'from end of its head. Therefore, the worm is ' ...
                'coiled, laid an egg, and/or is significantly ' ...
                'obscured and cannot be segmented.'];
            
            % Organize the available worm information.
            if verbose
                warning('segWorm:DoubleHeadWidth', ...
                    ['Frame %d: ' errMsg], frame_number);
                vWorm = worm2struct(frame_number, cPixels, [], [], [], ...
                    lfCAngles, headI, tailI, cCCLengths, [], [], ...
                    [], [], [], [], [], [], [], [], [], [], [], [], ...
                    [], [], [], [], [], [], [], [], [], [], [], [], ...
                    [], [], [], [], [], [], [], [], [], 0, [], [], ...
                    0, [], []);
            else
                return;
            end
        end
    end
    
    % Find concavities near the tail. If there are any concavities near
    % the tail, the tail may be portruding from a coil; in which case,
    % the width at the end of the tail may be inaccurate.
    if trcBounds(1) < tlcBounds(2)
        tBendI = lfCBendI(lfCBendI > trcBounds(1) & lfCBendI < tlcBounds(2));
    else
        tBendI = lfCBendI(lfCBendI > trcBounds(1) | lfCBendI < tlcBounds(2));
    end
    
    % Does the worm more than double its width from the tail?
    % If the worm coils, its width will grow to more than double
    % that at the end of the tail.
    if isempty(tBendI)
        if maxWidth / cWidths(tsBounds(1)) > 2
            errNum = 108;
            errMsg = ['The worm more than doubles its width ' ...
                'from end of its tail. Therefore, the worm is ' ...
                'coiled, laid an egg, and/or is significantly ' ...
                'obscured and cannot be segmented.'];
            
            % Organize the available worm information.
            if verbose
                warning('segWorm:DoubleTailWidth', ...
                    ['Frame %d: ' errMsg], frame_number);
                vWorm = worm2struct(frame_number, cPixels, [], [], [], ...
                    lfCAngles, headI, tailI, cCCLengths, [], [], ...
                    [], [], [], [], [], [], [], [], [], [], [], [], ...
                    [], [], [], [], [], [], [], [], [], [], [], [], ...
                    [], [], [], [], [], [], [], [], [], 0, [], [], ...
                    0, [], []);
            else
                return;
            end
        end
    end
    
    lfCAngles = contour_obj.lf_angles;
    
    % Use the most accurate estimate of head/tail width to
    % determine whether the width of the body is more than double
    % that at the end of the head/tail; in which case; the worm is
    % coiled.
    if ~(isempty(hBendI) && isempty(tBendI))
        
        % Find the distances of bends near the head.
        hBendDist = abs(headI - hBendI);
        hBendDist = min(hBendDist, abs(hBendDist - length(lfCAngles)));
        
        % Find the distances of bends near the tail.
        tBendDist = abs(tailI - tBendI);
        tBendDist = min(tBendDist, abs(tBendDist - length(lfCAngles)));
        
        % The bend near the head is furthest and, therefore, the
        % width at the end of the head is our most accurate
        % estimate of the worm's width.
        if min(hBendDist) >= min(tBendDist)
            if maxWidth / cWidths(hsBounds(2)) > 2
                errNum = 107;
                errMsg = ['The worm more than doubles its width ' ...
                    'from end of its head. Therefore, the worm is ' ...
                    'coiled, laid an egg, and/or is significantly ' ...
                    'obscured and cannot be segmented.'];
                
                % Organize the available worm information.
                if verbose
                    warning('segWorm:DoubleHeadWidth', ...
                        ['Frame %d: ' errMsg], frame_number);
                    vWorm = worm2struct(frame_number, cPixels, [], [], [], ...
                        lfCAngles, headI, tailI, cCCLengths, [], ...
                        [], [], [], [], [], [], [], [], [], [], [], ...
                        [], [], [], [], [], [], [], [], [], [], [], ...
                        [], [], [], [], [], [], [], [], [], [], [], ...
                        [], 0, [], [], 0, [], []);
                else
                    return;
                end
            end
            
            % The bend near the tail is furthest and, therefore, the
            % width at the end of the tail is our most accurate
            % estimate of the worm's width.
        else
            if maxWidth / cWidths(tsBounds(1)) > 2
                errNum = 108;
                errMsg = ['The worm more than doubles its width ' ...
                    'from end of its tail. Therefore, the worm is ' ...
                    'coiled, laid an egg, and/or is significantly ' ...
                    'obscured and cannot be segmented.'];
                
                % Organize the available worm information.
                if verbose
                    warning('segWorm:DoubleTailWidth', ...
                        ['Frame %d: ' errMsg], frame_number);
                    vWorm = worm2struct(frame_number, cPixels, [], [], [], ...
                        lfCAngles, headI, tailI, cCCLengths, [], ...
                        [], [], [], [], [], [], [], [], [], [], [], ...
                        [], [], [], [], [], [], [], [], [], [], [], ...
                        [], [], [], [], [], [], [], [], [], [], [], ...
                        [], 0, [], [], 0, [], []);
                else
                    return;
                end
            end
        end
    end
end


end