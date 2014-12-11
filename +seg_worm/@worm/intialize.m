function intialize(obj,img, is_normalized, varargin)
%x Segment the worm in an image
%
%
%   These have been lost in translation. Build an options object which
%   exposes these values ...
%   *NUMERODE, *NUMDILATE, *SAMPLES, *ISINTERP)
%
%   Inputs:
%   -------
%   img : 
%       The image to segment
%   is_normalized : logical
%       Is the image already normalized (i.e., all pixel values are 
%       between 0 to 1, inclusive)?
%
%   
%   These have been lost in translation and eventually should be
%   reintroduced.
%   
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
%
%   See Also:
%   seg_worm.segWormFrames
%   seg_worm.worm.contour

%Old Code:
%---------
%/SegWorm/Pipeline/segmentationMain.m


import seg_worm.cv.*

%Placeholders. Eventually we should pass these in ...
n_erode  = [];
n_dilate = [];


%Type: seg_worm.parse_error
error_handler = obj.error_handler;

%1) Create the contour
%--------------------------------------------------------------------------
%Lines - 93 - 537 of:
%https://github.com/openworm/SegWorm/blob/master/Worms/Segmentation/segWorm.m
contour_obj = seg_worm.worm.contour(obj,img,n_erode,n_dilate,is_normalized);
if error_handler.error_found, return; end

obj.contour = contour_obj;

%2) Create skeleton from contour
%--------------------------------------------------------------------------
skeleton_obj = seg_worm.worm.skeleton(contour_obj);
obj.skeleton = skeleton_obj;

%3) Parse 
%Lines 552 - 796:
%https://github.com/openworm/SegWorm/blob/master/Worms/Segmentation/segWorm.m
obj.head = seg_worm.worm.head(obj);
obj.tail = seg_worm.worm.tail(obj);

%Lines 
%https://github.com/openworm/SegWorm/blob/master/Worms/Segmentation/segWorm.m#L798
error_handler.checkHeadTailArea(obj);
if error_handler.error_found, return; end

%Lines 884 - 921
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

%This looks like it may be out of order:
%Looks like it comes in at line 579
obj.checkIfWormIsCoiled();
if error_handler.error_found, return; end

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
