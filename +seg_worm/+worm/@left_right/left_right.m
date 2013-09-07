classdef left_right < sl.obj.handle_light
    %
    %   Class:
    %   seg_worm.worm.left_right
    %
    %   TODO: Inherit from a shared class
    %   with head_tail 
    %
    %   body_side?
    %
    %   Then move shared code to parent
    
    properties
       parent
       contour
       skeleton
    end
        
    properties
        contour_pixels
%         left_contour_bounds
%         right_contour_bounds
        skeleton_bounds %[2 x 1]
        %NOTE: This varies from head to tail
        %
        %   Should I change this so that 1
        %   is inner and 2 is outer????
        %
        %   This would effect the left_right 
        %   bound assignment
        
        %TODO: Hold onto color or mask???
        
        pixel_area
        pixel_cdf %Colors at each percentile
        %NOTE: Ideally we have a lot of dark values, 0
        %255 is the worst ...
        %This data is not currently being used ...
        
        %TODO: Make constant
        pixel_std_dev 
    end
    
    properties (Abstract)
       type
    end
    properties (Hidden)
       is_right 
    end
    methods
        function value = get.is_right(obj)
            value = strcmp(obj.type,'right');
        end
    end
    
    
    methods (Static)
        function [left,right] = createSides(parent)
            %
            %   [left,right] = seg_worm.worm.left_right.createSides(parent)
           left  = seg_worm.worm.left;
           right = seg_worm.worm.right;
           left.initRefs(parent);
           right.initRefs(parent);
           %left.initBounds(right);
        
           head = parent.head;
           tail = parent.tail;
           %tsBounds and hsBounds
           %.skeleton_bounds()
           
           %TODO: This should be changed to be more encapsulated
           %Also, the 1,2 is confusing, might make 1 inner and 2 outer
           c = parent.contour;
           s = parent.skeleton;
           inner_skeleton_bound_head = head.skeleton_bounds(2);
           inner_skeleton_bound_tail = tail.skeleton_bounds(1);
           sPixels  = parent.skeleton.pixels;
           cPixels  = c.pixels;
           headI    = c.head_I;
           tailI    = c.tail_I;
           cCCLengths = c.cc_lengths;
           sCCLengths = s.cc_lengths;
           
           IS_SPLIT = true;
           
           RIGHT_SIDE = 1; %From previous code, not sure why this is ...
           LEFT_SIDE  = 2;
           %When split is true, the skeleton creates a divide
           %such that the first output is a cell array of size 2
           
           %[polygon,lcBounds,rcBounds,sBounds] = worm2poly(startSI, endSI, ...
                %skeleton, headCI, tailCI, contour, isSplit, varargin)
            [sides,lcBounds,rcBounds,sBounds] = worm2poly(...
                inner_skeleton_bound_head, ...
                inner_skeleton_bound_tail, ...
                sPixels, headI, tailI, cPixels, ...
                IS_SPLIT, sCCLengths, cCCLengths);  %#ok<ASGLU>
            
            
            left.initBoundsAndContour(sides{LEFT_SIDE},sBounds);
            right.initBoundsAndContour(sides{RIGHT_SIDE},sBounds);
        
        end
    end
    methods
        function initBoundsAndContour(obj,c_pixels,s_bounds)
            obj.contour_pixels  = c_pixels;
            obj.skeleton_bounds = s_bounds;
        end
        function initRefs(obj,parent)
            obj.parent   = parent;
            obj.contour  = parent.contour;
            obj.skeleton = parent.skeleton;
        end
        function initStats(obj)
            
            p = obj.contour_pixels;
            minY = min(p(:,1));
            maxY = max(p(:,1));
            minX = min(p(:,2));
            maxX = max(p(:,2));
            
            original_image = obj.parent.original_image;
            
            subImg     = original_image(minY:maxY, minX:maxX);
            pixels     = obj.contour_pixels;
            adjustedContour = [pixels(:,1) - minY + 1, pixels(:,2) - minX + 1];
            
            [mask,~] = seg_worm.cv.inPolyMask(adjustedContour, [], size(subImg));
            
            colors            = single(subImg(mask));
            obj.pixel_area    = length(colors);
            obj.pixel_cdf     = prctile(colors,obj.CDF_EVAL_LEVELS);
            obj.pixel_std_dev = std(colors);
        end
    end
    
end

%{

% Determine the left-side's MER (minimum enclosing rectangle).
[sides,lcBounds,rcBounds,sBounds] = worm2poly(hsBounds(2), ...
    tsBounds(1), sPixels, headI, tailI, cPixels, true, ...
    sCCLengths, cCCLengths);
lSide = sides{2};
lMinY = min(lSide(:,1));
lMaxY = max(lSide(:,1));
lMinX = min(lSide(:,2));
lMaxX = max(lSide(:,2));

% Measure the left side (counter clockwise from the head) statistics.
%lCDF = [];
%lStdev = [];
merLImg       = oImg(lMinY:lMaxY, lMinX:lMaxX);
merLSide      = [lSide(:,1) - lMinY + 1, lSide(:,2) - lMinX + 1];
[merLMask, ~] = inPolyMask(merLSide, [], size(merLImg));
lColors       = single(merLImg(merLMask));
lArea         = length(lColors);
lCDF          = prctile(lColors,[2.5 25 50 75 97.5]);
lStdev        = std(lColors);

% Determine the right-side's MER (minimum enclosing rectangle).
rSide = sides{1};
rMinY = min(rSide(:,1));
rMaxY = max(rSide(:,1));
rMinX = min(rSide(:,2));
rMaxX = max(rSide(:,2));


% Measure the right side (clockwise from the head) statistics.
%rCDF = [];
%rStdev = [];
merRImg       = oImg(rMinY:rMaxY, rMinX:rMaxX);
merRSide      = [rSide(:,1) - rMinY + 1, rSide(:,2) - rMinX + 1];
[merRMask, ~] = inPolyMask(merRSide, [], size(merRImg));
rColors       = single(merRImg(merRMask));
rArea         = length(rColors);
rCDF          = prctile(rColors,[2.5 25 50 75 97.5]);
rStdev        = std(rColors);



%}