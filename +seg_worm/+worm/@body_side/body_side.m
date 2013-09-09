classdef body_side < sl.obj.handle_light
    %
    %   Class:
    %   seg_worm.worm.body_side
    
    properties (Constant)
        CDF_EVAL_LEVELS = [2.5 25 50 75 97.5]
    end
    
    properties (Hidden)
        parent   %seg_worm.worm
        contour  %seg_worm.worm.contour
        skeleton %seg_worm.worm.skeleton
    end
    
    properties (Abstract)
        type
    end
    
    properties
        contour_pixels
        skeleton_bounds %[2 x 1]
        %NOTE: This varies from head to tail
        %
        %   Should I change this so that 1
        %   is inner and 2 is outer????
        %
        %   This would effect the left_right
        %   bound assignment
        
        pixel_area
        pixel_cdf %Colors at each percentile
        %NOTE: Ideally we have a lot of dark values, 0
        %255 is the worst ...
        %This data is not currently being used ...
        
        %TODO: Make constant
        pixel_std_dev
        
    end
    
    methods
        function initRefs(obj,parent)
            obj.parent   = parent;
            obj.contour  = parent.contour;
            obj.skeleton = parent.skeleton;
        end
        function initStats(obj)
            
            p    = obj.contour_pixels;
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

