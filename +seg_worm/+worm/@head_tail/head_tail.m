classdef head_tail < sl.obj.handle_light
    %
    %   Class:
    %   seg_worm.worm.head_tail
    
    properties (Constant)
        CDF_EVAL_LEVELS = [2.5 25 50 75 97.5]
    end
    
    properties (Abstract)
        type
    end
    
    
    properties (Hidden)
        parent
        contour  %seg_worm.worm.contour
        skeleton %seg_worm.worm.skeleton
        is_head %We'll use this, for now we'll assume
        %that the opposite of head is tail. I liked the
        %type display better but wanted the logical for processing
    end
    methods
        function value = get.is_head(obj)
            value = strcmp(obj.type,'head');
        end
    end
    
    properties
        contour_pixels %Specific to head or tail
        left_contour_bounds
        right_contour_bounds
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
    
    methods
        function obj = head_tail(parent)
            
            obj.parent   = parent;
            obj.contour  = parent.contour;
            obj.skeleton = parent.skeleton;
            obj.initializeBounds();
            obj.initStats();
        end
        function initializeBounds(obj)
            % Compute the worm's head and tail (at this point, we cannot
            % distinguish between the two). The worm's head and tail occupy,
            % approximately, 4 muscle segments each, on the skeleton and either
            % side of the contour.
            %
            % Note: "The first two muscle cells in the two ventral and two dorsal
            % rows [of the head] are smaller than their lateral counterparts,
            % giving a stagger to the packing of the two rows of cells in a
            % quadrant. The first four muscles in each quadrant are innervated
            % exclusively by motoneurons in the nerve ring. The second block of
            % four muscles is dually innervated, receiving synaptic input from
            % motoneurons in the nerve ring and the anterior ventral cord. The rest
            % of the muscles in the body are exclusively innervated by NMJs in the
            % dorsal and ventral cords." - The Structure of the Nervous System of
            % the Nematode C. elegans, on www.wormatlas.org
            
            s = obj.skeleton;
            c = obj.contour;
            
            ht_s_seg_length = s.cc_lengths(end)*4/s.N_SEGS;
            
            if obj.is_head
                startSI = 1;
                endSI   = seg_worm.cv.chainCodeLength2Index(ht_s_seg_length, s.cc_lengths);
            else
                startSI = size(s.pixels,1);
                endSI   = seg_worm.cv.chainCodeLength2Index(s.cc_lengths(end) - ht_s_seg_length, s.cc_lengths);
            end
            
            [obj.contour_pixels,obj.left_contour_bounds,...
                obj.right_contour_bounds,obj.skeleton_bounds] = ...
                worm2poly(startSI,endSI,s.pixels,c.head_I,c.tail_I,c.pixels,s.cc_lengths,c.cc_lengths);
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

% How much confidence do we have in our head-to-tail orientation?
% Note: generally, the head is less angled, and contains more white
% pixels (a higher 50% and 75% CDF for color) and less gray pixels (a higher
% variance and 25% to 75% interquartile range) than the tail. We give
% each probability equal weight, then compare.
isHeadTailFlipped = 0; % default orientation
hConfidenceScale = 134217728; % 2^26
hConfidence = ((180 - lfCAngles(headI)) * hCDF(3) * hCDF(4) * ...
    hStdev * (hCDF(4) - hCDF(2))) / hConfidenceScale;
tConfidence = ((180 - lfCAngles(tailI)) * tCDF(3) * tCDF(4) * ...
    tStdev * (tCDF(4) - tCDF(2))) / hConfidenceScale;



%}