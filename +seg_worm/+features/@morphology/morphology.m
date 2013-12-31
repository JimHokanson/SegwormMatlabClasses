classdef morphology < handle
    %
    %   Class:
    %   seg_worm.features.morphology
    %
    %
    %
    %   Nature Methods Description
    %   =======================================================================
    %
    %   Morphology Features 
    %
    %   1. Length. Worm length is computed from the segmented skeleton by
    %   converting the chain-code pixel length to microns.
    %   
    %   2. Widths. Worm width is computed from the segmented skeleton. The
    %   head, midbody, and tail widths are measured as the mean of the widths
    %   associated with the skeleton points covering their respective sections.
    %   These widths are converted to microns.
    %   
    %   3. Area. The worm area is computed from the number of pixels within the
    %   segmented contour. The sum of the pixels is converted to microns2.
    %
    %   4. Area/Length.
    %
    %   5. Midbody Width/Length.
    
    properties (Hidden)
       nw 
    end
    
    properties
       length
       width
       %    .head
       %    .midbody
       %    .tail
       %
       area
       areaPerLength
       widthPerLength
    end
    
    methods
        function obj = morphology(nw,debug_options)
            %
            %   seg_worm.features.morphology
            %
            %   nw : seg_worm.normalized_worm
            
            SI = seg_worm.skeleton_indices;
            
            obj.nw = nw; %For plotting ...
            
            obj.length         = nw.lengths;
            obj.width.head     = mean(nw.widths(SI.HEAD_INDICES,:),1);
            obj.width.midbody  = mean(nw.widths(SI.MID_INDICES,:),1);
            obj.width.tail     = mean(nw.widths(SI.TAIL_INDICES,:),1);
            obj.area           = nw.head_areas + nw.tail_areas + nw.vulva_areas + nw.non_vulva_areas;
            obj.areaPerLength  = obj.area./obj.length;
            obj.widthPerLength = obj.width.midbody./obj.length;
           
            
        end
    end
    
end

