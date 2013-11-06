classdef skeleton < sl.obj.handle_light
    %
    %   Class:
    %   seg_worm.worm.skeleton
    
    
    properties (Constant)
        % The worm is roughly divided into 24 segments of musculature (i.e., hinges
        % that represent degrees of freedom) on each side. Therefore, 48 segments
        % around a 2-D contour.
        % Note: "In C. elegans the 95 rhomboid-shaped body wall muscle cells are
        % arranged as staggered pairs in four longitudinal bundles located in four
        % quadrants. Three of these bundles (DL, DR, VR) contain 24 cells each,
        % whereas VL bundle contains 23 cells." - www.wormatlas.org 
       N_SEGS = 24 
    end
    
    properties
        pixels     %(sPixels) - the worm's continuous skeleton oriented from head to tail
        cc_lengths %%(sCCLengths) - chain code lengths of the skeleton pixels
        
        %.computeChainCodeLengths()
        c_widths   %(sWidths) - the contour's (worm's) widths, from head to tail, at
        % each skeleton point
        
        %.computeAngles()
        angles %(sAngles) - the skeleton's angles (curvature) per point
        % Note 1: NaNs indicate the end pieces where there is
        % insufficient information to compute the angle
        % Note 2: positive skeleton angles bulge towards the side
        % clockwise from the worm's head (unless the worm is flipped)
    end
    properties (Dependent)
       total_length %(sLength) - the total chain-code length
    end
    methods
        function value = get.total_length(obj)
           value = obj.cc_lengths(end); 
        end
        function set.pixels(obj,value)
           obj.pixels = value;
           obj.computeChainCodeLengths();
           obj.computeAngles();
        end
    end
     
    properties
        
        %When are these populated????
        
        touch_I %(sTouchI) - the paired pairs of indices marking, clockwise, the
        % start and end of the touching skeleton points
        % Note: if the worm isn't coiled, this value is empty.
        inner_I %(sInI) - the paired indices marking, clockwise, the start and
        % end of the inner skeleton points
        % Note: if the worm isn't coiled, this value is empty.
        outer_I %(sOutI) - the paired indices marking, clockwise, the start and
        % end of the outer skeleton points
        % Note: if the worm isn't coiled, this value is empty.
        inner_outer_I %(sInOutI) - the pairs of indices marking, from head to tail,
        % the start and end of the dual inner/outer skeleton points
        % Note: if the worm isn't coiled, this value is empty.
        
        


    end
    
    methods
        function obj = skeleton(contour)
           %seg_worm.worm.skeleton.initialize
           obj.initialize(contour) 
        end
        function computeChainCodeLengths(obj)
           obj.cc_lengths = seg_worm.cv.computeChainCodeLengths(obj.pixels); 
        end
        function computeAngles(obj)
           % Measure the skeleton angles (curvature).
            lfAngleEdgeLength = obj.cc_lengths(end) * 2 / obj.N_SEGS;
            obj.angles  = seg_worm.cv.curvature(obj.pixels, lfAngleEdgeLength, obj.cc_lengths);
        end
    end
    
end

