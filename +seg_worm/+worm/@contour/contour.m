classdef contour < sl.obj.handle_light
    %
    %   Class:
    %   seg_worm.worm.contour
    
    %IMPROVEMENTS:
    %---------------------------------------------------------------------
    %Add a log which documents the various fixes which were necessary
    
    %MLINT:
    %------------------------------
    %#ok<*MCSUP> Ignore set warnings
    
    properties (Constant)
        N_SEGS = 48 %NOTE: This should be 2x
        LF_MAX_ANGLE_LIMIT = 90;
        HF_MAX_ANGLE_LIMIT = 60;
        MIN_ANGLE_LIMIT    = -60;
    end
    
    %Error properties =====================================================
    properties (Hidden)
        parent
    end
    
    properties (Dependent)
        error_handler %seg_worm.parse_error
        parse_error
    end
    methods
        function value = get.error_handler(obj)
            value = obj.parent.error_handler;
        end
        function value = get.parse_error(obj)
            value = obj.error_handler.error_found;
        end
    end
    
    %Input & Image properties =============================================
    properties
        d0 = '----  Input properties  ----'
        input_image
        frame_number
        
        d1 = '----  Image properties  ----'
        final_image %
        %set: .initializeFirstPassContour()
        
        pixels %[n x 2], I,J location pairs of each contour location
        %Each row entry corresponds to a pixel that is next to the previous
        %pixel.
        %Modified/Set in:
        %.initializeFirstPassContour()
        %.cleanContour()
        
        cc_lengths %(cCCLengths)[n x 1] - chain-code length, more accurate
        % than just using an index since diagonals are apart by sqrt(2)
        % instead of just 1. First value is the distance from the first to the last point.
        %
        %   Updated automatically when setting 'pixels'
        
        old_contour_pixels %{1 x n} Whenever pixels is updated, this is added
        %with the previous value of pixels
    end
    properties (Dependent)
        n_pixels
        avg_segment_length
    end
    methods
        function value = get.n_pixels(obj)
           value = size(obj.pixels,1); 
        end
        function value = get.avg_segment_length(obj)
            value = obj.cc_lengths(end)/obj.N_SEGS;
        end
        function set.pixels(obj,value)
            if ~isempty(obj.pixels)
                %Store old pixel values
                obj.old_contour_pixels = [obj.old_contour_pixels {obj.pixels}];
            end
            obj.pixels = value;
            obj.computeCCLengths();
        end
    end
    
    properties
        d2 = '----  Angle Information  ----'
        lf_angles  %(cAngles) the contour's angles (curvature) at each index
        hf_angles_raw
        hf_angles
    end
    
    properties
        d3 = '----  Angle Peak (ap) Info  ----'
        lf_ap_min
        lf_ap_min_I
        lf_ap_max
        lf_ap_max_I
        
        hf_ap_min
        hf_ap_min_I %smoothed
        hf_ap_max
        hf_ap_max_I
    end
    
    properties
        head_I  %(cHeadI)  the contour index for the worm's head
        tail_I  %(cTailI)  the contour index for the worm's tail
        
        %More advanced properties
        %--------------------------------------------------------------
        touch_points_I  %(cTouchI) the paired pairs of indices marking, clockwise, the
        % start and end of the touching contour points Note: if the worm isn't
        % coiled, this value is empty.
        inner_I %(cInI) the paired indices marking, clockwise, the start and
        % end of the inner contour points
        outer_I %(cOutI) the paired indices marking, clockwise, the start and
        % end of the outer contour points
    end
    
    %Initialization methods  ==============================================
    methods
        function obj = contour(parent,img,num_erode,num_dilate,is_normalized)
            %For debugging allow no inputs
            if ~nargin
                return
            end
            
            obj.parent = parent;
            obj.input_image  = img;
            
            obj.initializeFirstPassContour(num_erode,num_dilate,is_normalized)
            if obj.parse_error, return; end
            
            %This corresponds to line 250 at:
            %https://github.com/openworm/SegWorm/blob/master/Worms/Segmentation/segWorm.m#L250
            obj.cleanWorm();
            obj.error_handler.contourTooSmall(obj.pixels,obj.N_SEGS,img);
            if obj.parse_error, return; end
            
            obj.computeAngleInfo();
            
            obj.initTailHeadIndices();
        end
        function initializeFirstPassContour(obj,num_erode,num_dilate,is_normalized)
            %
            %
            %   1) Grayscale
            %   2) OTSU
            %   3) Dilate and erode
            %   4) Identification of biggest blob
            %   5) Outer contour tracing
            %
            %   Populates:
            %   -------------------------
            %   .rough_pixels
            %
            %   FULL PATH:
            %   seg_worm.worm.contour.initializeFirstPassContour
            
            pe  = obj.error_handler;
            img = obj.input_image;
            
            % Convert the image to grayscale.
            if (size(img,3) == 3)
                img = rgb2gray(img);
            end
            
            % Store the original then binarize the image.
            img  = seg_worm.cv.otsuImg(img, is_normalized);
            img  = ~img; %invert image, normally the pixels are white, the background black
            %which is different from
            
            % Erode and dilate the binary image.
            %-------------------------------------------
            if ~isempty(num_erode) && num_erode > 0
                img = imerode(img, strel('disk', num_erode));
            end
            if ~isempty(num_dilate) && num_dilate > 0
                img = imdilate(img, strel('disk', num_dilate));
            end
            
            obj.final_image = img;
            
            %Find the worm blob
            %----------------------------------------
            %The largest connected blob is the worm
            cc = bwconncomp(img); %img tlbx call
            worm_pixels = [];
            if ~isempty(cc.PixelIdxList)
                %PixelIdxList - cell array of arrays, each array
                %contains indices of connected pixels
                [~,I] = max(cellfun('length',cc.PixelIdxList));
                worm_pixels = cc.PixelIdxList{I};
            end
            
            pe.noWormFound(worm_pixels);
            if obj.parse_error, return; end
            
            % Find a point on the contour.
            %----------------------------------------
            [y, x] = ind2sub(size(img), min(worm_pixels));
            seed   = [x y];
            
            % Trace the contour clockwise.
            %----------------------------------------
            contour_pixels = seg_worm.cv.bwClockTrace(img,seed,true);
            
            pe.contourTouchesBoundary(contour_pixels,img);
            if obj.parse_error, return; end
            
            pe.contourTooSmall(contour_pixels,obj.N_SEGS,img);
            if obj.parse_error, return; end
            
            obj.pixels = contour_pixels;
        end
        function computeAngleInfo(obj)
            % Compute the contour's local high/low-frequency curvature.
            %
            % Note: worm body muscles are arranged and innervated as staggered pairs.
            % Therefore, 2 segments have one theoretical degree of freedom (i.e. one
            % approximation of a hinge). In the head, muscles are innervated
            % individually. Therefore, we sample the worm head's curvature at twice the
            % frequency of its body.
            %
            % Note 2: we ignore Nyquist sampling theorem (sampling at twice the
            % frequency) since the worm's cuticle constrains its mobility and practical
            % degrees of freedom. JAH NOTE: What spatial resolution would
            % not be appropriate?
            
            
            
            
            USE_MAX = true;
            USE_MIN = ~USE_MAX;
            
            %NOTE: I later decided I didn't want to do cutoffs in the
            %function so we set the bounds outside the working range
            MAX_CUTOFF = -360;
            MIN_CUTOFF = 360;
            
            pixels_local = obj.pixels;
            avg_worm_segment_length = obj.cc_lengths(end)/obj.N_SEGS;
            
            hf_angle_edge_length = avg_worm_segment_length;
            lf_angle_edge_length = 2 * hf_angle_edge_length;
            
            obj.hf_angles_raw = seg_worm.cv.circCurvature(pixels_local, hf_angle_edge_length, obj.cc_lengths);
            obj.lf_angles     = seg_worm.cv.circCurvature(pixels_local, lf_angle_edge_length, obj.cc_lengths);
            
            filter_width  = ceil(avg_worm_segment_length/2);
            obj.hf_angles = seg_worm.util.circConv(obj.hf_angles_raw ,[],filter_width);
            
            FH = @seg_worm.util.peaksCircDist;
            [obj.lf_ap_max,obj.lf_ap_max_I] = FH(obj.lf_angles,lf_angle_edge_length,USE_MAX,MAX_CUTOFF,obj.cc_lengths);
            [obj.lf_ap_min,obj.lf_ap_min_I] = FH(obj.lf_angles,lf_angle_edge_length,USE_MIN,MIN_CUTOFF,obj.cc_lengths);
            [obj.hf_ap_max,obj.hf_ap_max_I] = FH(obj.hf_angles,hf_angle_edge_length,USE_MAX,MAX_CUTOFF,obj.cc_lengths);
            [obj.hf_ap_min,obj.hf_ap_min_I] = FH(obj.hf_angles,hf_angle_edge_length,USE_MIN,MIN_CUTOFF,obj.cc_lengths);
            
        end
        function computeCCLengths(obj)
            obj.cc_lengths = seg_worm.cv.circComputeChainCodeLengths(obj.pixels);
        end
        function plot(obj)
            %TODO: Add a third image which shows the subtraction
            
            pixels_local = obj.pixels;
            cMinY   = min(pixels_local(:,1));
            cMaxY   = max(pixels_local(:,1));
            cMinX   = min(pixels_local(:,2));
            cMaxX   = max(pixels_local(:,2));
            cHeight = cMaxY - cMinY + 1;
            cWidth  = cMaxX - cMinX + 1;
            cImg    = zeros(cHeight, cWidth);
            cImg(sub2ind(size(cImg), ...
                pixels_local(:,1) - cMinY + 1, ...
                pixels_local(:,2) - cMinX + 1)) = 255;
            imshow(cImg)
        end
    end
    
end

