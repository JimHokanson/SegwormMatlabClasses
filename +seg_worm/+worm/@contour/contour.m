classdef contour < handle
    %
    %   Class:
    %   seg_worm.worm.contour
    
    properties
        parse_error = false
        error_num = 0
        error_msg = ''
    end
    
    properties (Constant)
        N_SEGS = 48 %NOTE: This should be 2x
    end
    
    properties
        image_used
        rough_pixels %(cPixels)[n x 2] the worm's circularly continuous contour pixels,
        %ordered clockwise
        %
        %see: .initialize()
        cleaned_pixels %(cPixels)[n x 2] rough pixels after running 
        %cleaning algorithm
        %
        %see: .cleanContour()
        
        
        
        touch_points_I  %(cTouchI) the paired pairs of indices marking, clockwise, the
        % start and end of the touching contour points Note: if the worm isn't
        % coiled, this value is empty.
        inner_I %(cInI) the paired indices marking, clockwise, the start and
        % end of the inner contour points
        outer_I %(cOutI) the paired indices marking, clockwise, the start and
        % end of the outer contour points
        lf_angles  %(cAngles) the contour's angles (curvature) at each index
        hf_angles_raw
        hf_angles_smoothed
        
        worm_segment_length
        hf_edge_length
        lf_edge_length
        
        
        head_I  %(cHeadI)  the contour index for the worm's head
        tail_I  %(cTailI)  the contour index for the worm's tail
        cc_lengths %(cCCLengths) - the contour's circular chain-code pixel length, from
        % its vector's start to end, up to each contour point Note: this is a more
        % accurate representation of locations along the worm's contour than pixel
        % indices
    end
    
    methods
        function obj = contour(img,frame_number,verbose,num_erode,num_dilate,is_normalized)
            %seg_worm.worm.contour(img,frame_number,verbose,num_erode,num_dilate,is_normalized)
            obj.initialize(img,frame_number,verbose,num_erode,num_dilate,is_normalized)
        end
        function plot(obj)
            %TODO: Add a third image which shows the subtraction
            
            for iType = 1:2
                
                if iType == 1
                    pixels_local = obj.rough_pixels;
                else
                    pixels_local = obj.cleaned_pixels;
                end
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
                subplot(1,2,iType)
                imshow(cImg)
            end
        end
    end
    
end

