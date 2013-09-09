classdef parse_error < handle
    %
    %   Class:
    %   seg_worm.parse_error
    
    properties
        error_found   = false
        error_number  = 0
        error_message = ''
        
        %TODO: We might want to make this a link to the
        %seg_worm.worm object to reduce memory usage, especially
        %if we are going to clear it later when hold onto all of the worms
        %...
        original_image 
        frame_number
        verbose
    end
    
    %Contour ==============================================================
    methods
        function obj = parse_error(original_image,frame_number,verbose)
            obj.original_image = original_image;
            obj.frame_number = frame_number;
            obj.verbose = verbose;
        end
        function noWormFound(obj,worm_pixels)
            % No worm found.
            if isempty(worm_pixels)
                obj.error_found  = true;
                obj.error_number = 101;
                obj.error_message = 'No worm was found.';
                
                % Show the failure.
                if obj.verbose
                    warning('segWorm:NoWormFound', ['Frame %d: ' obj.error_message],obj.frame_number);
                    
                    % Open a big figure.
                    figure('OuterPosition', [50 50 1280 960]);
                    set(gcf, 'Color', [1 .5 .5]);
                    
                    % Show the original image.
                    imshow(obj.original_image);
                    title('Original Image');
                end
            end
        end
        function contourTouchesBoundary(obj,contour_pixels,img)
            % The contour touches a boundary.
            [m,n] = size(obj.original_image);
            if min(contour_pixels(:,1)) == 1 || min(contour_pixels(:,2)) == 1 || ...
                    max(contour_pixels(:,1)) == m || max(contour_pixels(:,2)) == n
                obj.error_found  = true;
                obj.error_number = 102;
                obj.error_message = 'The worm contour touches the image boundary.';
                warnID = 'segWorm:ContourTouchesBoundary';
                helper__showOrigAndNewVerbose(obj,warnID,obj.error_message,img);
            end
        end
        function contourTooSmall(obj,contour,C_WORM_SEGS,img)
            % The contour is too small.
            if size(contour, 1) < C_WORM_SEGS
                obj.error_found   = true;
                obj.error_number = 103;
                obj.error_message = 'The worm contour is too small.';
                warnID = 'segWorm:ContourTooSmall';
                helper__showOrigAndNewVerbose(obj,warnID,obj.error_message,img);
            end
        end
        function insufficientHTOptions(obj,n_lf_HT_I,n_hf_HT_I)
            if n_lf_HT_I > 2
                obj.error_found   = true;
                obj.error_number  = 104;
                obj.error_message = ['The worm has 3 or more low-frequency sampled convexities ' ...
                    'sharper than 90 degrees (possible head/tail points).'];
                return
            end
            
            if n_hf_HT_I < 2
                obj.error_found   = true;
                obj.error_number  = 105;
                obj.error_message = ['The worm contour has less than 2 high-frequency sampled '...
                    'convexities sharper than 60 degrees (the head and tail). ' ...
                    'Therefore, the worm is coiled or obscured and cannot be segmented.'];
                return
            end
        end
        function lopsidedSides(obj,c_obj)
            
            %TODO: Add verbose ...
            
            cc_lengths_local = c_obj.cc_lengths;
            head_I = c_obj.head_I;
            tail_I = c_obj.tail_I;
            
            % Find the length of each side.
            %--------------------------------------------------------------------------
            if head_I > tail_I
                size_1 = cc_lengths_local(head_I) - cc_lengths_local(tail_I);
                size_2 = cc_lengths_local(end)    - cc_lengths_local(head_I) + cc_lengths_local(tail_I);
            else
                size_1 = cc_lengths_local(tail_I) - cc_lengths_local(head_I);
                size_2 = cc_lengths_local(end)    - cc_lengths_local(tail_I) + cc_lengths_local(head_I);
            end
            
            % Are the sides within 50% of each others size?
            % Note: if a worm's length from head to tail is at least twice larger
            % on one side (relative to the other), than the worm must be touching
            % itself.
            if min(size_1, size_2)/ max(size_1, size_2) <= .5
                obj.error_found   = true;
                obj.error_number  = 106;
                obj.error_message = ['The worm length, from head to tail, is more than ' ...
                    'twice as large on one side than it is on the other. ' ...
                    'Therefore, the worm is coiled or obscured and cannot be segmented.'];
            end
        end
    end
    
    %Head/Tail/Left/Right ============================================================
    methods
        function checkHeadTailArea(obj,worm_obj)
            
            hArea = worm_obj.head.pixel_area;
            tArea = worm_obj.tail.pixel_area;
            
            % Is the tail too small (or the head too large)?
            % Note: the area of the head and tail should be roughly the same size.
            % A 2-fold difference is huge!
            if hArea > 2 * tArea
                obj.error_found  = true;
                obj.error_number = 109;
                obj.error_message = ['The worm tail is less than half the size of its ' ...
                    'head. Therefore, the worm is significantly obscured and ' ...
                    'cannot be segmented.'];
                
                % Defer organizing the available worm information.
                if obj.verbose
                    warning('segWorm:SmallTail', 'Frame %d: %s',obj.frame_number,obj.error_message);
                end
            elseif tArea > 2 * hArea
                obj.error_found  = true;
                obj.error_number = 110;
                obj.error_message = ['The worm head is less than half the size of its ' ...
                    'tail. Therefore, the worm is significantly obscured and ' ...
                    'cannot be segmented.'];
                
                % Defer organizing the available worm information.
                if obj.verbose
                    warning('segWorm:SmallHead','Frame %d: %s',obj.frame_number,obj.error_message);
                else
                    return;
                end
            end
        end
        function headTailSmallOrBodyLarge(obj,worm_obj)
            
            hArea = worm_obj.head.pixel_area;
            tArea = worm_obj.tail.pixel_area;
            lArea = worm_obj.left_side.pixel_area;
            rArea = worm_obj.right_side.pixel_area;
            
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
                obj.error_found  = true;
                obj.error_number = 111;
                obj.error_message = ['The worm head and tail are less than 1/4 the size ' ...
                    'of its remaining body. Therefore, the worm is ' ...
                    'significantly obscured and cannot be segmented.'];
                
                % Defer organizing the available worm information.
                if obj.verbose
                    warning('segWorm:SmallHeadTail', ['Frame %d: ' obj.error_message], obj.frame_number);
                end
            end
            
        end
    end
    %Coiled Checks  =======================================================
    methods
        function wormTooWideBasedOnHeadWidth(obj,hBendI,max_c_width,innermost_head_contour_width)
            % Does the worm more than double its width from the head?
% Note: if the worm coils, its width will grow to more than
% double that at the end of the head.
            if isempty(hBendI) && max_c_width > 2*innermost_head_contour_width
                obj.error_found  = true;
                obj.error_number = 107;
                obj.error_message = ['The worm more than doubles its width ' ...
                    'from end of its head. Therefore, the worm is ' ...
                    'coiled, laid an egg, and/or is significantly ' ...
                    'obscured and cannot be segmented.'];
                
                % Organize the available worm information.
                if obj.verbose
                    warning('segWorm:DoubleHeadWidth', ...
                        ['Frame %d: ' obj.error_message], obj.frame_number);
                end
            end
        end
        function wormTooWideBasedOnTailWidth(obj,tBendI,max_c_width,innermost_tail_contour_width)
            % Does the worm more than double its width from the tail?
% If the worm coils, its width will grow to more than double
% that at the end of the tail.
            if isempty(tBendI) && max_c_width > 2*innermost_tail_contour_width
                obj.error_found  = true;
                obj.error_number = 108;
                obj.error_message = ['The worm more than doubles its width ' ...
                    'from end of its tail. Therefore, the worm is ' ...
                    'coiled, laid an egg, and/or is significantly ' ...
                    'obscured and cannot be segmented.'];
                
                % Organize the available worm information.
                if obj.verbose
                    warning('segWorm:DoubleTailWidth', ...
                        ['Frame %d: ' obj.error_message], obj.frame_number);
                end
            end
 
        end
    end
    
end

function helper__showOrigAndNewVerbose(obj,warnID,errMsg,img)
% Show the failure.
if obj.verbose
    warning(warnID, ['Frame %d: ' errMsg],obj.frame_number);
    
    % Open a big figure.
    figure('OuterPosition', [50 50 1280 960]);
    set(gcf, 'Color', [1 .5 .5]);
    
    % Show the original image.
    subplot(1,2,1), imshow(obj.original_image);
    title('Original Image');
    
    %??? - This might change
    % Show the thresholded image.
    subplot(1,2,2), imshow(img);
    title('Thresholded Image');
end
end
