classdef parse_error < handle
    %
    %   Class:
    %   seg_worm.parse_error
    
    properties
        error_found
        error_number
        error_message
        original_image
        frame_number
        verbose
    end
    
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
            errNum = 0;
            errMsg = '';
            if size(contour, 1) < C_WORM_SEGS
                errNum = 103;
                errMsg = 'The worm contour is too small.';
                warnID = 'segWorm:ContourTooSmall';
                helper__showOrigAndNewVerbose(obj,warnID,errMsg,img);
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
