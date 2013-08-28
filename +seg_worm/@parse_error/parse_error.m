classdef parse_error < handle
    %UNTITLED10 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
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
        function [errNum,errMsg] = noWormFound(obj,worm_pixels)
            % No worm found.
            if isempty(worm_pixels)
                errNum = 101;
                errMsg = 'No worm was found.';
                
                % Show the failure.
                if obj.verbose
                    warning('segWorm:NoWormFound', ['Frame %d: ' errMsg],obj.frame_number);
                    
                    % Open a big figure.
                    figure('OuterPosition', [50 50 1280 960]);
                    set(gcf, 'Color', [1 .5 .5]);
                    
                    % Show the original image.
                    imshow(obj.original_image);
                    title('Original Image');
                end
            else
                errNum = [];
                errMsg = [];
            end
        end
        function [errNum,errMsg] = contourTouchesBoundary(obj,contour,img)
            % The contour touches a boundary.
            errNum = [];
            errMsg = [];
            [m,n] = size(obj.img);
            if min(contour(:,1)) == 1 || min(contour(:,2)) == 1 || ...
                    max(contour(:,1)) == m || max(contour(:,2)) == n
                errNum = 102;
                errMsg = 'The worm contour touches the image boundary.';
                warnID = 'segWorm:ContourTouchesBoundary';
                helper__showOrigAndNewVerbose(obj,warnID,img);
            end
        end
        function [errNum,errMsg] = contourTooSmall(obj,contour,C_WORM_SEGS,img)
            % The contour is too small.
            errNum = [];
            errMsg = [];
            if size(contour, 1) < C_WORM_SEGS
                errNum = 103;
                errMsg = 'The worm contour is too small.';
                warnID = 'segWorm:ContourTooSmall';
                helper__showOrigAndNewVerbose(obj,warnID,img);
            end
        end
    end
    
end

function helper__showOrigAndNewVerbose(obj,warnID,img)
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
