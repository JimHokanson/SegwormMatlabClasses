classdef frame_errors < sl.obj.handle_light
    %
    %   Class:
    %   seg_worm.parsing.frame_errors
    %
    %   There is a file that ends with the extension:
    %
    %   _failedFrames
    %
    %   It contains the variable 'failedFrames'
    %
    %   Column 1: frame number
    %   Column 2: warning id
    %
    %
    %   NOTE: I don't think I am going to use this to hold values
    %   instead I think I am just going to hold the numeric value and
    %   build an interface 
    
    properties
      %TODO: Replace with IDs ???
       %frame_numbers
       %id
    end
    
    %Segmentation Status:
    %----------------------------------------------------
    %    s = segmented
    %    f = segmentation failed
    %    m = stage movement
    %    d = dropped frame
    %    n??? - there is reference in some old code to this which 
    
    
    methods (Static)
        function error_codes = segmentationStatusToCodes(seg_status,failed_frame_info)
            %
            %   seg_worm.parsing.frame_errors.segmentationStatusToCodes
            %
            %   failed_frames: [n 2], col 1: frame #, col 2: error code
            %   - currently saved to disk by GUI
            
            ERROR_CODE_OFFSET = 100; %When the failed_frame_info is saved
            %to disk 100 is subtracted for some reason ...
            %We readd it here ...
            
            error_codes = zeros(size(seg_status));
            
            %Handling failed frames
            %--------------------------------------------------------------
            failed_frame_numbers = failed_frame_info(:,1);
            failed_frame_errors  = failed_frame_info(:,2);
            
            failed_I = find(seg_status == 'f');
            
            if ~isequal(failed_I(:),sort(failed_frame_numbers))
               error('Mismatch between saved failed frames and the segmentation status indicating a failed frame') 
            end
            error_codes(failed_frame_numbers) = failed_frame_errors + ERROR_CODE_OFFSET;
            
            %Set the rest ...
            %--------------------------------------------------------------
            mask_char  = 'smd';
            mask_value = [1 2 3];
            
            for iVar = 1:3
               error_codes(seg_status == mask_char(iVar)) = mask_value(iVar);
            end

            if any(error_codes == 0)
                %unique(seg_status(error_codes == 0))
                error('Some frame error codes haven''t been fixed')
            end
            
        end
    end
    
end

%{
    { 1 1       'segWorm:Success'                   'The worm was successfully segmented.' ...
      2 2       'findStageMovement:StageMovement'   'The video frame contains stage motion.' ...
      3 3       'segWorm:DroppedFrame'              'The video frame was dropped.'
      4 101     'segWorm:NoWormFound'               'No worm was found in the video frame.'
      5 102     'segWorm:ContourTouchesBoundary'    'The worm contour touches the image boundary.'
      6 103     'segWorm:ContourTooSmall'           'The worm contour is too small.'
      7 104     'segWorm:TooManyEnds'               'The worm has 3 or more low-frequency sampled convexities sharper than 90 degrees (possible head/tail points).'
      8 105, 	'segWorm:TooFewEnds'                'The worm contour has less than 2 high-frequency sampled convexities sharper than 60 degrees (the head and tail). Therefore, the worm is coiled or obscured and cannot be segmented.', ...
      9 106,    'segWorm:DoubleLengthSide'          'The worm length, from head to tail, is more than twice as large on one side than it is on the other. Therefore, the worm is coiled or obscured and cannot be segmented.', ...
      10 107,   'segWorm:DoubleHeadWidth'           'The worm more than doubles its width from the end of its head. Therefore, the worm is coiled, laid an egg, and/or is significantly obscured and cannot be segmented'], ...
      11 108,   'segWorm:DoubleTailWidth'           'The worm more than doubles its width from the end of its tail. Therefore, the worm is coiled, laid an egg, and/or is significantly obscured and cannot be segmented.', ...
      12 109,   'segWorm:SmallTail'                 'The worm tail is less than half the size of its head. Therefore, the worm is significantly obscured and cannot be segmented.'
      13 110, 	'segWorm:SmallHead'                 'The worm head is less than half the size of its tail. Therefore, the worm is significantly obscured and cannot be segmented.'
      14 111,   'segWorm:SmallHeadTail'             'The worm head and tail are less than 1/4 the size of its remaining body. Therefore, the worm is significantly obscured and cannot be segmented.'
      15 1001   'normWorms:TooShort'                'The worm is shorter than the sampling points requested.'});

%}