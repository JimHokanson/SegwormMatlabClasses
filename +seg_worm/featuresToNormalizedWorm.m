function featuresToNormalizedWorm(feature_file_path,normalized_worm)

h = load(feature_file_path);

keyboard

%This is difficult

segmentation_status  %[1 n],char
%    s = segmented
%    f = segmentation failed
%    m = stage movement
%    d = dropped frame
%    n??? - there is reference in some old code to this which 
        
        
        
        frame_codes         %see comments in seg_worm.parsing.frame_errors
        %near the bottom, I haven't yet coded in the values as constants
        %... :/
        %ex.
        %1    = ok
        %2    = stage movement
        %101  = no worm found
        %
        %   NOTE: This is needed for processing coils ...
        %   
        vulva_contours       %[49 2 n] double
        non_vulva_contours   %[49 2 n] double
        skeletons            %[49 2 n] double
        angles               %[49 n] double, degrees
        %??? Why are angles at edge undefined ??????
        %??? - NaN values ...
        %First and last 5 of 49 values ...
        in_out_touches       %[49 n] double
        lengths              %[1  n] double
        widths               %[49 n] double
        head_areas           %[1  n] double
        tail_areas           %[1  n] double
        vulva_areas          %[1  n] double
        non_vulva_areas      %[1  n] double

end