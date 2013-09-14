classdef orientation < sl.obj.handle_light
    %
    %   Class:
    %   seg_worm.worm.orientation
    
    %NOT YET IMPLEMENTED
    %See comments below for content that might go in this class ...
    
    properties
       ORIENTATION_PCT = [1:5 7:11] / 12; %???? - why these values? 
       HEAD_TAIL_SCALE = 134217728; % 2^26
       VULVA_SCALE     = 1048576;   % 2^20
    end
    
    properties
       parent %seg_worm.worm
       
       %Head/tail
       %----------------------------------------------
       head_and_tail_flipped = false
       head_confidence
       tail_confidence
       
       %Vulva
       %-----------------------------------------------
       vulva_clockwise_from_head = false
       vulva_confidence
       non_vulva_confidence
       
       %EDIT TO BELOW: I'm not actually sure what these variables
       %represent at this point ...
       %
       %    
       %
       %??? Why don't we normalize the confidence levels????
       %NOTE: This would make the scales irrelevant
    end
    
    methods
        function obj = orientation(parent)
           obj.parent = parent;
           obj.populateHeadTailConfidence();
           obj.populateVulvaConfidence();
        end
        function populateHeadTailConfidence(obj)
            
            worm = obj.parent;
            
            lfCAngles = worm.contour.lf_angles;
            
            cdfs   = {worm.head.pixel_cdf       worm.tail.pixel_cdf};
            Is     = [worm.contour.head_I       worm.contour.tail_I];
            stdevs = [worm.head.pixel_std_dev   worm.tail.pixel_std_dev];
            
            % How much confidence do we have in our head-to-tail orientation?
            % Note: generally, the head is less angled, and contains more white
            % pixels (a higher 50% and 75% CDF for color) and less gray pixels (a higher
            % variance and 25% to 75% interquartile range) than the tail. We give
            % each probability equal weight, then compare.
            
            %JAH NOTE: I'm not sure I understand this ...
            
            for iConfidence = 1:2
               cur_cdf = cdfs{iConfidence};
               cur_I   = Is(iConfidence);
               stdev   = stdevs(iConfidence);
               
               a = 180 - lfCAngles(cur_I);
               b = cur_cdf(3)*cur_cdf(4)*stdev*(cur_cdf(4)-cur_cdf(2));
               
               c_value = a*b/obj.HEAD_TAIL_SCALE;
                
               if iConfidence == 1
                   obj.head_confidence = c_value;
               else
                   obj.tail_confidence = c_value;
               end
            end
        end
        function populateVulvaConfidence(obj)           
            worm = obj.parent;

            cdfs   = {worm.right_side.pixel_cdf       worm.left_side.pixel_cdf};
            stdevs = [worm.right_side.pixel_std_dev   worm.left_side.pixel_std_dev];

            % How much confidence do we have in our vulva orientation?
            % Note: generally, the vulval side contains less white pixels (a lower
            % 50% and 75% CDF for color) and more gray pixels (a lower variance and
            % 25% to 75% interquartile range) than the opposing side. We give each
            % probability equal weight, then compare. Also, in the absence of
            % information, we assume the vulva is on the left side (and use a trick
            % to avoid reciprocals in our equations).
            for iConfidence = 1:2
               cur_cdf = cdfs{iConfidence};
               stdev   = stdevs(iConfidence);
               
               b = cur_cdf(3)*cur_cdf(4)*stdev*(cur_cdf(4)-cur_cdf(2));
               
               c_value = b/obj.VULVA_SCALE;
                
               if iConfidence == 2
                   obj.vulva_confidence = c_value;
               else
                   obj.non_vulva_confidence = c_value;
               end
            end
        end
    end
    
    %OLD FUNCTIONS     ====================================================
    methods
        function isFlipped = isWormCellHeadFlipped(worm)
        %   Input:
        %       worm - the worm information organized in a cell array
        %
        %   Output:
        %       isFlipped - are the worm's head and tail flipped?

        isFlipped = worm{8}{1}{1};
        end
    end
    
    %HELPER FUNCTIONS    ==================================================
    methods (Hidden)
       function points = ohelper__getPointsAlongSkeleton(obj,pct_samples)

           worm       = obj.parent;
           cc_lengths = worm.skeleton.cc_lengths;
           is_flipped = obj.head_and_tail_flipped;
           
           % Sample worm 1.
           if is_flipped
               flippedSamples = 1 - pct_samples;
               s1 = cc_lengths * flippedSamples;
           else
               s1 = cc_lengths * pct_samples;
           end
           points = seg_worm.cv.chainCodeLengthInterp(skeleton1.pixels, s1, skeleton1.chainCodeLengths);
       end
    end
    
end




%{
headConfidence = struct('head', hConfidence, 'tail', tConfidence);
headOrientation = struct('isFlipped', isHeadTailFlipped, 'confidence', headConfidence);
vulvaConfidence = struct('vulva', vConfidence, 'nonVulva', nvConfidence);
vulvaOrientation = struct('isClockwiseFromHead', isVulvaClockwiseFromHead, 'confidence', vulvaConfidence);
orientation = struct('head', headOrientation, 'vulva', vulvaOrientation);
%}
