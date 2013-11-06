classdef orientation < sl.obj.handle_light
    %
    %   Class:
    %   seg_worm.worm.orientation
    
    %NOT YET IMPLEMENTED
    %See comments below for content that might go in this class ...
    
    %NOTE: This is really tricky as we need to be pretty explicit
    %as to when these variables are evaluated ...
    
    
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
        head_confidence %confidence that current head is really the head
        tail_confidence %confidence that current tail is really the head
        
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
        %ohelper - orientation helper, only really meant to be used
        %internally, thus being hidden
        function points = ohelper__getPointsAlongSkeleton(obj,pct_samples)
            %
            %   
            %
            
            worm       = obj.parent;
            skeleton   = worm.skeleton;
            pixels     = skeleton.pixels;
            total_skeleton_length = skeleton.total_length;
            cc_lengths = skeleton.cc_lengths;
            is_flipped = obj.head_and_tail_flipped;
            
            % Sample worm 1.
            if is_flipped
                flippedSamples = 1 - pct_samples;
                s1 = total_skeleton_length * flippedSamples;
            else
                s1 = total_skeleton_length * pct_samples;
            end
            points = seg_worm.cv.chainCodeLengthInterp(pixels, s1, cc_lengths);
        end
        function ohelper__flipWormToMatchRef(cur_orientation,pRef,pCur)
            %
            %
            %   pRef : [n x 2], points along reference skeleton
            %   pCur : [n x 2], points along current skeleton
            %
            %   See Also:
            %   orientToRefWorm
            %
            
            % Compute the distance between worm 1 and 2 in both orientations.
            % Note: the skeleton samples are interpolated and, therefore, distances
            % less than 1 are likely to be, in reality, either 1 or 0. We set these to
            % 0 to avoid computing extreme confidences.
            proximity_same = sqrt(sum((pRef - pCur) .^ 2, 2));
            proximity_flip = sqrt(sum((pRef - pCur(end:-1:1,:)) .^ 2, 2));
            
            %??? - why is this done ????
            %Does this basically just say we can't really resolve
            %differences less than 1?
            %I think so ...
            proximity_same(proximity_same < 1) = 0;
            proximity_flip(proximity_flip < 1) = 0;
            
            % Compute the indicators for both orientations.
            %NOTE: In general we want the points to be close
            %between reference and current, so being less than the other
            %version is better, or is indicative of being that type
            indicators_same    = sum(proximity_same < proximity_flip);
            indicators_flipped = sum(proximity_same > proximity_flip);
            
            if indicators_flipped > indicators_same
                is_flipped = true;
            elseif indicators_same > indicators_flipped
                is_flipped = false;
            else
                
                %JAH NOTE: I haven't done anything with this bit of code
                %...
                
                % Compute the confidence for both orientations.
                % Note: to avoid dividing by 0, when the magnitude's distance denominator
                % is 0 we substitute the numerator's distance instead of the ratio.
                magnitudes1     = proximity_flip ./ proximity_same;
                zeroConfidence  = (proximity_same < 1);
                magnitudes1(zeroConfidence) = proximity_flip(zeroConfidence);
                confidence1     = mean(magnitudes1);
                magnitudes2     = proximity_same ./ proximity_flip;
                zeroConfidence  = (proximity_flip < 1);
                magnitudes2(zeroConfidence) = proximity_same(zeroConfidence);
                confidence2     = mean(magnitudes2);
                
                is_flipped = confidence1 >= confidence2;
            end
            
            if is_flipped
                cur_orientation.flipWormHead();
            end            
        end
    end
    methods
        function flipWormHead(obj)
           
            
            %TODO: addManipulation(obj,code,description,details)
            
            obj.head_and_tail_flipped = ~obj.head_and_tail_flipped;
           [obj.head_confidence,obj.tail_confidence] = deal(...
               obj.tail_confidence,obj.head_confidence);
           
           %???? - why isn't a call made to the flipWormVulva code???
           %
           %
           %Doesn't this change the confidences?????
           %Note: since the vulva is specified relative to the head, its location
           %flips to preserve its orientation.
           obj.vulva_clockwise_from_head = ~obj.vulva_clockwise_from_head;
        end
        function flipWormVulva(obj)
            
           %TODO: addManipulation(obj,code,description,details) 
            
            
           obj.vulva_clockwise_from_head = ~obj.vulva_clockwise_from_head;
           [obj.vulva_confidence,obj.non_vulva_confidence] = deal(...
               obj.non_vulva_confidence,obj.vulva_confidence);           
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
