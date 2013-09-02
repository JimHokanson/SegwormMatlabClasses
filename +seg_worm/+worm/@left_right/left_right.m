classdef left_right < sl.obj.handle_light
    %
    %   Class:
    %   seg_worm.worm.left_right
    
    properties
       parent
       contour
       skeleton
    end
    
    properties (Abstract)
       type
    end
    properties (Hidden)
       is_right 
    end
    methods
        function value = get.is_right(obj)
            value = strcmp(obj.type,'right');
        end
    end
    
    
    methods (Static)
        function [left,right] = createSides(parent)
            %
            %   [left,right] = seg_worm.worm.left_right.createSides(parent)
           left  = seg_worm.worm.left;
           right = seg_worm.worm.right;
           left.initRefs(parent);
           right.initRefs(parent);
           %left.initBounds(right);
        
           head = parent.head;
           tail = parent.tail;
           %tsBounds and hsBounds
           %.skeleton_bounds()
           keyboard
        
           %{
            [sides,lcBounds,rcBounds,sBounds] = worm2poly(hsBounds(2), ...
            tsBounds(1), sPixels, headI, tailI, cPixels, true, ...
            sCCLengths, cCCLengths); 
            %}
        
        end
    end
    methods
        function initBounds(obj)
            

        end
        function initRefs(obj,parent)
            obj.parent   = parent;
            obj.contour  = parent.contour;
            obj.skeleton = parent.skeleton;
        end
    end
    
end

%{

% Determine the left-side's MER (minimum enclosing rectangle).
[sides,lcBounds,rcBounds,sBounds] = worm2poly(hsBounds(2), ...
    tsBounds(1), sPixels, headI, tailI, cPixels, true, ...
    sCCLengths, cCCLengths);
lSide = sides{2};
lMinY = min(lSide(:,1));
lMaxY = max(lSide(:,1));
lMinX = min(lSide(:,2));
lMaxX = max(lSide(:,2));

% Measure the left side (counter clockwise from the head) statistics.
%lCDF = [];
%lStdev = [];
merLImg       = oImg(lMinY:lMaxY, lMinX:lMaxX);
merLSide      = [lSide(:,1) - lMinY + 1, lSide(:,2) - lMinX + 1];
[merLMask, ~] = inPolyMask(merLSide, [], size(merLImg));
lColors       = single(merLImg(merLMask));
lArea         = length(lColors);
lCDF          = prctile(lColors,[2.5 25 50 75 97.5]);
lStdev        = std(lColors);

% Determine the right-side's MER (minimum enclosing rectangle).
rSide = sides{1};
rMinY = min(rSide(:,1));
rMaxY = max(rSide(:,1));
rMinX = min(rSide(:,2));
rMaxX = max(rSide(:,2));


% Measure the right side (clockwise from the head) statistics.
%rCDF = [];
%rStdev = [];
merRImg       = oImg(rMinY:rMaxY, rMinX:rMaxX);
merRSide      = [rSide(:,1) - rMinY + 1, rSide(:,2) - rMinX + 1];
[merRMask, ~] = inPolyMask(merRSide, [], size(merRImg));
rColors       = single(merRImg(merRMask));
rArea         = length(rColors);
rCDF          = prctile(rColors,[2.5 25 50 75 97.5]);
rStdev        = std(rColors);



%}