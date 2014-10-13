function checkIfWormIsCoiled(obj)
%
%
%https://github.com/openworm/SegWorm/blob/master/Worms/Segmentation/segWorm.m#L579

contour_obj = obj.contour;

BEND_ANGLE_CUTOFF = -30;

% Is the worm coiled?
%--------------------------------------------------------------------------
% If there are no large concavities, the worm is not coiled.
lf_contour_bend_I = contour_obj.lf_ap_min_I(contour_obj.lf_ap_min < BEND_ANGLE_CUTOFF);

%Possible early return ...
%-------------------------------------------------------
if isempty(lf_contour_bend_I)
    return
end

error_handler     = obj.error_handler;
skeleton_c_widths = obj.skeleton.c_widths;
max_c_width = max(skeleton_c_widths);


head   = obj.head;
hBendI = helper__getBendsNearHeadOrTail(lf_contour_bend_I,head.left_contour_bounds,head.right_contour_bounds,true);
error_handler.wormTooWideBasedOnHeadWidth(hBendI,max_c_width,obj.head.innermost_contour_width);
if error_handler.error_found, return; end;

tail   = obj.tail;
tBendI = helper__getBendsNearHeadOrTail(lf_contour_bend_I,tail.left_contour_bounds,tail.right_contour_bounds,false);
error_handler.wormTooWideBasedOnTailWidth(tBendI,max_c_width,obj.tail.innermost_contour_width);
if error_handler.error_found, return; end;


lfCAngles = contour_obj.lf_angles;

% Use the most accurate estimate of head/tail width to
% determine whether the width of the body is more than double
% that at the end of the head/tail; in which case; the worm is
% coiled.
if ~(isempty(hBendI) && isempty(tBendI))
    
    % Find the distances of bends near the head.
    %----------------------------------------------------------------------
    headI = obj.contour.head_I;
    hBendDist = abs(headI - hBendI);
    %???  - what is this doing ????
    hBendDist = min(hBendDist, abs(hBendDist - length(lfCAngles)));
    
    % Find the distances of bends near the tail.
    %----------------------------------------------------------------------
    tailI = obj.contour.tail_I;
    tBendDist = abs(tailI - tBendI);
    %???  - what is this doing ????
    tBendDist = min(tBendDist, abs(tBendDist - length(lfCAngles)));
    
    %Why do we need another min, aren't the above values already
    %scalars????
    if min(hBendDist) >= min(tBendDist)
        % The bend near the head is furthest and, therefore, the
        % width at the end of the head is our most accurate
        % estimate of the worm's width.
        error_handler.wormTooWideBasedOnHeadWidth(hBendI,max_c_width,obj.head.innermost_contour_width);
    else
        % The bend near the tail is furthest and, therefore, the
        % width at the end of the tail is our most accurate
        % estimate of the worm's width.
        error_handler.wormTooWideBasedOnTailWidth(tBendI,max_c_width,obj.tail.innermost_contour_width);
    end
end


end

function ht_bend_I = helper__getBendsNearHeadOrTail(lf_contour_bend_I,left_bounds,right_bounds,is_head) %#ok<INUSD>

%Not sure why the difference in handling between head and tail
%This involves getting a better understanding of the worm2poly() function

%This approach (as opposed to the old code) seems better IF the inputs are sorted ...
%We'll do a check to make sure this is the case ...
%i.e. that value at index 2 is greater than the value at index 1
if left_bounds(1) < left_bounds(2) && right_bounds(1) < right_bounds(2)
    mask = (lf_contour_bend_I >= left_bounds(1) & lf_contour_bend_I <= left_bounds(2)) | ...
        (lf_contour_bend_I >= right_bounds(1) & lf_contour_bend_I <= right_bounds(2));
    ht_bend_I = lf_contour_bend_I(mask);
else
    %If this is not true it would indicate a circular pattern
    %in which case the code below might be better
    error('Assumption violated')
end

%{
%OLD CODE
%AAAAAAAAH
%------------------------------------------------------------------


%hhhhhh      <-    tttttt <-
%START OF CONTOUR         / \
%->hhhhhh            tttttt|


%NOTE: Earlier code ensured the head is at 1, but the original
%code tried to make sure to handle cases in which this was not the
%case. We'll mimic the original code ....

%FOR THE HEAD

% Find concavities near the head. If there are any concavities
% near the tail, the head may be portruding from a coil; in
% which case, the width at the end of the head may be
% inaccurate.
if is_head
    hlcBounds = left_bounds;
    hrcBounds = right_bounds;
    if hlcBounds(1) < hrcBounds(2)
        ht_bend_I = lf_contour_bend_I(lf_contour_bend_I > hlcBounds(1) & lf_contour_bend_I < hrcBounds(2));
    else
        ht_bend_I = lf_contour_bend_I(lf_contour_bend_I > hlcBounds(1) | lf_contour_bend_I < hrcBounds(2));
    end
    
else
    %FOR THE TAIL
    
    % Find concavities near the tail. If there are any concavities near
    % the tail, the tail may be portruding from a coil; in which case,
    % the width at the end of the tail may be inaccurate.
    
    trcBounds = right_bounds;
    tlcBounds = left_bounds;
    if trcBounds(1) < tlcBounds(2)
        ht_bend_I = lf_contour_bend_I(lf_contour_bend_I > trcBounds(1) & lf_contour_bend_I < tlcBounds(2));
    else
        ht_bend_I = lf_contour_bend_I(lf_contour_bend_I > trcBounds(1) | lf_contour_bend_I < tlcBounds(2));
    end
    
end
%}

end