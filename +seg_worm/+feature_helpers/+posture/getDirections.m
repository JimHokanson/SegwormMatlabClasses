function directions = getDirections(sx,sy)
%
%
%   seg_worm.feature_helpers.posture.getDirections
%
%   INPUTS
%   ===========================================================
%   skeleton_x : [49 x n_frames]
%   skeleton_y : [49 x n_frames]
%
%   skeleton_XY : [49 (x y) frames]
%
%   OUTPUTS
%   ===========================================================
%   directions : 
%       .tail2head 
%       .head
%       .tail
%  
%   Old Name: featureProcess.m
%
%   Nature Methods Description
%   ===========================================================
%   Orientation. 
%   ------------------
%   The worm’s orientation is measured overall (from tail to head) as well
%   as for the head and tail individually. The overall orientation is
%   measured as the angular direction from the tail to the head centroid.
%   The head and tail centroids are computed as the mean of their
%   respective skeleton points. The head and tail direction are computed by
%   splitting these regions in two, then computing the centroid of each
%   half. The head direction is measured as the angular direction from the
%   its second half (the centroid of points 5-8) to its first half (the
%   centroid of points 1-4). The tail direction is measured as the angular
%   direction from the its second half (the centroid of points 42-45) to
%   its first half (the centroid of points 46-49).
%

SI = seg_worm.skeleton_indices;

%For each set of indices, compute the centroids of the tip and tail then
%compute a direction vector between them (tip - tail)
TIP_INDICES  = {SI.HEAD_INDICES     SI.HEAD_TIP_INDICES     SI.TAIL_TIP_INDICES};
TAIL_INDICES = {SI.TAIL_INDICES     SI.HEAD_BASE_INDICES    SI.TAIL_BASE_INDICES};

%These are the names of the final fields
NAMES = {'tail2head' 'head' 'tail'};

directions = struct;
for iVector = 1:3
    
   tip_x   = mean(sx(TIP_INDICES{iVector},:),1);
   tip_y   = mean(sy(TIP_INDICES{iVector},:),1);
   tail_x  = mean(sx(TAIL_INDICES{iVector},:),1);   
   tail_y  = mean(sx(TAIL_INDICES{iVector},:),1);  
   
   directions.(NAMES{iVector}) = 180/pi*atan2(tip_y - tail_y, tip_x - tail_x);
end

end