function directions = getDirections(skeleton_XY)
%
%
%   seg_worm.feature_helpers.posture.getDirections
%
%   INPUTS
%   ===========================================================
%   skeleton_x : [49 n frames]
%   skeleton_y : [49 n frames]
%
%   skeleton_XY : [49 (x y) frames]
%
%   OUTPUTS
%   ===========================================================
%   directions : (struct)
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
   %Take mean over segments (dim 1)
   tip_centroid  = squeeze(mean(skeleton_XY(TIP_INDICES{iVector},:,:),1));
   tail_centroid = squeeze(mean(skeleton_XY(TAIL_INDICES{iVector},:,:),1));
   directions.(NAMES{iVector}) = 180/pi*atan2(tip_centroid(2,:) - tail_centroid(2,:), tip_centroid(1,:) - tip_centroid(1,:));
end

end