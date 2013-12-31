function bends = getPostureBends(bend_angles)
%
%
%   Output
%   =======================================================================
%   bends
%       .head
%           .mean
%           .stdDev
%       .neck
%           .mean
%           .stdDev
%       .midbody
%           .mean
%           .stdDev
%       .hips
%           .mean
%           .stdDev
%       .tail
%           .mean
%           .stdDev
%
%
%   Old Names: featureProcess.m, schaferFeatures_process.m
%
%   Nature Methods Description
%   =======================================================================
%   Bends. 
%   -------------------
%   Worm bending is measured using the supplementary angles to the bends
%   formed along the skeleton, with each skeleton point serving as the
%   vertex to its respective bend (Supplementary Fig. 4b). The
%   supplementary angle can also be expressed as the difference in tangent
%   angles at the skeleton point. The supplementary angle provides an
%   intuitive measurement. Straight, unbent worms have an angle of 0°.
%   Right angles are 90°. And the largest angle theoretically possible, a
%   worm bending back on itself, would measure 180°. The supplementary
%   angle is determined, per skeleton point, using edges 1/12 the
%   skeleton’s chaincode length, in opposing directions, along the
%   skeleton. When insufficient skeleton points are present, the angle
%   remains undefined (i.e., the first and last 1/12 of the skeleton have
%   no bending angle defined). The mean and standard deviation are measured
%   for each body segment. The angle is signed to provide the bend’s
%   dorsal-ventral orientation. When the worm has its ventral side internal
%   to the bend, the bending angle is signed negatively.

SI = seg_worm.skeleton_indices;

ALL_NORMAL_INDICES = SI.ALL_NORMAL_INDICES;
FIELDS             = SI.ALL_NORMAL_NAMES;

n_fields = length(FIELDS);

bends = struct;
for iField = 1:n_fields
    cur_indices = ALL_NORMAL_INDICES{iField};
    cur_name    = FIELDS{iField};
    bends.(cur_name).mean   = nanmean(bend_angles(cur_indices,:));
    bends.(cur_name).stdDev = nanstd(bend_angles(cur_indices,:));
    
    %Sign the standard deviation ...
    %----------------------------------------------------------------------
    mask = bends.(cur_name).mean < 0;
    bends.(cur_name).stdDev(mask) = -1*bends.(cur_name).stdDev(mask);
end

end