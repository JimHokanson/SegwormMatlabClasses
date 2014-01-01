function getPostureBends(obj, bend_angles, d_opts)
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
%   Old Names: 
%   - featureProcess.m
%   - schaferFeatures_process.m
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
%
%   Jim's Summary
%   -----------------------------------------------------------------------
%   Given the "bend angles" for a given frame, calculate the mean and
%   standrard deviation of a section of bend angles (e.g. average all bend
%   angles at the head)

SI = seg_worm.skeleton_indices;
ALL_NORMAL_INDICES = SI.ALL_NORMAL_INDICES;
FIELDS             = SI.ALL_NORMAL_NAMES;



if d_opts.mimic_old_behavior
    %indices_for_mean = {1:9 9:17 17:32 31:39 39:48};  %From code
    indices_for_mean = {1:9 9:17 17:33 32:40 40:49};   %From results
    indices_for_std  = indices_for_mean;
    indices_for_std{3} = 17:32;
else
    indices_for_mean = ALL_NORMAL_INDICES;
    indices_for_std  = ALL_NORMAL_INDICES;
end

n_fields = length(FIELDS);

bends = struct;
for iField = 1:n_fields
    cur_mean_indices = indices_for_mean{iField};
    cur_std_indices  = indices_for_std{iField};
    cur_name    = FIELDS{iField};
    bends.(cur_name).mean   = nanmean(bend_angles(cur_mean_indices,:));
    bends.(cur_name).stdDev = nanstd(bend_angles(cur_std_indices,:));
    
    %Sign the standard deviation ...
    %----------------------------------------------------------------------
    mask = bends.(cur_name).mean < 0;
    bends.(cur_name).stdDev(mask) = -1*bends.(cur_name).stdDev(mask);
end

obj.bends = bends;

end

%{

Indices Mismatch
%            OLD                        NEW
%---------------------------------------------
%head     : 1:9                         1:8
%neck     : 9:17                        9:16
%midbody  : 17:32 (mean) 17:31 (std)    17:33
%hip      : 31:39                       34:41
%tail     : 39:48                       42:49

???? - had to change hip to 32:40 to get it to match ...
Perhaps the code has changed since it was last run ...

%??? - tail at 40:48 (49 would work as well since ends are NaN)


%}