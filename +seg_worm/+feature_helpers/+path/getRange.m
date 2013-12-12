function range = getRange(contour_x,contour_y)
%
%   seg_worm.feature_helpers.path.getRange
%
%   Old Name: wormPathRange.m
%
%   Nature Methods Description
%   =======================================================================
%   The centroid of the worm’s entire path is computed. The range is
%   defined as the distance of the worm’s midbody from this overall
%   centroid, in each frame (Supplementary Fig. 4h).

%Average by frame
%---------------------------------------
mean_cx = mean(contour_x,1);
mean_cy = mean(contour_y,1);

%Averages over all frames for subtracting
%--------------------------------------------------------
x_centroid_cx = nanmean(mean_cx);
y_centroid_cy = nanmean(mean_cy);

range = sqrt((mean_cx - x_centroid_cx).^2 + (mean_cy - y_centroid_cy).^2);

%Skeleton version
%---------------------------------------------------
%{
s_x = nw.x;
s_y = nw.y;

mean_sx = mean(s_x,1);
mean_sy = mean(s_y,1);

x_centroid_sx = nanmean(mean_sx);
y_centroid_sy = nanmean(mean_sy);

range2 = sqrt((mean_sx - x_centroid_sx).^2 + (mean_sy - y_centroid_sy).^2);
%}





end