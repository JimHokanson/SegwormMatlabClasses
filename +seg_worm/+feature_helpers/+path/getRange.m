function range = getRange(nw)
%
%
%   range = seg_worm.feature_helpers.path.getRange(nw)
%
% The range is defined, per frame, as the distance of the worm’s midbody
% from its final path centroid. The central dot displays the final path
% centroid. The two arrows display the range at early and late times within
% the experiment. (i) The locations of worm dwelling are shown as a
% heatmap. A single location of dwelling dominates faint traces of the
% worm’s path during motion.

%NOTE: This code uses the contours, but the skeleton is probably just as
%good and is only slightly different ...
contour_x = nw.contour_x;
contour_y = nw.contour_y;

%Average first by frame ...
%---------------------------------------
mean_cx = mean(contour_x,1);
mean_cy = mean(contour_y,1);

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