function initialize(obj,contour)

lfCMaxP = contour.lf_ap_max;
lfCMaxI = contour.lf_ap_max_I;
lfCMinP = contour.lf_ap_min;
lfCMinI = contour.lf_ap_min_I;

headI = contour.head_I;
tailI = contour.tail_I;
contour_pixels = contour.pixels;
cc_lengths     = contour.cc_lengths;
wormSegLength  = contour.avg_segment_length(true);

% Compute the worm's skeleton.
[obj.pixels,obj.c_widths] = obj.linearSkeleton(...
    headI, tailI, lfCMinP, lfCMinI, ...
    lfCMaxP, lfCMaxI, contour_pixels, wormSegLength, cc_lengths);

end