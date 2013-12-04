function morphology = getMorphologyFeatures(nw)
%
%   morphology = seg_worm.feature_calculator.getMorphologyFeatures(nw)
%
%   INPUTS
%   =======================================================================
%   nw : seg_worm.normalized_worm
%
%
%   Nature Methods Description
%   =======================================================================
%
%   Morphology Features 
%
%   1. Length. Worm length is computed from the segmented skeleton by
%   converting the chain-code pixel length to microns.
%   
%   2. Widths. Worm width is computed from the segmented skeleton. The
%   head, midbody, and tail widths are measured as the mean of the widths
%   associated with the skeleton points covering their respective sections.
%   These widths are converted to microns.
%   
%   3. Area. The worm area is computed from the number of pixels within the
%   segmented contour. The sum of the pixels is converted to microns2.
%
%   4. Area/Length.
%
%   5. Midbody Width/Length.


%Old files that served as a reference ...
%------------------------------------------------------------
%morphology_process.m
%schaferFeatures_process.m

SI = seg_worm.skeleton_indices;

morphology.length         = nw.lengths;
morphology.width.head     = mean(nw.widths(SI.HEAD_INDICES,:),1);
morphology.width.midbody  = mean(nw.widths(SI.MID_INDICES,:),1);
morphology.width.tail     = mean(nw.widths(SI.TAIL_INDICES,:),1);
morphology.area           = nw.head_areas + nw.tail_areas + nw.vulva_areas + nw.non_vulva_areas;
morphology.areaPerLength  = morphology.area./morphology.length;
morphology.widthPerLength = morphology.width.midbody./morphology.length;
