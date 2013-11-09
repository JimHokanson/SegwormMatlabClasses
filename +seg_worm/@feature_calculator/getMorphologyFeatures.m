function morphology = getMorphologyFeatures(nw)
%
%   morphology = seg_worm.feature_calculator.getMorphologyFeatures(nw)
%
%   INPUTS
%   =======================================================================
%   nw : seg_worm.normalized_worm

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
