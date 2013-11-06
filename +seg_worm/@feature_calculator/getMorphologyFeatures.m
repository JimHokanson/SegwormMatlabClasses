function morphology = getMorphologyFeatures(nw)
%
%   INPUTS
%   ===========================================================
%   nw : seg_worm.normalized_worm
%
%   TODO:
%   It is not clear that unsegmented frames have NaN values. This should
%   be verified.
%

%Old files that served as a reference ...
%------------------------------------------------------------
%morphology_process.m
%schaferFeatures_process.m


MIDBODY_INDICES = 17:33;    %round((2/6*NOOFPOINTS))+1:round((4/6*NOOFPOINTS)
HEAD_INDICES    = 1:8;      %1:round(1/6*NOOFPOINTS)
TAIL_INDICES    = 42:49;    %round(5/6*NOOFPOINTS)+1:end

%TODO: Do I need to replace some with NaN where not segmented
%or has this already been done ?????
%TODO: Update normalized_worm to indicate the result of the question above

%These look to be equivalent, should do formal check ...
% I1 = find(not_segmented_mask);
% I2 = find(isnan(morphology.length));

morphology.length = nw.lengths;

morphology.width.head     = mean(nw.widths(HEAD_INDICES,:),1);
morphology.width.midbody  = mean(nw.widths(MIDBODY_INDICES,:),1);
morphology.width.tail     = mean(nw.widths(TAIL_INDICES,:),1);
morphology.area           = nw.head_areas + nw.tail_areas + nw.vulva_areas + nw.non_vulva_areas;
morphology.areaPerLength  = morphology.area./morphology.length;
morphology.widthPerLength = morphology.width.midbody./morphology.length;
