function morphology = get_morphology_features(nw)


keyboard

%TODO: Do I need to replace some with NaN where not segmented
%or has this already been done ?????
%TODO: Update normalized_worm to indicate the result of the question above
morphology.length = nw.lengths;





%INPUTS SO FAR
%------------------------------------------
%Normalized worm skeleton - 49 points ...

total_skeleton_length = []; %NOTE: Maybe this should be obtained from the
%skeleton class instead of the chain code lengths?
widths_along_skeleton = []; %How do we get this from the contour?


% worm.morphology.length
% worm.morphology.width.head
% worm.morphology.width.midbody
% worm.morphology.width.tail
% worm.morphology.area
% worm.morphology.areaPerLength
% worm.morphology.widthPerLength


%     % Morphology feature set
%     [wormAreaBlock, wormLenBlock, wormWidthBlock, wormThicknessBlock,...
%         wormFatnessBlock] = morphology_process(hObject, eventdata, handles, fileInfo, mainBlock);
%     
%     featureData.area        = [featureData.area,        wormAreaBlock];
%     featureData.wormLength  = [featureData.wormLength,  wormLenBlock];
%     featureData.width       = [featureData.width,       wormWidthBlock];
%     featureData.thickness   = [featureData.thickness,   wormThicknessBlock];
%     featureData.fatness     = [featureData.fatness,     wormFatnessBlock];
% 
% % The worm widths.
% widths = struct( ...
%     'head',     headWidths,     ...
%     'midbody',  midbodyWidths,  ...
%     'tail',     tailWidths);
% 
% %% The worm morphology.
% morphology = struct( ...
%     'length',           lengths,    ...
%     'width',            widths,     ....
%     'area',             areas,      ...
%     'areaPerLength',    fatness,    ...
%     'widthPerLength',   thickness);

%https://github.com/JimHokanson/mrc_wormtracker_gui/blob/master/Features/morphology_process.m

[wormAreaBlock, wormLenBlock, wormWidthBlock, wormThicknessBlock,wormFatnessBlock] = ...
    morphology_process(hObject, eventdata, handles, fileInfo, mainBlock);

%wormWidthBlock - width at midpoint of body
%


%worm.morphology.length - length of the skeleton (when segemented),
%otherwise NaN
%-> call to morphology_process

%Skeleton lengths - when segmented, otherwise NaN
lengths   = featureData.wormLength;
worm.morhphology.length = lengths;


%Relies on: norm_worm.m
%From morphology_process.m
%--------------------------------------------------------
%worm.morphology.length -> lengths

%width - mean(widths(17:33)) %assuming 49 points

%area - sum of 9,10,11,12 - norm_worm
%thickness = width/length
%fatness   = area/length
%-------------------------------------------------------

%From: schaferFeatures_process.m

%    [widthsArrayBlock, eccentricityArrayBlock,...
        trackLengthBlock, amplitudeRatioBlock, ...
        meanOfBendAnglesBlock, stdOfBendAnglesBlock, croninAmplitudeBlock,...
        croninWavelength1Block, croninWavelength2Block, numKinksBlock,...
        skeletonAnglesBlock, skeletonMeanAnglesBlock, projectedAmpsBlock]...
        = schaferFeatures_process(hObject, eventdata, handles, fileInfo, mainBlock, eigenWorms);

headWidths    = featureData.widthsAtTips(1,:);
midbodyWidths = featureData.width;
tailWidths    = featureData.widthsAtTips(2,:);
centroidPathX = featureData.outlineCentroid(1,:);
centroidPathY = featureData.outlineCentroid(2,:);
wormSamples   = NUMBER_OF_POINTS;



% X worm.morphology.length
% worm.morphology.width.head
% worm.morphology.width.midbody
% worm.morphology.width.tail
% X worm.morphology.area
% X worm.morphology.areaPerLength
% X worm.morphology.widthPerLength
