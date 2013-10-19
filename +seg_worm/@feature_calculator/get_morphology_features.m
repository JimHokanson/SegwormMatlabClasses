%get_morphology_features

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


[wormAreaBlock, wormLenBlock, wormWidthBlock, wormThicknessBlock,wormFatnessBlock] = ...
    morphology_process(hObject, eventdata, handles, fileInfo, mainBlock);


%worm.morphology.length - length of the skeleton (when segemented),
%otherwise NaN
%-> call to morphology_process

%Skeleton lengths - when segmented, otherwise NaN
lengths   = featureData.wormLength;
worm.morhphology.length = lengths;

%worm.morphology.length -> lengths