function posture = getPostureFeatures(nw)
%
%   posture = seg_worm.feature_calculator.getPostureFeatures(nw)
%
%   Old Files
%   - schaferFeatures_process
%
%   NOTES:
%   - Indices were inconsistently defined for bends relative to other code
%   - stdDev for bends is signed as well, based on means ...
%
%
%   MISSING INPUTS - framecodes when parsing ...


%MISSING FILE
%-------------------------------
%getAmpWavelength
%
% BMC Genetics, 2005
% C.J. Cronin, J.E. Mendel, S. Mukhtar, Young-Mee Kim, R.C. Stirb, J. Bruck,
% P.W. Sternberg
% "An automated system for measuring parameters of nematode
% sinusoidal movement" BMC Genetics 2005, 6:5
%
%   Code supposedly at:
%   http://wormlab.caltech.edu/publications/CaltechTracker(BMCGenetics2005).zip
%
%   Doesn't contain getAmpWavelength

%keyboard

%TODO: Also should be dependent method of object ...
%segmented_mask = nw.segmentation_status == 's'

N_ECCENTRICITY = 50;


%{
worm.posture.amplitude.max
worm.posture.amplitude.ratio

worm.posture.wavelength.primary
worm.posture.wavelength.secondary

worm.posture.tracklength
worm.posture.eccentricity
worm.posture.kinks

worm.posture.coils

worm.posture.directions.tail2head
worm.posture.directions.head
worm.posture.directions.tail

DONE worm.posture.skeleton
DONE worm.posture.eigenProjection

%}



%Bends
%---------------------------------------------------------------------

%NOTE: Most of these were defined slighly wrong (in my opinion)
%in schaferFeatures_process.m
%I've changed them, although it might be good to push these up a level
%to the calculator class as constants ... (or even higher)


SI          = seg_worm.skeleton_indices;
ALL_INDICES = SI.ALL_NORMAL_INDICES;
FIELDS      = SI.ALL_NORMAL_NAMES;

n_fields = length(FIELDS);

bends = struct;
for iField = 1:n_fields
    cur_indices = ALL_INDICES{iField};
    cur_name    = FIELDS{iField};
    bends.(cur_name).mean = nanmean(nw.angles(cur_indices,:));
    bends.(cur_name).stdDev = nanstd(nw.angles(cur_indices,:));
    mask = bends.(cur_name).mean < 0;
    bends.(cur_name).stdDev(mask) = -1*bends.(cur_name).stdDev(mask);
end

posture.bends = bends;

%Eccentricity & Orientation
%--------------------------------------------------------------------------
%This essentially reconstructs the contour from the side contours
%This should probably just be a dependent method on the object ...
%
%   The concatentation order doesn't matter
contour_x = squeeze([nw.vulva_contours(:,1,:); nw.non_vulva_contours(end-1:-1:2,1,:);]);
contour_y = squeeze([nw.vulva_contours(:,2,:); nw.non_vulva_contours(end-1:-1:2,2,:);]);

%x,y - > 96 x n (assuming 49*2 - 2)


%THIS IS REALLY SLOW !!!!!!
%
%   Requires filling in the worm
%
%[posture.eccentricity, worm_orientation] = seg_worm.feature_helpers.posture.getEccentricity(contour_x, contour_y, N_ECCENTRICITY);
%??? How would pca compare to the algorithm used???




%Amplitude, Wavelengths, TrackLength, Amplitude Ratio
%--------------------------------------------------------------------------
%
%   STATUS: This requires a function which I am missing
%

%From Sternberg code, might use:
%tracks.m           - ampt, wavelnth
%
%   ??? What is trackLen and ampRatio

%   -> ampt - maximum amplitude along the worm body
%   -> ampRatio - maximum on opposiding sides, smaller is numerator
%   -> uses orientation from ellipse
%   -> trackLen - horizontal length of projection
%       -> I think this is the opposite of ampt
%
%       trackLen/ampt -> ampRatio???

%wavelength - x must be monotonic after rotation
%primary - peak of fft
%secondary - 2nd peak of fft, must exceed 0.5 * 1st peak
%
%   NOTE: Wavelength can never exceed 2x worm length
%   -> If greater then what? NaN?
%
%
%   ???? - velocity based, but their code isn't ...

%{

%MISSING FUNCTION
[ampt, wavelnths, trackLen, ampRatio] = getAmpWavelength(thetaValRad, skCoords, wormLen, mydata.mainImg, guiItem, 0, 0);

worm.posture.amplitude.max        = ampt;
worm.posture.amplitude.ratio      = ampRatio;
worm.posture.wavelength.primary   = wavelnths(1);
worm.posture.wavelength.secondary = wavelnths(2);
worm.posture.tracklength = trackLen;

%}

%Kinks - CODE NOT YET EXAMINED ... (But it works)
%--------------------------------------------------------------------------
posture.kinks = seg_worm.feature_helpers.posture.wormKinks(nw.angles);

%Coils - NOT YET FINISHED ...
%--------------------------------------------------------------------------

%{

This is really complicated because it relies on parsing errors to detect
coiling events:

wormTouchFrames - main function


%}




%   worm.posture.coils
%
%        frames: [1x5 struct]
%     frequency: 0.0278
%     timeRatio: 0.0291
%
%   frames(1)
%             start: 267
%               end: 300
%              time: 1.3158
%         interTime: 22.4073
%     interDistance: 5.1523e+03


%INPUTS
%==========================
%-> parsing error codes ...
%-> distance
%-> fps


%{

% Compute the coiled frames.
%
%   Is this run on a frame by frame basis????


frameLabels = [frameLabels, data{1}]; %#ok<AGROW>

% This section will compute warning frame labels
failedFrames = [];
load(fileInfo.expList.failedFramesFile, 'failedFrames');

% Older version <3 of segmentation had a bug in saving the failed frames.
% They have been indexed starting 0 not 1 (because of frame number being
% generated from time stamp rather than globalFrameCounter). To counter act
% it the indices for the frames that failed need to be added 1 to shift the
% failed frames by one and re-allign them. Here we will make a check for
% that and will raise a flag to add 1 in the upcoming loop
shiftFailedFrames = 0;
if ~isempty(failedFrames) && length(failedFrames(:,1)) > 2
    if sum(frameLabels(failedFrames(2:end,1))~='f') ~= 0
        shiftFailedFrames = 1;
    end
end

for i=1:length(failedFrames(:,1))
    % Here for each failed frame we add a warning id, this warning ID will
    % always strat from a value 101 and end 1000. Everything greater than
    % 100 is a failed frame. If 100 is subtracted then the number
    % corresponds to the array element in handles.warningList
    
    if shiftFailedFrames
        failedFrameIndex = failedFrames(i,1) + 1;
    else
        failedFrameIndex = failedFrames(i,1);
    end
    
    if failedFrameIndex > 0 && failedFrameIndex < length(numFrameLabel)
        numFrameLabel(failedFrameIndex) = failedFrames(i,2) + 100;
    else
        warningStr = strcat('Failed frame index is:', {' '}, num2str(failedFrameIndex), {' '},'and the length of the video is:', {' '}, num2str(numFrameLabel),'.');
        warning('features:failedFrameLabel',warningStr{1});
    end
end

% We have three kinds of labels each of them will have a number value assigned.
% 1. from 1-10 general labels
% 2. from 11-100 normWorm labels
% 3. from 101-1000 segmentation labels

numFrameLabel = zeros(1, length(frameLabels));

numFrameLabel(frameLabels == 's') = 1;
numFrameLabel(frameLabels == 'm') = 2;
numFrameLabel(frameLabels == 'd') = 3;
numFrameLabel(frameLabels == 'n') = 1001;


numFrameLabel

coilFrames = seg_worm.feature_helpers.posture.wormTouchFrames(frameCodes, fps);


%coilFrames (struct array)
%.start - start frame index

% Compute the coiled statistics.
[coilEventStats, coiledStats] = events2stats(coilFrames, fps, distance, [], 'interDistance');

% Reorganize everything for the feature file.
coilFrames = coilEventStats;
coilFrequency = [];
coilTimeRatio = [];
if ~isempty(coiledStats)
    coilFrequency = coiledStats.frequency;
    coilTimeRatio = coiledStats.ratio.time;
end
coils = struct( ...
    'frames',       coilFrames, ...
    'frequency',    coilFrequency, ...
    'timeRatio',    coilTimeRatio);

%}




%Directions
%--------------------------------------------------------------------------
posture.directions = seg_worm.feature_helpers.posture.getDirections(nw.skeletons);

%Skeleton
%--------------------------------------------------------------------------
posture.skeleton.x = nw.x;
posture.skeleton.y = nw.y;

%EigenProjection
%--------------------------------------------------------------------------
posture.eigenProjection = seg_worm.feature_helpers.posture.getEigenWorms(nw);



