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
coilFrames = wormTouchFrames(frameCodes, fps);

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



