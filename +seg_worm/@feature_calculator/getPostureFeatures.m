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
%   UNFINISHED STATUS:
%   - seg_worm.feature_helpers.posture.wormKinks - not yet examined
%   - distance - missing input to function, need to process locomotion
%   first
%



FPS = 20; %TODO: get this from higher up ...
N_ECCENTRICITY = 50;

%Bends
%---------------------------------------------------------------------
SI          = seg_worm.skeleton_indices;
ALL_INDICES = SI.ALL_NORMAL_INDICES;
FIELDS      = SI.ALL_NORMAL_NAMES;

n_fields = length(FIELDS);

bends = struct;
for iField = 1:n_fields
    cur_indices = ALL_INDICES{iField};
    cur_name    = FIELDS{iField};
    bends.(cur_name).mean   = nanmean(nw.angles(cur_indices,:));
    bends.(cur_name).stdDev = nanstd(nw.angles(cur_indices,:));
    
    %Sign the standard deviation ...
    %----------------------------------------------------------------------
    mask = bends.(cur_name).mean < 0;
    bends.(cur_name).stdDev(mask) = -1*bends.(cur_name).stdDev(mask);
end

posture.bends = bends;



%Eccentricity & Orientation
%--------------------------------------------------------------------------
%This is relatively slow ...
[posture.eccentricity, worm_orientation] = seg_worm.feature_helpers.posture.getEccentricity(nw.contour_x, nw.contour_y, N_ECCENTRICITY);



%Amplitude, Wavelengths, TrackLength, Amplitude Ratio
%--------------------------------------------------------------------------
[posture.amplitude,posture.wavelength,posture.trackLength] = ...
  seg_worm.feature_helpers.posture.getAmplitudeAndWavelength(worm_orientation,nw.x,nw.y,nw.lengths);



%Kinks - CODE NOT YET EXAMINED ... (But it works)
%--------------------------------------------------------------------------
posture.kinks = seg_worm.feature_helpers.posture.wormKinks(nw.angles);

%Coils - NOT YET FINISHED ...
%--------------------------------------------------------------------------

%NOTE: we need distance!, moving onto locomotion to get it ...

distance = rand(1,nw.n_frames); %This is temporary ...

%NOTE: Output is zero based currently :/
coiled_frames = seg_worm.feature_helpers.posture.wormTouchFrames(nw.frame_codes, FPS);

%TODO: I might want to create a static method of the class so I can merge
%these two lines ...
coiled_events = seg_worm.feature.event(coiled_frames,FPS,distance,[],'interDistance');
posture.coils = coiled_events.getStruct;


%Directions
%--------------------------------------------------------------------------
%seg_worm.feature_helpers.posture.getDirections
posture.directions = seg_worm.feature_helpers.posture.getDirections(nw.skeletons);



%Skeleton
%--------------------------------------------------------------------------
posture.skeleton.x = nw.x;
posture.skeleton.y = nw.y;



%EigenProjection
%--------------------------------------------------------------------------
posture.eigenProjection = seg_worm.feature_helpers.posture.getEigenWorms(nw);



