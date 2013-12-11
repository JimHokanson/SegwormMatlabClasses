function posture = getPostureFeatures(nw)
%
%   posture = seg_worm.feature_calculator.getPostureFeatures(nw)
%
%   Inputs
%   =======================================================================
%   nw : seg_worm.normalized_worm
%
%
%   UNFINISHED STATUS:
%   - seg_worm.feature_helpers.posture.wormKinks - not yet examined
%   - distance - missing input to function, need to process locomotion
%   first
%


FPS = 20; %TODO: get these from higher up ...
N_ECCENTRICITY = 50; %grid size for estimating eccentricity, this is the
%max # of points that will fill the wide dimension

%Bends - DONE
%---------------------------------------------------------------------
posture.bends = seg_worm.feature_helpers.posture.getPostureBends(nw.angles);


%Eccentricity & Orientation - DONE
%--------------------------------------------------------------------------
[posture.eccentricity, worm_orientation] = ...
    seg_worm.feature_helpers.posture.getEccentricity(...
    nw.contour_x, nw.contour_y, N_ECCENTRICITY);


%Amplitude, Wavelengths, TrackLength - DONE
%--------------------------------------------------------------------------
[posture.amplitude,posture.wavelength,posture.trackLength] = ...
  seg_worm.feature_helpers.posture.getAmplitudeAndWavelength(...
  worm_orientation,nw.x,nw.y,nw.lengths);


%Kinks (aka bend counts) - CODE NOT YET EXAMINED ... (But it works)
%--------------------------------------------------------------------------
posture.kinks = seg_worm.feature_helpers.posture.wormKinks(nw.angles);


%Coils - NOT YET FINISHED ... - needs distance from locomotion ...
%--------------------------------------------------------------------------
posture.coils = seg_worm.feature_helpers.posture.getCoils(nw.n_frames,nw.frame_codes);


%Directions (AKA orientation) - DONE
%--------------------------------------------------------------------------
%seg_worm.feature_helpers.posture.getDirections
posture.directions = seg_worm.feature_helpers.posture.getDirections(nw.x,nw.y);



%Skeleton - DONE
%--------------------------------------------------------------------------
posture.skeleton.x = nw.x;
posture.skeleton.y = nw.y;


%EigenProjection - DONE
%--------------------------------------------------------------------------
posture.eigenProjection = seg_worm.feature_helpers.posture.getEigenWorms(nw.x,nw.y,nw.eigen_worms);



