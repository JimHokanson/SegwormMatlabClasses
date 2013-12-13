function posture = getPostureFeatures(nw,midbody_distance,FPS)
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

N_ECCENTRICITY = 50; %grid size for estimating eccentricity, this is the
%max # of points that will fill the wide dimension

%Bends - DONE
%---------------------------------------------------------------------
posture.bends = seg_worm.feature_helpers.posture.getPostureBends(nw.angles);


%Eccentricity & Orientation - DONE
%--------------------------------------------------------------------------
%??? - Would we be able to identify simple worms by integrating the angles?
%cumsum(angles) - ensure cumsum never exceeds some value
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


%Coils - NOT YET FINISHED - just added distance, needs documentation
%--------------------------------------------------------------------------
posture.coils = seg_worm.feature_helpers.posture.getCoils(nw.frame_codes,midbody_distance,FPS);


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



