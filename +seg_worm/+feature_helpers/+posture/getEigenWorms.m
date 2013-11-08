function projected_amps = getEigenWorms(nw)
%
%   projected_amps = seg_worm.feature_helpers.posture.getEigenWorms(nw);
%
%   OUTPUTS
%   =====================================================
%   projected_amps : [6 n]

NUMBER_EIGENWORMS_USE = 6;
INTERP_NAN = false; %Whether to interpolate over missing frames ...

eigenAngles = helper__getAngleArray(nw.x,nw.y);

[projected_amps, ~] = helper__eigenWormProject(nw.eigen_worms, eigenAngles, NUMBER_EIGENWORMS_USE, INTERP_NAN);

projected_amps = projected_amps';

end


function angles = helper__getAngleArray(x,y)
%
%
%   [angles, mean_angles] = helper__getAngleArray(x,y)
%
%   INPUTS
%   ==============================================
%   x : [49 n]
%   y : [49 n]
%
%   OUTPUTS
%   ===============================================
%   angles : [n 48]
%   mean_angles : [n 48]
%
%
%NOTE ON IMPLEMENTATION
%These angles are not calculated
%in the same way that the other angles are
%They are relative to the x,y coordinate system, although they
%are mean subtracted ...

%i.e.
%
%   OLD ANGLES
%        
%    2   
%   1 3   - angle is computed with 2 as the vertex, value might be 90 here
%   ...
%
%   NEW ANGLES
%
%    2
%   1     - angle is from 1 to 2, relative to x-y coordinates, angle would
%   be 45 (although later it is mean subtracted)

angles = atan2(diff(y,1,1),diff(x,1,1));

n_frames = size(x,2);

% need to deal with cases where angle changes discontinuously from -pi
% to pi and pi to -pi.  In these cases, subtract 2pi and add 2pi
% respectively to all remaining points.  This effectively extends the
% range outside the -pi to pi range.  Everything is re-centred later
% when we subtract off the mean.
false_start = false(1,n_frames);

%NOTE: By adding the row of falses, we shift the trues
%to the next value
mask_pos = [false_start; diff(angles,1,1) > pi];
mask_neg = [false_start; diff(angles,1,1) < -pi];
%[49 x n]

%Only fix the frames we need to
fix_frames_I = find(any(mask_pos | mask_neg,1));

for iFrame = 1:length(fix_frames_I)
    cur_frame = fix_frames_I(iFrame);
   
    positive_jump_I = find(mask_pos(:,cur_frame));
    negative_jump_I = find(mask_neg(:,cur_frame));
    
    % subtract 2pi from remainging data after positive jumps
    for j = 1:length(positive_jump_I)
        angles(positive_jump_I(j):end,cur_frame) = angles(positive_jump_I(j):end,cur_frame) - 2*pi;
    end
    
    % add 2pi to remaining data after negative jumps
    for j = 1:length(negative_jump_I)
        angles(negative_jump_I(j):end,cur_frame) = angles(negative_jump_I(j):end,cur_frame) + 2*pi;
    end
end

mean_angles = mean(angles,1);

angles = bsxfun(@minus,angles,mean_angles);

%:/ to be the same as they were previously ...
angles = angles';

end

function [projectedAmps, isInterpolated] = helper__eigenWormProject(eigenWorms, angleArray, numEigWorms, interpNaN)

%EIGENWORMPROJECT Project a series of skeleton angles onto the given
%                 eigenworms
%
%   [PROJECTEDAMPS ISINTERPOLATED] = EIGENWORMPROJECT(EIGENWORMS, ANGLEARRAY, NUMEIGWORMS, INTERPNAN)
%
%   Input:
%       eigenWorms    - the basis eigenWorms that the skeleton angles will
%                       be projected onto. Must be produced from skeletons
%                       of the same length as those used for angleArray.
%       angleArray    - a numFrames by numSkelPoints - 1 array of skeleton
%                       tangent angles that have been rotated to have a mean
%                       angle of zero.
%       numEigWorms   - the number of eigenworms to use in the projection.
%                       We simply take the first numEigWorms dimensions from
%                       eigenWorms and project onto those.
%       interpNaN     - Logical. Should missing data (e.g. stage motions, coiled
%                       shapes, dropped frames) be interpolated over?
%
%   Output:
%       projectedAmps - a numEigWorms dimensional time series of length
%                       numFrames containing the projected amplitudes for
%                       each eigenworm in each frame of the video
%       isInterpolated - a vector with length equal to the number of frames
%                        indicating whether the given value was
%                        interpolated (1) or not (0).

% check for empty or all-NaN input data and a request for interpolation.
% If interpNaN is 0 then result will be all NaNs and no error is necessary.
if isempty(~isnan(angleArray)) && interpNaN
    error('angleArray contains only NaN values so there is nothing to interpolate.')
end

%eigenWorms
%[D 48]
%
%angleArray [n x 48]

projectedAmps = angleArray*eigenWorms(1:numEigWorms,:)';

%if any of the projected amplitudes is NaN at a given frame, they all
%should be, so simply set interpolated flag based on first eigenworm.
isInterpolated = isnan(projectedAmps(:,1));

if interpNaN
    
    %interpolate over NaN values.  'extrap' is necessary in case the
    %first values are NaN which happens occasionally
    projectedAmps(isInterpolated,:) = interp1(...
        repmat(find(~isInterpolated),[1 numEigWorms]),...
        projectedAmps(~isInterpolated,:),...
        repmat(find(isInterpolated),[1 numEigWorms]),...
        'linear','extrap');
end

end




%{
function [angleArray, meanAngles] = helper__getAngleArray(x,y)




x = x'; %[49 n] to [n 49]
y = y';


[numFrames, lengthX] = size(x);


% initialize arrays
angleArray = zeros(numFrames, lengthX-1);
meanAngles = zeros(numFrames, 1);

tic
% for each video frame
for i = 1:numFrames
    
    % calculate the x and y differences
    dX = diff(x(i,:));
    dY = diff(y(i,:));
    
    % calculate tangent angles.  atan2 uses angles from -pi to pi instead...
    % of atan which uses the range -pi/2 to pi/2.
    angles = atan2(dY, dX);
    
    % need to deal with cases where angle changes discontinuously from -pi
    % to pi and pi to -pi.  In these cases, subtract 2pi and add 2pi
    % respectively to all remaining points.  This effectively extends the
    % range outside the -pi to pi range.  Everything is re-centred later
    % when we subtract off the mean.
    
    % find discontinuities larger than pi (skeleton cannot change direction
    % more than pi from one segment to the next)
    positiveJumps = find(diff(angles) > pi) + 1; %+1 to cancel shift of diff
    negativeJumps = find(diff(angles) < -pi) + 1;
    
    % subtract 2pi from remainging data after positive jumps
    for j = 1:length(positiveJumps)
        angles(positiveJumps(j):end) = angles(positiveJumps(j):end) - 2*pi;
    end
    
    % add 2pi to remaining data after negative jumps
    for j = 1:length(negativeJumps)
        angles(negativeJumps(j):end) = angles(negativeJumps(j):end) + 2*pi;
    end
    
    % rotate skeleton angles so that mean orientation is zero
    meanAngle     = mean(angles(:));
    meanAngles(i) = meanAngle;
    angles        = angles - meanAngle;
    
    % append to angle array
    angleArray(i,:) = angles;
end
toc



end


function [projectedAmps, isInterpolated] = helper__eigenWormProject(eigenWorms, angleArray, numEigWorms, interpNaN)

%helper__eigenWormProject Project a series of skeleton angles onto the given eigenworms
%
%   [projectedAmps isInterpolated] = helper__eigenWormProject(EIGENWORMS, ANGLEARRAY, NUMEIGWORMS, INTERPNAN)
%
%   Input:
%       eigenWorms    - the basis eigenWorms that the skeleton angles will
%                       be projected onto. Must be produced from skeletons
%                       of the same length as those used for angleArray.
%       angleArray    - a numFrames by numSkelPoints - 1 array of skeleton
%                       tangent angles that have been rotated to have a mean
%                       angle of zero.
%       numEigWorms   - the number of eigenworms to use in the projection.
%                       We simply take the first numEigWorms dimensions from
%                       eigenWorms and project onto those.
%       interpNaN     - Logical. Should missing data (e.g. stage motions, coiled
%                       shapes, dropped frames) be interpolated over?
%
%   Output:
%       projectedAmps - a numEigWorms dimensional time series of length
%                       numFrames containing the projected amplitudes for
%                       each eigenworm in each frame of the video
%       isInterpolated - a vector with length equal to the number of frames
%                        indicating whether the given value was
%                        interpolated (1) or not (0).

% check for empty or all-NaN input data and a request for interpolation.
% If interpNaN is 0 then result will be all NaNs and no error is necessary.
if isempty(~isnan(angleArray)) && interpNaN
    error('angleArray contains only NaN values so there is nothing to interpolate.')
end

projectedAmps  = NaN(size(angleArray, 1), numEigWorms);
isInterpolated = zeros(size(angleArray, 1), 1);

%Calculate time series of projections onto eigenworms
%for each frame
for i = 1:size(angleArray, 1)
    rawAngles = angleArray(i,:);
    
    %for each eigenworm
    for j = 1:numEigWorms
        projectedAmps(i,j) = sum(eigenWorms(j,:).*rawAngles);
    end
end

if interpNaN
    %interpolate over NaNs in projectedAmps
    projectedAmpsNoNaN = projectedAmps;
    %interpolate over NaN values.  'extrap' is necessary in case the
    %first values are NaN which happens occasionally
    for j = 1:numEigWorms
        pAmp = projectedAmps(:,j);
        pAmp(isnan(pAmp)) = interp1(find(~isnan(pAmp)),...
            pAmp(~isnan(pAmp)), find(isnan(pAmp)),'linear', 'extrap');
        projectedAmpsNoNaN(:,j) = pAmp;
    end
    %if any of the projected amplitudes is NaN at a given frame, they all
    %should be, so simply set interpolated flag based on first eigenworm.
    isInterpolated = isnan(projectedAmps(:,1));
    
    %now set projectedAmps to the interpolated values with no NaNs
    projectedAmps = projectedAmpsNoNaN;
end

end

%}
