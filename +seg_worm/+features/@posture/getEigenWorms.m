function getEigenWorms(obj,sx,sy,eigen_worms,N_EIGENWORMS_USE)
%
%
%   seg_worm.features.posture.getEigenWorms
%
%   Inputs
%   =======================================================================
%   sx  : [49 x n_frames]
%   sy  : [49 x n_frames]
%   eigen_worms : [7 x 48]
%
%   Outputs
%   =======================================================================
%   projected_amps : [6 n_frames]
%
%
%   Old Names: 
%   - schaferFeatures_process.m  
%   - eigenWormProject.m
%
%   Nature Methods Description
%   =======================================================================
%   Eigen Projections. 
%   ------------------------------------------------
%   The eigenworm amplitudes are a measure of worm posture.
%   They are the projections onto the first six eigenworms which together account for
%   97% of the variance in posture. The eigenworms were computed from 15 N2
%   videos (roughly 3 hours of video, 1/3 of a million frames) as previously
%   described8.
%
%   Briefly, 48 tangent angles are calculated along the skeleton and rotated to have a
%   mean angle of zero. Principal components analysis is performed on the pooled
%   angle data and we keep the 6 principal components (or eigenworms) that capture
%   the most variance. The first eigenworm roughly corresponds to body curvature.
%   The next two eigenworms are akin to sine and cosine waves encoding the
%   travelling wave during crawling. The fourth eigenworm captures most of the
%   remaining variance at the head and tail. Projected amplitudes are calculated from
%   the posture in each frame. Even for the mutants, the data is always projected onto
%   the N2-derived eigenworms.



angles = h__getAngleArray(sx,sy);

obj.eigenProjection = (eigen_worms(1:N_EIGENWORMS_USE,:)*angles);

end


function angles = h__getAngleArray(sx,sy)
%
%
%   INPUTS
%   ==============================================
%   sx : [49 n]
%   sy : [49 n]
%
%   OUTPUTS
%   ===============================================
%   angles : [n 48]
%
%NOTE ON IMPLEMENTATION
%These angles are not calculated in the same way that the other angles
%(bend angles) are They are relative to the x,y coordinate system and mean
%subtracted

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

angles = atan2(diff(sy,1,1),diff(sx,1,1));

n_frames = size(sx,2);

% need to deal with cases where angle changes discontinuously from -pi
% to pi and pi to -pi.  In these cases, subtract 2pi and add 2pi
% respectively to all remaining points.  This effectively extends the
% range outside the -pi to pi range.  Everything is re-centred later
% when we subtract off the mean.
false_row = false(1,n_frames);

%NOTE: By adding the row of falses, we shift the trues
%to the next value, which allows indices to match. Otherwise after every
%find statement we would need to add 1, I think this is a bit faster ...
mask_pos = [false_row; diff(angles,1,1) > pi];
mask_neg = [false_row; diff(angles,1,1) < -pi];
%[49 x n]

%Only fix the frames we need to, in which there is a jump in going from one
%segment to the next ...
fix_frames_I = find(any(mask_pos | mask_neg,1));

for iFrame = 1:length(fix_frames_I)
    cur_frame = fix_frames_I(iFrame);
   
    positive_jump_I = find(mask_pos(:,cur_frame));
    negative_jump_I = find(mask_neg(:,cur_frame));
    
    % subtract 2pi from remainging data after positive jumps
    % Note that the jumps impact all subsequent frames
    for j = 1:length(positive_jump_I)
        angles(positive_jump_I(j):end,cur_frame) = angles(positive_jump_I(j):end,cur_frame) - 2*pi;
    end
    
    % add 2pi to remaining data after negative jumps
    for j = 1:length(negative_jump_I)
        angles(negative_jump_I(j):end,cur_frame) = angles(negative_jump_I(j):end,cur_frame) + 2*pi;
    end
end

angles = bsxfun(@minus,angles,mean(angles,1));

%:/ to be the same as they were previously ...
%angles = angles';

end