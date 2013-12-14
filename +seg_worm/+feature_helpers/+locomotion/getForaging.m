function foraging = getForaging(sx,sy,is_segmented_mask,ventral_mode,FPS)
%
%
%   seg_worm.feature_helpers.locomotion.getForaging
%
%
%   Old Name: 
%   - part of wormBends.m
%
%
%   Inputs
%   =======================================================================
%   sx  :
%   sy  :
%   is_segmented_mask : [1 x n_frames]
%   ventral_mode      : (scalar)
%   FPS :
%
%   Nature Methods Description
%   =======================================================================
%   Foraging
%   -----------------------
%   Foraging. Worm foraging is expressed as both an amplitude and an
%   angular speed (Supplementary Fig. 4g). Foraging is signed negatively
%   whenever it is oriented towards the ventral side. In other words, if
%   the nose is bent ventrally, the amplitude is signed negatively.
%   Similarly, if the nose is moving ventrally, the angular speed is signed
%   negatively. As a result, the amplitude and angular speed share the same
%   sign roughly only half the time. Foraging is an ambiguous term in
%   previous literature, encompassing both fine movements of the nose as
%   well as larger swings associated with the head. Empirically we have
%   observed that the nose movements are aperiodic while the head swings
%   have periodicity. Therefore, we measure the aperiodic nose movements
%   and term these foraging whereas the head swings are referred to as
%   measures of head crawling (described earlier in this section).
% 
%   Foraging movements can exceed 6Hz7 and, at 20-30fps, our video frame
%   rates are just high enough to resolve the fastest movements. By
%   contrast, the slowest foraging movements are simply a continuation of
%   the crawling wave and present similar bounds on their dynamics.
%   Therefore, we bound foraging between 1/30Hz (the lower bound used for
%   crawling) and 10Hz.
% 
%   To measure foraging, we split the head in two (skeleton points 1-4 and
%   5-8) and measure the angle between these sections. To do so, we measure
%   the mean of the angle between subsequent skeleton points along each
%   section, in the tail-to-head direction. The foraging angle is the
%   difference between the mean of the angles of both sections. In other
%   words, the foraging angle is simply the bend at the head.
% 
%   Missing frames are linearly interpolated, per each skeleton point, for
%   fragments up to 0.2 seconds long (4-6 frames at 20-30fps – twice the
%   upper foraging bound). When larger fragments are missing, foraging is
%   marked undefined. Segmentation of the head at very small time scales
%   can be noisy. Therefore, we smooth the foraging angles by convolving
%   with a Gaussian filter 1/5 of a second long (for similar reasons to
%   those mentioned in frame interpolation), with a width defined by the
%   Matlab “gausswin” function’s default ? of 2.5 and normalized such that
%   the filter integrates to 1.
% 
%   The foraging amplitude is defined as the largest foraging angle
%   measured, prior to crossing 0°. In other words, the largest nose bend
%   prior to returning to a straight, unbent position. Therefore, the
%   foraging amplitude time series follows a discrete, stair-step pattern.
%   The amplitude is signed negatively whenever the nose points towards the
%   worm’s ventral side. The foraging angular speed is measured as the
%   foraging angle difference between subsequent frames divided by the time
%   between these frames. To center the foraging angular speed at the frame
%   of interest and eliminate noise, each frame is assigned the mean of the
%   angular speed computed between the previous frame and itself and
%   between itself and the next frame. The angular speed is signed
%   negatively whenever its vector points towards the worm’s ventral side.


% Initialize the function state.
%
% Note: empirically I've found the values below achieve good signal.
%
% Furthermore ...
% Huang et al. in 2006, measure foraging frequencies for several worms and
% find the signal centered at roughly 4Hz. For N2 worms, they see a second
% signal at 10Hz but I find this value too close to the background noise
% present in segmentation. Visually inspecting the foraging signal, as the
% bend between the nose and neck, corroborates a roughly 4Hz signal. But,
% foraging usually encompasses only half to a quarter cycle. In other
% words, the worm bends it nose sharply and sometimes bends it back but a
% full wave, akin to a body bend, occurs far less frequently. Therefore I
% chose to measure angular speed for foraging.

SI = seg_worm.skeleton_indices;
NOSE_I = SI.HEAD_TIP_INDICES(end:-1:1); %flip to maintain orientation for
%angles and consistency with old code ...
NECK_I = SI.HEAD_BASE_INDICES(end:-1:1);

MIN_NOSE_WINDOW = round(0.1*FPS);
MAX_NOSE_INTERP = 2*MIN_NOSE_WINDOW - 1;

nose_x = sx(NOSE_I,:);
nose_y = sy(NOSE_I,:);
neck_x = sx(NECK_I,:);
neck_y = sy(NECK_I,:);

%Step 1: Interpolation of skeleton indices
%--------------------------------------------------------------------------
interp_nose_mask = h__getNoseInterpolationIndices(is_segmented_mask,MAX_NOSE_INTERP);

dataI       = find(is_segmented_mask);
noseInterpI = find(interp_nose_mask);

nose_xi = h__interpData(nose_x,dataI,noseInterpI);
nose_yi = h__interpData(nose_y,dataI,noseInterpI);
neck_xi = h__interpData(neck_x,dataI,noseInterpI);
neck_yi = h__interpData(neck_y,dataI,noseInterpI);

%Step 2: Calculation of the bend angles
%--------------------------------------------------------------------------
nose_bends = h__computeNoseBends(nose_xi,nose_yi,neck_xi,neck_yi);

%Step 3: 
[nose_amps,nose_freqs] = h__foragingData(nose_bends, MIN_NOSE_WINDOW, FPS);

if ventral_mode > 1
    nose_amps  = -nose_amps;
    nose_freqs = -nose_freqs;
end

foraging.amplitude  = nose_amps;
foraging.angleSpeed = nose_freqs;


end

function nose_bends_d = h__computeNoseBends(nose_x,nose_y,neck_x,neck_y)
%
%   Compute the difference in angles between the nose and neck (really the
%   head tip and head base).
%   
%
%   Inputs
%   ======================================================
%   nose_x: [4 x n_frames]
%   nose_y: [4 x n_frames]
%   neck_x: [4 x n_frames]
%   neck_y: [4 x n_frames]
%
%   Outputs
%   ======================================================
%   nose_bends_d

noseAngles = h__computeAvgAngles(nose_x,nose_y);
neckAngles = h__computeAvgAngles(neck_x,neck_y);

%TODO: These three should be a method, calculating the difference
%in angles and ensuring all results are within +/- 180
nose_bends_d  = (noseAngles - neckAngles)'*180/pi;

nose_bends_d(nose_bends_d > 180)  = nose_bends_d(nose_bends_d > 180) - 360;
nose_bends_d(nose_bends_d < -180) = nose_bends_d(nose_bends_d < -180) + 360;
end

function angles = h__computeAvgAngles(x,y)
%
%   Take average difference between successive x and y skeleton points, the
%   compute the arc tangent from those averages.
%
%   Simple helper for h__computeNoseBends

    avg_diff_x = mean(diff(x,1,1));
    avg_diff_y = mean(diff(y,1,1));
    
    angles = atan2(avg_diff_y,avg_diff_x);
end


function x_new = h__interpData(x_old,good_data_I,fix_data_I)
%
%
%   TODO: move to: seg_worm.feature_helpers.interpolateNanData
%
%   Inputs
%   =======================================================================
%   x_old       : [4 x n_frames]
%   good_data_I : [1 x m]
%   fix_data_I  : [1 x n]

x_new = x_old;

%NOTE: This version is a bit weird because the size of y is not 1d
for i1 = 1:size(x_old,1)
   x_new(i1,fix_data_I) = interp1(good_data_I,x_old(i1,good_data_I),fix_data_I,'linear',NaN); 
end

end


function interp_data_mask = h__getNoseInterpolationIndices(is_segmented_mask,max_nose_interp_sample_width)
%
%   
%   Interpolate data:
%   - but don't extrapolate
%   - don't interpolate if the gap is too large
%
%   Inputs
%   =======================================================================
%   is_segmented_mask :
%   max_nose_interp_sample_width : Maximum # of frames that 
%
%   Outputs
%   =======================================================================
%   interp_data_mask : [1 x n_frames} whether or not to interpolate a data point

%JAH NOTE: I'm considering replacing this and other interpolation code
%with a specific function for all feature processing

% Find the start and end indices for missing data chunks.
isNotData        = ~is_segmented_mask;
interp_data_mask = isNotData;
diffIsNotData    = diff(isNotData);

startNotDataI    = find(diffIsNotData == 1);
endNotDataI      = find(diffIsNotData == -1) - 1;

% Don't interpolate missing data at the very start and end.
%--------------------------------------------------------------------------
if ~isempty(startNotDataI) && (isempty(endNotDataI) || startNotDataI(end) > endNotDataI(end))
    interp_data_mask(startNotDataI(end):end) = false;
    startNotDataI(end) = [];
end

if ~isempty(endNotDataI) && (isempty(startNotDataI) || startNotDataI(1) > endNotDataI(1))
    interp_data_mask(1:endNotDataI(1)) = false;
    endNotDataI(1) = [];
end

% Don't interpolate large missing chunks of data.
for i = 1:length(startNotDataI)
    if endNotDataI(i) - startNotDataI(i) > max_nose_interp_sample_width
        interp_data_mask(startNotDataI(i):endNotDataI(i)) = false;
    end
end



end



%% Compute the foraging amplitude and angular speed.
function [amps,speeds] = h__foragingData(nose_bend_angle_d, min_win_size, fps)
%
%
%
%   Inputs
%   =======================================================================
%   nose_bend_angle_d : [n_frames x 1]
%   min_win_size : (scalar)
%   fps : (scalar)
%
%   Outputs
%   =======================================================================
%   amps   : [1 x n_frames]
%   speeds : [1 x n_frames]
%
%

% Clean up the signal with a gaussian filter.
%--------------------------------------------------------------------------
if min_win_size > 0
    gaussFilter       = gausswin(2 * min_win_size + 1) / min_win_size;
    nose_bend_angle_d = conv(nose_bend_angle_d, gaussFilter, 'same');
    
    %Remove partial data frames ...
    nose_bend_angle_d(1:min_win_size) = NaN;
    nose_bend_angle_d((end - min_win_size + 1):end) = NaN;
end

%Calculate amplitudes
%--------------------------------------------------------------------------
amps = h__getAmps(nose_bend_angle_d);


%Calculate angular speed
%--------------------------------------------------------------------------
% Compute the speed centered between the back and front foraging movements.
%
%  1     2    3
%    d1    d2     d1 = 2 - 1,   d2 = 3 - 2
%        x        assign to x, avg of d1 and d2

%???? - why multiply and not divide by fps????

d_data = diff(nose_bend_angle_d) * fps;
speeds = NaN(size(amps));
speeds(2:end-1) = (d_data(1:(end - 1)) + d_data(2:end)) / 2;


%Propagate NaN for speeds to amps
%--------------------------------------------------------------------------
amps(isnan(speeds)) = NaN;

end

function amps = h__getAmps(nose_bend_angle_d)
%
%In between all sign changes, get the maximum or minimum value and
%apply to all indices that have the same sign within the stretch
%
%i.e. 
%
%1 2 3 2 1 -1 -2 -1 1 2 2 5 becomes
%3 3 3 3 3 -2 -2 -2 5 5 5 5
%
%   NOTE: This code is very similar to wormKinks

n_frames = length(nose_bend_angle_d);

dataSign      = sign(nose_bend_angle_d);
sign_change_I = find(dataSign(2:end) ~= dataSign(1:end-1));

end_I   = [sign_change_I; n_frames];
start_I = [1; sign_change_I+1];

%All Nan values are considered sign changes, remove these ...
mask = isnan(nose_bend_angle_d(start_I));
start_I(mask) = [];
end_I(mask)   = [];

amps = NaN(1,n_frames);

%For each chunk, get max or min, depending on whether the data is positive
%or negative ...
for iChunk = 1:length(start_I)
   cur_start = start_I(iChunk);
   cur_end   = end_I(iChunk);
   
   if nose_bend_angle_d(cur_start) > 0
       amps(cur_start:cur_end) = max(nose_bend_angle_d(cur_start:cur_end));
   else
       amps(cur_start:cur_end) = min(nose_bend_angle_d(cur_start:cur_end));
   end
end


end